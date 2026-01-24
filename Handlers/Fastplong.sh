#!/bin/env bash

#   Handler for long-read quality filtering using fastplong
#   For ONT and PacBio reads
#   SUPPORTS RESUME functionality via checkpoint logging

set -o pipefail

#   shellcheck disable=SC2034
declare -a Fastplong_Handler_Dependencies=(fastplong parallel)

function runFastplong() {
    local sampleName="$1"
    local input="$2"
    local outdir="$3"
    local adapters="$4"
    local skipAdapter="$5"
    local out="${outdir}/${sampleName}"

    mkdir -p "${out}"

    local -a adapterArgs=()
    if [[ "${skipAdapter}" == "true" ]]; then
        adapterArgs=(--disable_adapter_trimming)
    else
        if [[ -n "${adapters}" ]]; then
            adapterArgs=(--adapter_fasta "${adapters}")
        fi
    fi

    echo "  [RUNNING] fastplong for: ${sampleName}..."

    # Run fastplong and capture exit code
    if fastplong \
        -i "${input}" \
        -o "${out}/${sampleName}_trimmed.fastq.gz" \
        --thread 4 \
        "${adapterArgs[@]}" \
        --html "${out}/${sampleName}_fastp.html" \
        --json "${out}/${sampleName}_fastp.json" \
        --report_title "${sampleName} fastp report"; then
        return 0
    else
        local exit_code=$?
        echo "ERROR: fastplong failed for ${sampleName} with exit code ${exit_code}" >&2
        return ${exit_code}
    fi
}

export -f runFastplong

function Fastplong() {
    local sampleList="$1"
    local outPrefix="$2"/Fastplong
    local project="$3"
    local adapters="$4"
    local configPath="${5:-${CONFIG_FASTP:-}}"

    echo "Starting Fastplong Handler..." >&2

    if [[ ! -f "$sampleList" ]]; then
        echo "ERROR: sampleList file not found: $sampleList" >&2
        return 1
    fi

    mkdir -p "$outPrefix"

    # --- CHECKPOINT FILES ---
    # 1. Manifest: The full paths of resulting files (for next steps)
    local manifestFile="${outPrefix}/${project}_fastplong_trimmed.txt"
    # 2. Checkpoint Log: List of Sample Names that have finished successfully
    local checkpointFile="${outPrefix}/${project}_checkpoint.log"

    # --- RESUME LOGIC ---
    # Load completed samples into an associative array for fast lookup
    declare -A completedSamples
    if [[ -f "$checkpointFile" ]]; then
        echo "Found checkpoint file. Loading completed samples..."
        while read -r doneSample; do
            [[ -z "$doneSample" ]] && continue
            completedSamples["$doneSample"]=1
        done < "$checkpointFile"
        echo "Resuming: ${#completedSamples[@]} samples already finished."
    else
        # If no checkpoint, clear the manifest to avoid duplicates from a partial previous run
        : > "${manifestFile}"
        : > "${checkpointFile}"
    fi

    # --- ADAPTER LOGIC ---
    local skipAdapter="false"
    if [[ "${FORCE_ADAPTER_TRIM:-}" == "true" ]]; then
        skipAdapter="false"
        echo "Adapter trimming FORCED ON via FORCE_ADAPTER_TRIM" >&2
    elif [[ "${SKIP_ADAPTER_TRIM:-}" == "true" ]]; then
        skipAdapter="true"
        echo "Adapter trimming FORCED OFF via SKIP_ADAPTER_TRIM" >&2
    elif [[ -n "$configPath" ]] && [[ -f "$configPath" ]]; then
        local seqPlatform
        local minimapType
        seqPlatform=$(grep -E "^SEQ_PLATFORM=" "$configPath" | tail -n1 | cut -d= -f2 | tr '[:lower:]' '[:upper:]')
        minimapType=$(grep -E "^MINIMAP2_READ_TYPE=" "$configPath" | tail -n1 | cut -d= -f2 | tr '[:lower:]' '[:upper:]')
        if [[ "$seqPlatform" == "PACBIO" ]] && [[ "$minimapType" =~ HIFI ]]; then
            skipAdapter="true"
            echo "Auto-detected PacBio HiFi - skipping adapter trimming" >&2
        fi
    fi

    echo "skipAdapter = $skipAdapter" >&2

    # --- PROCESS FUNCTION ---
    function process_sample_direct() {
        local sampleName="$1"
        shift
        local files=("$@")

        [[ "+#files[@]" -eq 0 ]] && return

        # 1. CHECKPOINT CHECK
        if [[ -n "${completedSamples[$sampleName]}" ]]; then
            echo "  [SKIP] ${sampleName} is already done."
            return 0
        fi

        local sampleDir="${outPrefix}/${sampleName}"
        local tempConcat="${sampleDir}/${sampleName}_concat_temp.fastq.gz"
        local expectedOutput="${sampleDir}/${sampleName}_trimmed.fastq.gz"

        mkdir -p "${sampleDir}"

        echo "  [CONCAT] Concatenating ${#files[@]} files for ${sampleName}..."

        # Concatenate all input files
        if ! cat "${files[@]}" > "${tempConcat}"; then
            echo "ERROR: Failed to concatenate files for ${sampleName}" >&2
            rm -f "${tempConcat}"
            return 1
        fi

        # Verify concatenated file exists and is not empty
        if [[ ! -s "${tempConcat}" ]]; then
            echo "ERROR: Concatenated file is empty or missing for ${sampleName}" >&2
            rm -f "${tempConcat}"
            return 1
        fi

        # Run fastplong on the concatenated file
        if runFastplong "${sampleName}" "${tempConcat}" "${outPrefix}" "${adapters}" "${skipAdapter}"; then
            local exit_code=0
        else
            local exit_code=$?
        fi

        # Clean up temporary concatenated file
        rm -f "${tempConcat}"

        # 2. SUCCESS HANDLING
        if [[ $exit_code -eq 0 ]]; then
            # Verify output file was created and is not empty
            if [[ ! -s "${expectedOutput}" ]]; then
                echo "ERROR: Expected output file not created or is empty: ${expectedOutput}" >&2
                return 1
            fi

            echo "  [COMPLETE] ${sampleName}"
            
            # Log to Checkpoint (Sample Name) - atomic append
            echo "${sampleName}" >> "${checkpointFile}"
        else
            echo "ERROR: fastplong failed for ${sampleName} with exit code ${exit_code}" >&2
            return ${exit_code}
        fi
    }

    # Export the function and necessary variables for potential parallel execution
    export -f process_sample_direct
    export outPrefix
    export adapters
    export skipAdapter

    # --- PARSING LOOP ---
    local currentSample=""
    local -a sampleFiles=()
    local -a failedSamples=()
    
    while read -r line; do
        [[ -z "$line" ]] && continue

        if [[ "$line" =~ \\.fastq|fq\\.gz$ ]]; then
            sampleFiles+=("$line")
        else
            if [[ -n "$currentSample" ]] && [[ "+#sampleFiles[@]" -gt 0 ]]; then
                if ! process_sample_direct "$currentSample" "${sampleFiles[@]}"; then
                    failedSamples+=("$currentSample")
                fi
            fi
            currentSample="$line"
            sampleFiles=()
        fi
    done < "${sampleList}"
    
    # Process the last sample
    if [[ -n "$currentSample" ]] && [[ "+#sampleFiles[@]" -gt 0 ]]; then
        if ! process_sample_direct "$currentSample" "${sampleFiles[@]}"; then
            failedSamples+=("$currentSample")
        fi
    fi
    
    # --- REGENERATE MANIFEST FROM CHECKPOINT ---
    # This ensures the manifest is always accurate, even after resume
    echo "Regenerating manifest from checkpoint..."
    : > "${manifestFile}"
    
    while read -r doneSample; do
        [[ -z "$doneSample" ]] && continue
        local trimmedPath="${outPrefix}/${doneSample}/${doneSample}_trimmed.fastq.gz"
        if [[ -f "${trimmedPath}" ]]; then
            echo "${trimmedPath}" >> "${manifestFile}"
        else
            echo "WARNING: Checkpoint lists ${doneSample} as done, but output file missing: ${trimmedPath}" >&2
        fi
    done < "$checkpointFile"

    # Sort the manifest for consistency
    sort -o "${manifestFile}" "${manifestFile}"
    
    echo ""
    echo "============================================"
    echo "Fastplong Handler Complete"
    echo "============================================"
    echo "Manifest:    ${manifestFile}"
    echo "Checkpoint:  ${checkpointFile}"
    echo "Total completed: $(wc -l < "${checkpointFile}") samples"
    
    if [[ "+#failedSamples[@]" -gt 0 ]]; then
        echo ""
        echo "WARNING: ${#failedSamples[@]} sample(s) failed:"
        printf '  - %s\n' "${failedSamples[@]}"
        echo ""
        echo "To retry failed samples, run this handler again."
        return 1
    fi
    
    echo "============================================"
    return 0
}

export -f Fastplong
