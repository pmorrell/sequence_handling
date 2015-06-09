# sequence_handling
### A series of scripts to automate genome sequence workflows
___
___
## installer.sh
The `installer.sh` script installs [_Seqqs_](https://github.com/morrelllab.seqqs), [_Sickle_](https://github.com/vsbuffalo/sickle), and [_Scythe_](https://github.com/vsbuffalo/scythe) for use with the `Quality_Triming.sh` script. All three of these are required for use with the `Quality_Triming.sh` as-is. Running `installer.sh install` runs the installation process, temporarily adding each program to the `PATH`. To permanently add these programs to the `PATH` for future use, please follow on-screen instructions or run `installer.sh echo`. `installer.sh` is designed to be run from the command line directly.

## sample\_list_generator.sh
The `sample_list_generator.sh` script creates a list of samples using a directory tree for its searching. This will find **all** samples in a given directory and its subdirectories. Only use this if you are using all samples within a directory tree. Running it with no arguments will give a detailed usage message, or one can edit the script to have variables hard-coded. `sample_list_generator.sh` is designed to be run from the command line directly.
___
## read_counts.sh
The `read_counts.sh` script calls _bioawk_ to get accurate counts for read number for a list of samples. The sample list is currently hard-coded into the script to permit qsub job submission. Output is written toa tab-delimited file file with sample name drawn from the file name for the list of samples.

_Bioawk_ is available through [Github](https://github.com/lh3/bioawk).

## FastQC.sh
The `FastQC.sh` script runs [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on the command line on a series of samples organized in a project directory for quality control. Our recommendation is using this both before and after quality trimming and before read mapping. This script is designed to be run using the Portable Batch System

## Quality\_Trimming.sh
The `Quality_Trimming.sh` script runs `trim_autoplot.sh` (part of the [_Seqqs_](https://github.com/morrelllab.seqqs) repository on GitHub) on a series of samples organized in a project directory.. In addition to requiring _Seqqs_ to be installed, this also requires [GNU Parallel](http://www.gnu.org/software/parallel/) to be installed on the system. This script is set up to be run using the [Portable Batch System](http://www.pbsworks.com/).

## Read\_Mapping.sh
The `Read_Mapping.sh` scripts uses the [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/) (BWA) to read map a series of sequences. These scripts are both designed to use the results from `Quality_Trimming.sh` for the read mapping. Both scripts find files organized into a project directory and automatically sort them by sample. `QSub_Read_Mapping_Parallel.sh` uses [GNU Parallel](http://www.gnu.org/software/parallel/) to run the read mapping for each sample in parallel. This script designed to be run using the Portable Batch System.
