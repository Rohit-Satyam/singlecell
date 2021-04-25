# My Single Cell Resources

### Video Lectures that helped me
1. Single-Cell RNA-Seq Analysis- Day 1 from [UCLA](https://www.youtube.com/watch?v=Cn5tI2oo1l0&t=10s).



### Installation of CellRanger for rawdata alignment

Follow along this [link](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) to install and setup cellranger pipeline.

#### 1. Convert .bcl files to fastq
The cellranger `mkfastq` pipeline is a 10x-enhanced wrapper around Illumina bcl2fastq, which demultiplexes BCL files from a sequencer into FASTQs for analysis. For example, I will download a small `.bcl` file to show the functionality:

```bash
wget http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
tar -xvzf cellranger-tiny-bcl-1.2.0.tar.gz
```
To run `mkfastq pipeline`, an Illumina Experiment Manager (IEM) sample sheet is also required. Note that his sample sheet is an example that is only valid for the `Single Cell 3′ v2 chemistry`.  enter command below to download the data sheet.

```bash
wget http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-samplesheet-1.2.0.csv
```
The samplesheet looks like:
```bash
Lane,Sample,Index
1,test_sample,SI-P03-C9
```

This file tells the pipeline that the library "test_sample" was sequenced on lane 1 of the flow cell and indexed using the SI-P03-C9 index set. The Sample column is used as the prefix for naming output files. This prefix also serves as the sample id in later Cell Ranger pipelines. In this case the prefix/sample id is test_sample.

```bash
## Install the bcl2fq first.This will be used to write fast files. Switch to this environment.

mamba create -n bcl2fq -c dranew bcl2fastq

## Also export the Path of cell ranger. To make it permanent, paste it in .bashrc and source
export PATH=/home/pichkari/rohit/tools/cellranger-6.0.1:$PATH
source ~/.bashrc

conda activate bcl2fastq
```

Running `mkfastq`

There are three  _arguments_  or inputs that are added to the  cellranger mkfastq  command:  `–-id`,  `--run`, and  `--csv`.

The  `--id`  can be anything. It is used by the pipeline to name the output directory that Cell Ranger is going to create to run in. This directory is called a  _pipestance_, which is short for pipeline instance.

The  `--run`  argument points to the Illumina run folder that contains the BCL files.

The  `--csv`  argument is a comma-separated values (CSV) file that describes how samples were indexed on the Illumina flow cell.

```bash
cellranger mkfastq --id=tutorial_walk_through \
--run=/home/pichkari/rohit/tools/cellranger-tiny-bcl-1.2.0 \
--csv=/home/pichkari/rohit/tools/cellranger-tiny-bcl-simple-1.2.0.csv
```
Run times vary based on the system resources, but it shouldn’t take more than a few minutes. Now go to the fastq_path directory.

```bash
H35KCBCXY  Reports  Stats  Undetermined_S0_L001_I1_001.fastq.gz  Undetermined_S0_L001_R1_001.fastq.gz  Undetermined_S0_L001_R2_001.fastq.gz
```
The Undetermined FASTQ files here at this level contain sequences that were unable to be assigned to valid index.

Demultiplexed FASTQ files with valid sequencing indices are found under the directory named after the flow cell id, in this case 'H35KCBCXY'.

![enter image description here](https://i.imgur.com/YuMcUUr.png)

>FastQ is the most raw form of scRNASeq data you will encounter. All scRNASeq protocols are sequenced with paired-end sequencing. Barcode sequences may occur in one or both reads depending on the protocol employed. However, protocols using unique molecular identifiers (UMIs) will generally contain one read with the cell and UMI barcodes plus adapters but without any transcript sequence. Thus reads will be mapped as if they are single-end sequenced despite actually being paired end.













### References
1. https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_fq
2. Obsolete tutorial still informative: https://bioinformatics.uconn.edu/single-cell-rna-sequencing-cell-ranger-2/#LibPrep
3. Broad Institute tutorial: https://broadinstitute.github.io/2019_scWorkshop/data-preprocessing.html
4. https://davetang.org/muse/2018/08/09/getting-started-with-cell-ranger/
