# My Single Cell Resources

### Video Lectures that helped me
1. Single-Cell RNA-Seq Analysis- Day 1 from [UCLA](https://www.youtube.com/watch?v=Cn5tI2oo1l0&t=10s).



### Installation of CellRanger for rawdata alignment

Follow along this [link](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) to install and setup cellranger pipeline.

#### Convert .bcl files to fastq
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

| \[Header\]        |                    |              |               |              |               |       |                 |             |
| ----------------- | ------------------ | ------------ | ------------- | ------------ | ------------- | ----- | --------------- | ----------- |
| IEMFileVersion    | 4                  |              |               |              |               |       |                 |             |
| Investigator Name | rjr                |              |               |              |               |       |                 |             |
| Experiment Name   | hiseq\_test        |              |               |              |               |       |                 |
| Date              | ########           |              |               |              |               |       |                 |             |
| Workflow          | GenerateFASTQ      |              |               |              |               |       |                 |
| Application       | HiSeq FASTQ Only   |              |               |              |               |       |                 |
| Assay             | TruSeq HT          |              |               |              |               |       |                 |
| Description       | hiseq sample sheet |              |               |              |               |       |                 |
| Chemistry         | Default            |              |               |              |               |       |                 |             |
|                   |                    |              |               |              |               |       |                 |             |
| \[Reads\]         |                    |              |               |              |               |       |                 |             |
| 26                |                    |              |               |              |               |       |                 |             |
| 98                |                    |              |               |              |               |       |                 |             |
|                   |                    |              |               |              |               |       |                 |             |
| \[Settings\]      |                    |              |               |              |               |       |                 |
|                   |                    |              |               |              |               |       |                 |             |
| \[Data\]          |                    |              |               |              |               |       |                 |             |
| Lane              | Sample\_ID         | Sample\_Name | Sample\_Plate | Sample\_Well | I7\_Index\_ID | index | Sample\_Project | Description |
| 1                 | s1                 | test\_sample |               | SI-P03-C9    | SI-P03-C9     | p1    |                 |
|                   |                    |              |               |              |               |       |                 |             |

