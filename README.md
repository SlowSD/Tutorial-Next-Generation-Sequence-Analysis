# Tutorial: Next Generation Sequence Analysis - Stepwise

**This repository attempts to provide a stepwise tutorial on the Next Generation Sequence (NGS) analysis.**

In biological contexts, sequencing refers to determining the order of nucleotides (A, T, G, C) in a DNA sample. The term next generation sequencing (NGS) refers to an advanced and high-throughput methodnology of determining nucleotides order.

> **What about sequencing of nucleotides from RNA sample?**
>
> Similar to the DNA sequencing, RNA sequencing also does exist. While dealing with RNA samples, for most biological interpretations RNA is converted into an intermediate form known as cDNA (complementary DNA). However, methods namely direct RNA sequencing (DRS) involves sequencing of RNA directly.

### The goal of NGS analyses can be

- Determining genes expression (a.k.a single cell RNA sequencing). Here, the complete mRNA pool of a cell is first converted into cDNAs using reverse transcriptase. Then cDNAs are sequenced to determine their counts. (cDNA counts directly correlates with gene expression levels).
- Identifying nucleotide mutations or variants.
- Generating a reference genome of a newly identified species. Upon encountering a novel species, NGS can be utilized to determine the complete genome sequence.
- In case of diseases, after symptoms based characterization, NGS can be utilized to identify and confirm the disease causing species or a species subtype.

---

## 0. Prerequisite for this tutorial

In this section, we will prepare to begin the analysis.

### 0.1. **Linux architecture and knowledge of unix commands**

To follow this tutorial, you will need a Linux terminal.

- If you are on mac, you only need to open command prompt.
- If you are on windows, you need to activate **Windows subsystem for linux (WSL)**. To activate WSL, you can follow this [youtube tutorial](https://www.youtube.com/watch?v=FiGLvQTsfC8).

>💡 If don’t have any prior experience with unix commands, I would suggest you to check out this article on [Basic Linux commands](https://www.notion.so/Basic-Linux-commands-680d073407834ad981dab47f9411eb08?pvs=21). However, for this analysis, you would only need to know the `ls`, `grep`, `cd`, `mkdir` commands.

### 0.2. **Follow along this tutorial**

To follow this tutorial, the best is to clone this repository in your local device.

```
git clone <>
```

### 0.2. Creating a working directory to save the analysis progress

```
mkdir ngs_analysis && cd ngs_analysis
```

### 0.3. Creating a conda environment for this whole NGS analysis session and activating it

This step generates an environment in which tools required to perform NGS analysis can be downloaded. Therefore, in future whenever you would be performing NGS analysis you will need to simply activate the environment and begin analysis

```
conda create --name NGS_analysis
conda activate NGS_analysis
```

#### 0.4.2. Add channels

```
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### 0.4.3. Installing required tools

```
conda install -c bioconda bcftools bedtools blast bwa fastqc igv igvtools samtools sra-tools trim-galore vcftools
```

---

## 1. Downloading sequencing data

As sequence analysis is a multi-step procedure, we will be using a smaller dataset for learning purposes. I am using example data featured by [Qiagen](https://digitalinsights.qiagen.com/downloads/example-data/).

To download this dataset, copy and paste each links in the url.

```
https://zenodo.org/records/11791175/files/SRR396636.sra_1.fastq
https://zenodo.org/records/11791175/files/SRR396636.sra_2.fastq
https://zenodo.org/records/11791175/files/SRR396637.sra_1.fastq
https://zenodo.org/records/11791175/files/SRR396637.sra_2.fastq
```

---

## 2. Information on .fastq files

Once the sequencing data has been downloaded, with the listing command the files can be seen in the current directory.

```
ls
```

### 2.1. Viewing .fastq files content

```
zcat SRR28409626_1.fastq.gz | head -n 8 # To view top 8 lines or top 2 reads
```

```
zcat SRR28409626_2.fastq.gz | head -n 8 # To view top 8 lines or top 2 reads
```

You should see that both files have same reads ID, i.e `SRR28409626_1` and `SRR28409626_2`

### 2.2 Determining the number of reads

```
zcat SRR28409626_1.fastq.gz | grep @SRR | wc -l
```

Check yourself whether the other file have the same no. of reads.

### 2.3. fastqc

During the sequecning process due to many possible reasons, the sequence process could be prone to errors. This would lead to

As the sequencing process is a time consuming process,

```
mkdir fastqc_results
fastqc SRR28409626_1.fastq.gz SRR28409626_2.fastq.gz -o /fastqc_results
```

> Reading fastqc report:
    > - Per base sequence quality:
    > - Per base sequence content
    > - Per base GC content
    > - Per base N content
    > - Sequence Length distribution
    > - Sequence Duplication levels
    > - Overrepresented sequences
    > - Adapter content

### 2.4 multiqc

While `fastqc` generates report for each .fastq files, `multiqc` takes fastqc reports as an input and generates one report.

```
mkdir multiqc_results
multiqc *.fastqc -o multiqc_results
```

---

## 3. Trimming reads

This step could be done by any of these three methods.

> Is it necessary?
Removal of low quality reads.
Reads containing adaptor sequences.

<details>
<summary><h3>Using Trimmomatic</h3></summary>

> [Important resource](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

```
java -jar /trimmomatic-0.39.jar\\
 PE -phred33\\
 <input_file1_1> <input_file1_2>\\
 <output_file1_1> <output_file1_2>\\
 ILLUMINACLIP:adapter.fasta TRAILING:3 LEADING:3 MINLEN:36 SLIDINGWINDOW:4:15
```

</details>

<details>
<summary><h3>Using Cutadapt</h3></summary>

```
cutadapt -a <adapter_seq_forw> -A <adapter_seq_rev> <file1_1.fastq.gz> <file1_2.fastq.gz> <file1_1.trimmed.fastq.gz> <file1_2.trimmed.fastq.gz>
```

</details>

<details>
<summary><h3>Using Trim Galore</h3></summary>

```
trim_galore --paired --phred33 --gzip <file1_1.fastq.gz> <file1_2.fastq.gz> 
```

</details>

---

## 4. Alignment of reads with the reference genome

### 4.1. Downloading reference genome

The information on *Mycobacterium tuberculosis* reference genome is [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000195955.2/).

**If you simply want to download the reference genome and would liek to move on to the next step of the analysis, use this code.**

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
```

However, if you want to learn the steps to find a reference genome assembly and download it, here below are the three methods.

<details>
    <summary><b>Method 1: With no use of command</b></summary>
<aside>

1. Go to NCBI homepage
2. Select assembly from the drop-down menu
3. Search for `Mycobacterium tuberculosis`
4. Click on `ASM19595v2` assembly
5. Click on `Download assembly`
6. Select `Refseq` and `Genomic Fasta (.fna)`

    ![image](https://github.com/SlowSD/Tutorial-Next-Generation-Sequence-Analysis/assets/111181145/f7eb9ecc-bbd5-42f1-9501-8b79e89d0191)

</aside>
</details>

<details>
    <summary><b>Method 2: Best practice</b></summary>
<aside>

1. Go to NCBI homepage
2. Select assembly from the drop-down menu
3. Search for `Mycobacterium tuberculosis`
4. Click on `ASM19595v2` assemblyRepeat the 1-4 steps of Method 1.
5. Click on '**FTP directory for RefSeq assembly**'

    ![image](https://github.com/SlowSD/Tutorial-Next-Generation-Sequence-Analysis/assets/111181145/4982efb3-9fa0-46a9-a923-348c17fbc958)
4. Right click on `GCF_000195955.2_ASM19595v2_genomic.fna.gz` and copy link.
5. Download the refernece sequence data `wget` command.

    ```
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
    ```

</aside>
</details>

<details>
<summary><b>Method 3: Downloading reference genome from NCBI genomes.</b></summary>
<aside>

1. Go to [NCBI genome collection](https://ncbi.nlm.nih.gov/datasets/genome/)
2. Select/search `Mycobacterium tuberculosis`
3. Select the reference genome assembly
4. Now you have three ways of downloading reference genome data.

   ![image](https://github.com/SlowSD/Tutorial-Next-Generation-Sequence-Analysis/assets/111181145/9dc3b1da-8a9a-4856-9fc4-63ec7ea34c7f)

    <details>
        <summary><b> Using Download</b></summary>
   <aside>

   1. Click on download button
   2. Select `Genome sequences (FASTA)`
   </aside>
   </details>

    <details>
        <summary><b> Using dataset command </b></summary>
   <aside>

    1. Installing datasets conda package

        ```
        conda install -c conda-forge ncbi-datasets-cli
        ```

    Now `datasets` function is ready to execute commands.

    2. Click on `datasets`
    3. Copy the command

        ```
        datasets download genome accession GCF_000195955.2 --include gff3,rna,cds,protein,genome,seq-report
        ```

    4. Paste it in command line to download refernce genome and some additonal files.
   </aside>
   </details>

   <details>
       <summary><b>Using curl command</b></summary>
   <aside>

   1. click on `curl`
   2. Copy the command
   3. Download the dataset with `curl` command.

       ```
       curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000195955.2/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED
       ```

   </aside>
   </details>

</details>

### 4.2. Indexing reference genome: one-time step

**This is a one-time step and can be performed using any of these tools.**

<details>
    <summary><h4>Using samtools</h4></summary>
<aside>

```
samtools faidx <ref_genome>
```

</aside>
</details>

<details>
    <summary><h4>Using bwa</h4></summary>
<aside>

```
bwa index <ref_genome>
```

> This command generates:
>
> .amb
>
> .ann
>
> .bwt: Binary file
>
> .pac: Binary file
>
> .sa: Binary file

</aside>
</details>

<details>
    <summary><h4>Using bowtie2</h4></summary>
<aside>

```
bowtie2-build <ref_genome> <prefix>
```

> This command generates 6 files
>
> <prefix>.1.bt2
>
> <prefix>.2.bt2
>
> <prefix>.3.bt2
>
> <prefix>.4.bt2
>
> <prefix>.rev.1.bt2
>
> <prefix>.rev.2.bt2
  
</aside>
</details>

### 4.3. Performing alignment generating SAM file

<details>
    <summary><h4>Using samtools</h4> </summary>
<aside>

```
samtools 
```

</aside>
</details>

<details>
    <summary><h4>Using bwa</h4></summary>
<aside>

```
bwa mem <ref_genome> <file1_1.fq> <file1_2.fq> > <read1.sam>
```

</aside>
</details>

<details>
    <summary><h4>Using bowtie2</h4></summary>
<aside>

```
bowtie2 -x <prefix> -1 <file1_1.fq> -2 <file1_2.fq> -S <read1.sam>
```

</aside>
</details>

### 4.4 Converting SAM into BAM file

Using #samtools

```
samtools view -Sb <file.sam> > -o <file1.bam>

```

> -Sb : Input in SAM format (S) and the output will be be BAM format(b)
>

### 4.5 Sorting BAM file

Using #samtools

```
samtools sort <file1.bam> -o <file1_sorted.bam>

```

---

### Samtools fixmate & Samtools markdup

Sorting reads by `read name` rather than `chromosomal coordinates`

```
samtools sort -n <file1.bam> -o <file1_namesort.bam>

```

Appending 'ms' and 'MC' tags

```
samtools fixmate -m <file1_namesort.bam> <file1_namesort_fixmate.bam>

```

Re-sort the BAM files by `chromosomal coordinates`

```
samtools sort <file1_namesort_fixmate.bam> -o <file1_namesort_fixmate_sort.bam>

```

Marking the duplicates.

```
samtools markdup -r <file1_namesort_fixmate_sort.bam> <file1_namesort_fixmate_sort_markdup.bam>

```

---

### 4.6 BAM file statistics

Using #samtools

```
samtools flagstat <file1_sorted.bam> -o <file1_sorted.txt>

```

### 4.7 Viewing SAM or BAM file

Using #samtools

```
samtools view <file1_sorted.sam/bam>

```

### 4.8 Creating a BAM file index to view in IGV

Using #samtools

```
samtools index <file1_sorted.bam>

```

## 4.9 Viewing whole alignment

Using #samtools

```
samtools tview <file1_sorted.bam> <ref_genome>

```

---

## Qualimap

```
qualimap bamqc -outdir <dir> \\
-bam <file1.bam> \\
-gff <file.gff>

```

---

## 5. Generating VCF file for variant calling

- Estimate of variant frequency
- Measure of confidence

Using #freebayes

```
freebayes -f <ref_genome> <file1.bam> > <file1.vcf>

```

Using #bcftools
**Step1: Generating pileup file**
Generates genotype likelihoods at each genomic position with coverage.

```
bcftools mpileup -O b -f <ref_genome> <file1_sorted.bam> -o <file1_raw.bcf

```

**Step2: Detecting SNVs**

```
bcftools call -v <file1_raw.bcf> -o <file1_variants.vcf>

```

> -v: outputs variant sites only.
>

*In a single step*

```
bcftools mpileup -Ou -f <ref_genome> <file1.bam> | bcftools call -mv -Ob -o <file1.bcf>

```

> -O: output type
b: compressed BCf, u: uncompressed BCf
z: compressed VCF, v: uncompressed VCF
>

**Determining the number of variant sites**

```
grep -v -c "^#" <file1.vcf>

```

---

## 6. Viewing VCF file

```
bcftools view <file1.vcf>

```

### 6.1 Generating an index file for VCF to view in IGV

Useful commands
*Extracting information from VCF file (can be saved as a #bedfile), ```

```
bcftools query -f '%COLa\\t%COLb\\t%COLc\\n' <file1.vcf>

```

*Converting BCF file to VCF file*

```
bcftools view <file1.bcf> > <file1.vcf>

```

*Sorting VCF file by chromosome → coordinate*

```
vcf-sort <file1.vcf> > <file1.sorted.vcf>

```

---

*Converting VCF to BED file*

Using #vcf2bed**

```
vcf2bed < <file1.vcf> > <file1.bed>

```

*identify annotated genes which are not covered by reads across their full length.*

```
bedtools coverage -a <ref_genome.gff> -b <file.bam> -o <file1_gene-coverage.txt>

```

[[Generating an assembly with unmapped reads]]

### Some important concepts

#### Sequencing depth
