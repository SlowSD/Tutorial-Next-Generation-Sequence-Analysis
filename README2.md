# Tutorial: Next Generation Sequence Analysis - Stepwise

**This repository attempts to provide a stepwise tutorial on the Next Generation Sequence (NGS) analysis.**

In biological contexts, sequencing refers to determining the order of nucleotides (A, T, G, C) in a DNA sample. The term next generation sequencing (NGS) refers to an advanced and high-throughput methodnology of determining nucleotides order.

> **What about sequencing of nucleotides from RNA sample?**
>
> Similar to the DNA sequencing, RNA sequencing also does exist. While dealing with RNA samples, for most biological interpretations RNA is converted into an intermediate form known as cDNA (complementary DNA). However, methods namely direct RNA sequencing (DRS) involves sequencing of RNA directly.

## The goals of NGS analysis can be:

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

```{bash}
git clone https://github.com/SlowSD/Tutorial-Next-Generation-Sequence-Analysis.git
```

### 0.3. Creating a conda environment for this whole NGS analysis session and activating it

This step generates a conda environment in which tools required to perform the NGS analysis can be downloaded. In future, you would only need to activate the conda environment to resume the current analysis or re-run the analysis on different sequencing datasets.

```{bash}
conda deactivate
conda create --name NGS_analysis        #Creating new conda environment and naming it
conda activate NGS_analysis             #Activating the conda environment
```

#### 0.4.2. Add channels

```{bash}
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### 0.4.3. Installing required tools

This step would install required tools in the conda environment **NGS_analysis**

```{bash}
conda install -c bioconda bcftools bedtools blast bwa fastqc igv igvtools samtools sra-tools trim-galore vcftools trimmomatic bowtie2
```

---

## 1. Downloading sequencing data

As the complete sequence analysis is a multi-step and time-consuming procedure, we will be using a smaller dataset for the learning purposes. In this tutorial, we will be using example dataset featured by the [Qiagen](https://digitalinsights.qiagen.com/downloads/example-data/).\
The sequencing data comes from sequencing sample of *Pseudomonas aeruginosa* species.

💡

```{bash}
wget -L -O SRR396636.sra_1.fastq https://zenodo.org/records/11791175/files/SRR396636.sra_1.fastq?download=1
wget -L -O SRR396636.sra_2.fastq https://zenodo.org/records/11791175/files/SRR396636.sra_2.fastq?download=1
wget -L -O SRR396637.sra_1.fastq https://zenodo.org/records/11791175/files/SRR396637.sra_1.fastq?download=1
wget -L -O SRR396637.sra_2.fastq https://zenodo.org/records/11791175/files/SRR396637.sra_2.fastq?download=1
```

### 1. Saving input files as variables

```{bash}
read1_1 = "SRR396636.sra_1.fastq"
read1_2 = "SRR396636.sra_2.fastq"
read2_1 = "SRR396637.sra_1.fastq"
read2_2 = "SRR396637.sra_1.fastq"
```

## 2. Information on .fastq files

Once the sequencing data has been downloaded, with the listing command the 4 `.fastq` files can be seen in the current directory.

```{bash}
ls
```

<Insert picture>

### 2.1. Viewing .fastq files content

To view content of the `.fastq` file, there are multiple ways

#### 2.1.1. Using `cat` command

```{bash}
cat $read1_1 | head       #To print top 5 lines in the terminal
```

#### 2.2.2. Using `less` command

```{bash}
less $read1_1 | head       #To view top 5 lines in the terminal (without printing)
```

### 2.2 Determining the number of reads

To determine the number of reads in a `fastq` file, use the below command.

```{bash}
cat $read1_1 | grep @SRR | wc -l
```

> Here,\
>`grep` command works as search tool and searches for occurence of `SRR` in the `.fastq` file.\
> `wc -l` command count the occurences of `SRR`.

🔍 **Check yourself whether the mate file (i,e, SRR396636.sra_2.fastq) has the same no. of reads.**

### 2.3. Fastqc

Fastqc is a tool developed by **Babraham Bioinformatics** for looking at quality of fastq reads.

```{bash}
mkdir -p fastqc_results        #Creates a folder to save the fastqc results
fastqc *.fastq -o /fastqc_results           #Here, *.fastq acts as a wildcard for us.
```

<details>
    <summary><b>Reading fastqc report:</b></summary>
    <aside>

        > - Per base sequence quality:
        > - Per base sequence content:
        > - Per base GC content:
        > - Per base N content:
        > - Sequence Length distribution:
        > - Sequence Duplication levels:
        > - Overrepresented sequences:
        > - Adapter content:

</aside>
</details>

---

### 2.4 multiqc

While `fastqc` generates report for each .fastq files, `multiqc` compiles the results of fastqc reports as a single report. For working with multiqc, either the individual `fastqc reports` or the directory containing `fastqc results` could be provided. In the latter case, multiqc automatically searches for reports for its working.

```{bash}
mkdir -p multiqc_results       #Creates a folder to save multiqc_results
multiqc /fastqc_results/ . -o multiqc_results
```

## 3. Trimming reads

This step of trimming low quality reads and removal of adapters could be done by any of these three tools.

<details>
<summary>Using <b>Trimmomatic</b></summary>

> [Important resource](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

```{bash}
trimmomatic PE -phred33 \
 $read1_1 $read1_2 \
 SRR396636.sra_1_trim.fastq SRR396636.sra_2_trim.fastq \
 ILLUMINACLIP:adapter.fasta TRAILING:3 LEADING:3 MINLEN:36 SLIDINGWINDOW:4:15
```

</details>

<details>
<summary>Using <b>Cutadapt</b></summary>

```{bash}
cutadapt -a <adapter_seq_forw> -A <adapter_seq_rev> $read1_1 $read1_2 SRR396636.sra_1_trim.fastq SRR396636.sra_2_trim.fastq
```

</details>

<details>
<summary>Using <b>Trim Galore</b></summary>

```{bash}
trim_galore --paired --phred33 --gzip $read1_1 $read1_2 
```

</details>

---

## 4. Alignment of reads with the reference genome

### 4.1. Downloading reference genome

There are multiple methods to download reference genome.

<details>
<summary><b>Method 1: Using FTP (Best practice)</b></summary>
    <aside>

1. Go to NCBI homepage
2. Select assembly from the drop-down menu
3. Search for `Pseudomonas aeruginosa`
4. From the filters column at left, select `Reference`
5. Click on `ASM676v1` assembly

    ![alt text](image.png)
6. Select `FTP` to access the directory containig records on `Pseudomonas aeruginosa`.

    ![alt text](image-1.png)
7. Right click on `GCF_000006765.1_ASM676v1_genomic.fna.gz` and copy the link.

```{bash}
mkdir -P pseudomonas_aeruginosa_rg     #Creating directory to save reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz
```

</aside>
</details>

<details>
<summary><b>Method 2: Downloading reference genome from NCBI genomes.</b></summary>
<aside>

1. Go to [NCBI genome collection](https://ncbi.nlm.nih.gov/datasets/genome/)
2. Search `Pseudomonas aeurignosa`.
3. Select the reference genome assembly
4. Now you have three methods of downloading reference genome data.

    <details>
        <summary>Using <b>Download</b></summary>
   <aside>

   1. Click on `Download` button.
   2. Select file source as `Refseq only` and file types as `Genome sequence (FASTA)`.
   3. This will download a zip folder that needs to be unzipped for further use.
   </aside>
   </details>

    ---

    <details>
        <summary>Using <b>datasets</b> </summary>
    <aside>

    1. Installing `datasets` conda package

        ```{bash}
        conda install -c conda-forge ncbi-datasets-cli
        ```

    2. Click on `datasets`
    3. Copy the command

        ```{bash}
        datasets download genome accession GCF_000195955.2 --include genome
        ```

    4. Paste it in command line to download refernce genome and some additonal files.
    </aside>
    </details>

    ---

   <details>
       <summary>Using <b>FTP (Best practice)</b></summary>
   <aside>

   1. click on `FTP`.
   2. Right click on the file `GCF_000006765.1_ASM676v1_genomic.fna.gz` and copy the link.
   3. Use the `wget` or `curl` command.

        ```{bash}
        mkdir reference_genome          #Create a directory to save reference genome file.
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz          #In place of wget, curl command can also be used
        ```

   </aside>
   </details>

    ---

</details>

---

### 4.2. Indexing reference genome: one-time step

**Indexing a reference genome is a one-time step and can be performed using any of below tools. However, keep in mind that the tool used for generating reference index should also be used for alignment.**

<details>
    <summary>Using <b>samtools</b></summary>
<aside>

```{bash}
samtools faidx pseudomonas_aeruginosa_rg/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz
```

</aside>
</details>

<details>
    <summary>Using <b>bwa</b></summary>
<aside>

```{bash}
bwa index <ref_genome>
```

> This command generates the following 5 files:
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
    <summary>Using <b>bowtie2</b></summary>
<aside>

```{bash}
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

---

### 4.3. Performing alignment generating SAM file

<details>
    <summary>Using <b>samtools</b> </summary>
<aside>

```{bash}
samtools 
```

</aside>
</details>

<details>
    <summary>Using <b>bwa</b></summary>
<aside>

```{bash}
bwa mem <ref_genome> <file1_1.fq> <file1_2.fq> > <read1.sam>
```

</aside>
</details>

<details>
    <summary>Using <b>bowtie2</b></summary>
<aside>

```{bash}
bowtie2 -x <prefix> -1 <file1_1.fq> -2 <file1_2.fq> -S <read1.sam>
```

</aside>
</details>

---

## 5 Analysis of alignment results

### 5.1 Converting SAM file into BAM file

Using `samtools`

```{bash}
samtools view -Sb <file.sam> > -o <file1.bam>
```

> -Sb : Input in SAM format (S) and the output will be be BAM format(b)

### 5.2 Sorting BAM file

Using `samtools`

```{bash}
samtools sort <file1.bam> -o <file1_sorted.bam>
```

---

### Samtools fixmate & Samtools markdup

Sorting reads by `read name` rather than `chromosomal coordinates`

```{bash}
samtools sort -n <file1.bam> -o <file1_namesort.bam>

```

Appending 'ms' and 'MC' tags

```{bash}
samtools fixmate -m <file1_namesort.bam> <file1_namesort_fixmate.bam>

```

Re-sort the BAM files by `chromosomal coordinates`

```{bash}
samtools sort <file1_namesort_fixmate.bam> -o <file1_namesort_fixmate_sort.bam>

```

Marking the duplicates.

```{bash}
samtools markdup -r <file1_namesort_fixmate_sort.bam> <file1_namesort_fixmate_sort_markdup.bam>

```

---

### 5.6 Generating BAM file statistics

Using `samtools`

```{bash}
samtools flagstat <file1_sorted.bam> -o <file1_sorted.txt>
```

### 5.7 Viewing SAM or BAM file

Using `samtools`

```{bash}
samtools view <file1_sorted.sam/bam>
```

### 4.8 Creating a BAM file index to view in IGV

Using `samtools`

```{bash}
samtools index <file1_sorted.bam>
```

## 4.9 Viewing whole alignment

Using #samtools

```{bash}
samtools tview <file1_sorted.bam> <ref_genome>

```

---

## Qualimap

```{bash}
qualimap bamqc -outdir <dir> \\
-bam <file1.bam> \\
-gff <file.gff>

```

---

## 5. Generating VCF file for variant calling

- Estimate of variant frequency
- Measure of confidence

Using #freebayes

```{bash}
freebayes -f <ref_genome> <file1.bam> > <file1.vcf>

```

Using #bcftools
**Step1: Generating pileup file**
Generates genotype likelihoods at each genomic position with coverage.

```{bash}
bcftools mpileup -O b -f <ref_genome> <file1_sorted.bam> -o <file1_raw.bcf

```

**Step2: Detecting SNVs**

```{bash}
bcftools call -v <file1_raw.bcf> -o <file1_variants.vcf>
```

> -v: outputs variant sites only.
>

*In a single step*

```{bash}
bcftools mpileup -Ou -f <ref_genome> <file1.bam> | bcftools call -mv -Ob -o <file1.bcf>

```

> -O: output type
b: compressed BCf, u: uncompressed BCf
z: compressed VCF, v: uncompressed VCF
>

**Determining the number of variant sites**

```{bash}
grep -v -c "^#" <file1.vcf>

```

---

## 6. Viewing VCF file

Using `bcftools`

```{bash}
bcftools view <file1.vcf>

```

### 6.1 Generating an index file for VCF to view in IGV

Useful commands
*Extracting information from VCF file (can be saved as a #bedfile), ```

```{bash}
bcftools query -f '%COLa\\t%COLb\\t%COLc\\n' <file1.vcf>

```

*Converting BCF file to VCF file*

```{bash}
bcftools view <file1.bcf> > <file1.vcf>

```

*Sorting VCF file by chromosome → coordinate*

```{bash}
vcf-sort <file1.vcf> > <file1.sorted.vcf>

```

---

*Converting VCF to BED file*

Using #vcf2bed**

```{bash}
vcf2bed < <file1.vcf> > <file1.bed>

```

*identify annotated genes which are not covered by reads across their full length.*

```{bash}
bedtools coverage -a <ref_genome.gff> -b <file.bam> -o <file1_gene-coverage.txt>

```

[[Generating an assembly with unmapped reads]]

### Some important concepts

#### Sequencing depth
