# Tutorial-Next-Generation-Sequence-Analysis

This repository attempts to provide a easy tutorial on the Next Generation Sequence (NGS) analysis.
What happens in a NGS analysis?

A next generation sequencing (NGS) involves the sequencing of genome of a species. The principal goal of NGS analysis is often to look and characterize the nucleotide variations or mutations in the sequenced genome.

Why performing NGS analysis is important?

A NGS analysis 

could provide many insights such as

- List of mutations in the genome.
- Identifying mutation
    - Location
    - Type

---

## 0. Prerequisite for this tutorial

1. **Linux architecture and knowledge of bash commands.**
    
    <aside>
    💡 If don’t have any prior experience with bash commands, I would suggest you to check out this article.
    
    [Basic Linux commands](https://www.notion.so/Basic-Linux-commands-680d073407834ad981dab47f9411eb08?pvs=21)
    
    </aside>
    
2. **Sequencing data to analyze.**
    
    To perform NGS analysis, we need sequencing reads to input. The largest public repository that stores the raw sequencing data is the [**Sequence Read Archive (SRA)**](https://www.ncbi.nlm.nih.gov/sra).
    
3. **Modules required to perform analysis.**
    - SRA toolkit
    - 
    

---

## 1. Downloading short reads from SRA

The whole genome sequencing process is time taking and involves reading each base position. Thus, the larger the species genome the bigger would be the sequenced data file. Due to this reason, for instructional purposes, it is ideal to use a smaller genome.

### 1.1. Information of data we will be using in this tutorial.

In this tutorial, we will analyze the result of whole genome sequences of Mycobacterium tuberculosis isolates collected in New York State.

<aside>
💡 **SRR23086706**

</aside>

**WGS of Mycobacterium tuberculosis: isolate (SRR23086706)**

### 1.2. Two methods of downloading data

You can choose any one of these two methods.

#### 1.2.1. Directly using cloud - Fast

```bash
wget s3://sra-pub-src-13/SRR28409626/IDR1900024110-01-01.R1.fastq.gz.1
wget s3://sra-pub-src-13/SRR28409626/IDR1900024110-01-01.R2.fastq.gz.1
```

#### 1.2.2. Using SRA toolkit

```bash
fastq-dump --split-3 --gzip SRR28409626
```

> --split-3
> 
> - single-end reads will end up in a single file.
> - paired-end reads will produce two files
> - unpaired reads (if any) will be placed into a third file

---

## 2. Information on .fastq files

### 2.1. Viewing .fastq files content
``` 
head <file1_1.fastq> # To view top 10 lines
```

``` 
tail <file1_1.fastq> # To view bottom 10 lines
```

### 2.2 Determining the number of reads
```
zcat read.fastq.gz | grep @SRR | wc -l
```

### 2.3. fastqc

```
fastqc <file1_1.fastq> <file1_2.fastq> -q
```

> Reading fastqc report
> 
> - Per base sequence quality:
> - Per base sequence content
> - Per base GC content
> - Per base N content
> - Sequence Length distribution
> - Sequence Duplication levels
> - Overrepresented sequences
> - Adapter content

### 2.2. multiqc

```
multiqc .<input_dir> -o <out_dir>

```

---

## 3. Trimming reads

> Is it necessary?
Removal of low quality reads.
Reads containing adaptor sequences.
> 

Using #Triimmomatic
[Important resource](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf))

```
java -jar /trimmomatic-0.39.jar\\
 PE -phred33\\
 <input_file1_1> <input_file1_2>\\
 <output_file1_1> <output_file1_2>\\
 ILLUMINACLIP:adapter.fasta TRAILING:3 LEADING:3 MINLEN:36 SLIDINGWINDOW:4:15

```

Using #Cutadapt

### 3.2.1 To remove 3' adapter

**Single-end reads**

```
cutadapt -a <adapter_seq> <input.fastq.gz> -o <output.fastq.gz>

```

**Paired-end reads**

```
cutadapt -a <adapter_seq_forw> -A <adapter_seq_rev> <file1_1.fastq.gz> <file1_2.fastq.gz> <file1_1.trimmed.fastq.gz> <file1_2.trimmed.fastq.gz>

```

Using #Trim_Galore

```
trim_galore --paired --phred33 --gzip <file1_1.fastq.gz> <file1_2.fastq.gz>

```

---

## 4. Alignment of reads with the reference genome

### 4.1. Downloading reference genome

### 4.2. Indexing reference genome: one-time step

Using #samtools

```
samtools faidx <ref_genome>

```

Using #bwa : Works better with short (~100bp) reads.

```
bwa index <ref_genome>

```

This command generates
.amb
.ann
.bwt: Binary file
.pac: Binary file
.sa: Binary file

Using #bowtie2

```
bowtie2-build <ref_genome> <prefix>

```

This commands generates 6 files

```
<prefix>.1.bt2
<prefix>.2.bt2
<prefix>.3.bt2
<prefix>.4.bt2
<prefix>.rev.1.bt2
<prefix>.rev.2.bt2

```

### 4.3. Performing alignment generating SAM file

Using #bowtie2

```
bowtie2 -x <prefix> -1 <file1_1.fq> -2 <file1_2.fq> -S <read1.sam>

```

or
Using #bwa

```
bwa mem <ref_genome> <file1_1.fq> <file1_2.fq> > <read1.sam>

```

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
*Extracting information from VCF file (can be saved as a #bedfile),	```

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

Using #vcf2bed

```
vcf2bed < <file1.vcf> > <file1.bed>

```

*identify annotated genes which are not covered by reads across their full length.*

```
bedtools coverage -a <ref_genome.gff> -b <file.bam> -o <file1_gene-coverage.txt>

```

[[Generating an assembly with unmapped reads]]
