# Tutorial: Next Generation Sequence Analysis - Stepwise

**This repository attempts to provide a stepwise tutorial on the Next Generation Sequence (NGS) analysis.**

In biological contexts, sequencing refers to determining the order of nucleotides (A, T, G, C) in a DNA sample. The term next generation sequencing (NGS) refers to an advanced and high-throughput methodnology of determining nucleotides order.

> **What about sequencing of nucleotides from RNA sample?**
>
> Similar to the DNA sequencing, RNA sequencing also does exist. While dealing with RNA samples, for most biological interpretations RNA is converted into an intermediate form known as cDNA (complementary DNA). However, methods namely direct RNA sequencing (DRS) involves sequencing of RNA directly.

## The goals of NGS analysis can be

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

The best way to follow this tutorial is to clone this repository in your local device.

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

#### 0.4.4. Setting working directory

Now to avoid writing long code chunks especially long directory path and file names, we will use variables.

```{bash}
working_dir=$(pwd)          #Here we have saved working directory as a variable
```

---

## 1. Downloading sequencing data

As the complete sequence analysis is a multi-step and time-consuming procedure, we will be using a smaller dataset for the learning purposes. In this tutorial, we will be using example dataset featured by the [Qiagen](https://digitalinsights.qiagen.com/downloads/example-data/).\
The sequencing data comes from sequencing sample of *Pseudomonas aeruginosa* species.

```{bash}
mkdir -p $working_dir/fastq_reads       #Creating directory to save fastq_reads

wget -L -O $working_dir/fastq_reads/SRR396636.sra_1.fastq https://zenodo.org/records/11791175/files/SRR396636.sra_1.fastq?download=1
wget -L -O $working_dir/fastq_reads/SRR396636.sra_2.fastq https://zenodo.org/records/11791175/files/SRR396636.sra_2.fastq?download=1
wget -L -O $working_dir/fastq_reads/SRR396637.sra_1.fastq https://zenodo.org/records/11791175/files/SRR396637.sra_1.fastq?download=1
wget -L -O $working_dir/fastq_reads/SRR396637.sra_2.fastq https://zenodo.org/records/11791175/files/SRR396637.sra_2.fastq?download=1
```

### 1. Saving fastq files as variables

```{bash}
read1_F=$working_dir/fastq_reads/SRR396636.sra_1.fastq
read1_R=$working_dir/fastq_reads/SRR396636.sra_2.fastq
read2_F=$working_dir/fastq_reads/SRR396637.sra_1.fastq
read2_R=$working_dir/fastq_reads/SRR396637.sra_2.fastq
```

## 2. Information on .fastq files

Once the sequencing data has been downloaded, with the listing command the 4 `.fastq` files can be seen in the current directory.

```{bash}
ls -1 $working_dir/fastq_reads
```

### 2.1. Viewing .fastq files content

To view content of the `.fastq` file, there are multiple ways

#### 2.1.1. Using `cat` or `less` command

Use `head` to view top 10 lines or `tail` to view bottom 10 lines of a `.fastq` file.

```{bash}
cat $read1_F | head -n 10       #cat command often print the output in the terminal
```

```{bash}
less $read1_F | head -n 10      #less command is used to view the output without printing in the terminal
```

>![alt text](<Screenshot 2024-06-24 123033.png>)
> In fastq files, there are four lines corresponding to a read.

### 2.2 Determining the number of reads

To determine the number of reads in a `fastq` file, use the below command.

```{bash}
cat $read1_R | grep @SRR | wc -l        #If you encounter compressed fastq files (such as file.fastq.gz or file.fq.gz), cat should be replaced with zcat
```

> Here,\
> Both `cat` and `less` command can be used. \
> `grep` command works as search tool and searches for occurence of word `@SRR` in the `.fastq` file.\
> `wc -l` command count the occurences of `@SRR`.

🔍 **Check yourself whether the mate file (i,e, SRR396636.sra_2.fastq) has the same no. of reads.**

### 2.3. Fastqc

Fastqc is a tool developed by **Babraham Bioinformatics** for looking at quality of fastq reads.

```{bash}
mkdir -p $working_dir/fastq_reads/fastqc_results        #Creates a folder to save the fastqc results
fastqc *.fastq -o $working_dir/fastq_reads/fastqc_results           #Here, *.fastq works as a wildcard for us.
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
mkdir -p $working_dir/fastq_reads/multiqc_results       #Creates a folder to save multiqc_results
multiqc /fastqc_results/ . -o $working_dir/fastq_reads/multiqc_results
```

## 3. Trimming reads

In sequencing, reads are flagged with adapter sequences (either at one or both ends). Being complementary to the adaper/oligo sequences attached to the flow cells, they help in initating the sequencing-by-synthesis (by acting as primers) using the polymerase enzyme.

However, after the sequencing adapter sequences should be removed for the following reasons;

1. In the next step of alignment of reads with the reference genome, the quality of alignment might be hamperred as adapter sequences does not belong to the species genome.

This step of trimming low quality reads and removal of adapters could be done by any of these three tools.

<details>
<summary>Using <b>Trimmomatic</b></summary>

> [Important resource](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

```{bash}
mkdir -p trimmomatic_results
trimmomatic PE -phred33 \
 $read1_R $read1_R \
 SRR396636.sra_1_trim.fastq SRR396636.sra_2_trim.fastq \
 ILLUMINACLIP:adapter.fasta TRAILING:3 LEADING:3 MINLEN:36 SLIDINGWINDOW:4:15
```

</details>

<details>
<summary>Using <b>Cutadapt</b></summary>

```{bash}
mkdir -p cutadapt_results
cutadapt -a <adapter_seq_forw> -A <adapter_seq_rev> $read1_R $read1_R SRR396636.sra_1_trim.fastq SRR396636.sra_2_trim.fastq -o cutadapt_results/
```

> [Documentation](https://cutadapt.readthedocs.io/en/v1.9/guide.html)

> Use: \
> -a [adapter] : To remove 3' [adapter] sequence
> -a file:[adapters.fasta] : To read adapter sequences from a fasta file \
> -g [adapter] : To remove 5' [adapter] sequence \
> -u [n] : To remove 'n' bases from the begining of each read. Replace by '-n' to remove bases from the end of each read. \
> -q [n] : To trim read from 3' end, if quality of read drops below [n] \
> -m and -M [n] : To remove reads having length smaller than [m] or greater than [M]

</details>

<details>
<summary>Using <b>Trim Galore</b></summary>

[Documentation](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)

```{bash}
mkdir -p trim_galore_results
trim_galore --paired --phred33 --gzip \
$read1_R $read1_R -o trim_galore_results/

mkdir -p fastqc_after_trim_galore           #Creating directory for running fastqc on trimmed reads
--fastqc_args --outdir /fastqc_after_trim_galore
```

</details>

---

## 4. Alignment of reads with the reference genome

As name indicates, a reference genome is the complete genome sequence of a species. Such a genome is annotated based on evidences. Reference genomes gets updated regularly according to the updated knowledge.

### 4.1. Downloading reference genome

For our project, there are multiple methods to download reference genome. If you are already aware about downloading reference genomes or simply not in the mood, use the below code and move to section 4.2

```
mkdir -p $working_dir/reference_genome       #Creating a directory to store reference genome 
ref_genome="GCF_000006765.1_ASM676v1_genomic.fna.gz"        #Saving the reference_genome file as a variable named `ref_genome` for ease

wget -O reference_genome/$ref_genome https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz      #Downloading reference genome in the directory we created.
```

#### 4.1.1. Learn how to download reference genome

1. Go to [NCBI homepage](https://ncbi.nlm.nih.gov).
2. Select `assembly` from the drop-down menu and search for `Pseudomonas aeruginosa`.

    > ![alt text](image-4.png)

3. From the filters, under `RefSeq category`, select `Reference`

    > ![alt text](image-3.png)

4. Click on `ASM676v1` assembly

    > ![alt text](image-5.png)

**Now you have four methods to download the reference_genome**.

> ![alt text](image-6.png)

<details>
    <summary>Method 1: Using <b>Download</b></summary>
<aside>

1. Click `Download`.
2. Among file sources select `RefSeq only`, and among file types make sure `Genome sequences (FASTA)` being selected.
3. This will download a zip folder that would needed be move to current working directory (for working ease) and would needed to be uncompressed.

</aside>
</details>

<details>
        <summary>Method 2: Using <b>Datasets</b></summary>
    <aside>

1. Installing `datasets` conda package

    ```{bash}
    conda install -c conda-forge ncbi-datasets-cli
    ```

2. Click on `datasets`
3. Copy the command

    ```{bash}
    datasets download genome accession GCF_000195955.2 --include genome         #Here I have modified the copied command to download only genome file.
    ```

4. Paste in the command line to download the refernce genome.

</aside>
</details>

<details>
<summary>Method 3: Using <b>URL</b></summary>
<aside>

1. This method allows you to simply copy the provided URL and paste in your browser to download all the files.
2. Certainly, you would need to unzip the folder to be able to use the files.

</aside>
</details>

<details>
    <summary>Method 4: Using <b>FTP (Best practice)</b></summary>
<aside>

1. click on `FTP`.
2. Right click on the file `GCF_000006765.1_ASM676v1_genomic.fna.gz` and copy the link.

    ![alt text](<Screenshot 2024-06-22 191142.png>)

3. Use the `wget` or `curl` command.

    ```{bash}
    mkdir -p reference_genome          #Create a directory to save the reference genome.
    wget -o reference_genome https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz          #In place of wget, curl command can also be used.
    ref_genome="GCF_000006765.1_ASM676v1_genomic.fna.gz"
    ```

</aside>
</details>

---

### 4.2. Indexing reference genome: one-time step

**Indexing a reference genome is a one-time step and can be performed using any of below tools. However, keep in mind that the tool used for generating reference index should also be used for alignment.**

<details>
<summary>Using <b>HISAT2</b></summary>
<aside>

```{bash}
mkdir -p reference_genome/hisat2_index
hisat2-build $ref_genome hisat2_index/[prefix]
```

> This commmand generates the 8 files. \
> [prefix].1-8.ht2

</aside>
</details>

<details>
    <summary>Using <b>samtools</b></summary>
<aside>

```{bash}
mkdir -p reference_genome/samtools_index
samtools faidx $ref_genome -o samtools_index/[prefix].fasta.fai
```

</aside>
</details>

<details>
    <summary>Using <b>bwa</b></summary>
<aside>

```{bash}
mkdir -p reference_genome/bwa_index
bwa index $ref_genome
```

> ![alt text](image-2.png)

> This command generates the following 5 files: \
> .amb \
> .ann \
> .bwt: Binary file \
> .pac: Binary file \
> .sa: Binary file

</aside>
</details>

<details>
    <summary>Using <b>bowtie2</b></summary>
<aside>

```{bash}
mkdir -p reference_genome/bowtie2_index
bowtie2-build $ref_genome <prefix>
```

> This command generates 6 files \
> [prefix].1.bt2 \
> [prefix].2.bt2 \
> [prefix].3.bt2 \
> [prefix].4.bt2 \
> [prefix].rev.1.bt2 \
> [prefix].rev.2.bt2
  
</aside>
</details>

---

### 4.3. Performing alignment generating sequence alignment map (SAM) file

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
mkdir -p bwa_alignment
bwa mem $ref_genome $read1_F $read1_R > ../bwa_alignment/SRR396636.sam
```

</aside>
</details>

<details>
    <summary>Using <b>bowtie2</b></summary>
<aside>

```{bash}
bowtie2 -x <prefix> -1 $read1_R -2 $read1_R -S alignment_results/SRR396636.sam
```

</aside>
</details>

---

## 5 Analysis of alignment results

### 5.1 Converting SAM file into binary alignment map (BAM) file

The SAM file is human readable form

```{bash}
samtools view alignment_results/SRR396636.sam
```

<details>
<summary><b>Understand SAM file format</b></summary>
<aside>

> SAM file has three header lines starting with @
>
>1. @HD \
>
- VN: Version of SAM \
- SO: Sorted/Unsorted \
- GO:
>
>2. @SQ: Reference seq \
>
- SN: Ref seq accession \
- LN: Length of Ref seq
>
>3. @PG: Programme used to generate SAM file \
>
- ID: ID of program \
- PN: Name of program \
- VN: Version of program \
- CL: Wrapper script used to generate SAM

> ![alt text](image-8.png)

</aside>
</details>

Using `samtools`

```{bash}
samtools view -Sb alignment_results/SRR396636.sam -o alignment_results/SRR396636.bam
```

> -Sb : Input in SAM format (S) and the output will be be BAM format(b)

### 5.2 Sorting BAM file

Using `samtools`

```{bash}
samtools sort SRR396636.bam -o alignment_results/SRR396636.sort.bam
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
samtools flagstat alignment_results/SRR396636.sort.bam > alignment_results/alignment_summary.txt
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

<details>
<summary>Using <b>Freebayes</b></summary>

```{bash}
freebayes -f <ref_genome> <file1.bam> > <file1.vcf>
```

</details>

<details>
<summary>Using <b>Bcftools</b></summary>
<aside>

- **Step 1: Generating pileup file**
Generates genotype likelihoods at each genomic position with coverage.

```{bash}
mkdir -p variant_call_results
bcftools mpileup -O b -f $ref_genome alignment_results/*.sort.bam -o variant_call_results/SRR396636.mpileup.bcf
```

- **Step 2: Detecting SNVs**

```{bash}
bcftools call -O b \            #Output file type will be bcf
--vc \
--ploidy 1 \
variant_call_results/SRR396636.mpileup.bcf -o variant_call_results/SRR396636.call.bcf
```

- **Step 3: Filtering variants**

```{bash}
bcftools filter -O v \          #Output file type will be vcf
-i '%QUAL>=n' \
variant_call_results/SRR396636.call.bcf -o variant_call_results/SRR396636.call.filtered.vcf

```

</aside>
</details>

In a single step

```{bash}
bcftools mpileup -Ou -f <ref_genome> <file1.bam> | bcftools call -mv -Ob -o <file1.bcf>
```

> **-O: output type**
>> -b: compressed BCf, -u: uncompressed BCf \
>> -z: compressed VCF, -v: uncompressed VCF

---

## 6. Analyzing VCF file

### 6.1: Viewing VCF file

```{bash}
bcftools view <file1.vcf>
```

### 6.2: Determining the number of variant sites

```{bash}
grep -v -c "^#" <file1.vcf>
```

### 6.3: Generating an index file for VCF to view in IGV

Useful commands
*Extracting information from VCF file (can be saved as a #bedfile)

```{bash}
bcftools query -f '%COLa\\t%COLb\\t%COLc\\n' <file1.vcf>
```

*Converting BCF file to VCF file*

```{bash}
bcftools view <file1.bcf> > <file1.vcf>
```

### 6.4: Sorting VCF file by chromosome → coordinate

```{bash}
vcf-sort <file1.vcf> > <file1.sorted.vcf>
```

---

### 6.5: Converting VCF to BED file

Using `vcf2bed`

```{bash}
vcf2bed < <file1.vcf> > <file1.bed>
```

*identify annotated genes which are not covered by reads across their full length.*

```{bash}
bedtools coverage -a <ref_genome.gff> -b <file.bam> -o <file1_gene-coverage.txt>
```

---
