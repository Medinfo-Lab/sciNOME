**fastp**：version 0.23.4

**STAR**：version 2.7.5a_2020-06-29

**featureCounts**：version 2.0.1

**bismark**：version 0.24.2

**bowtie2**：version 2.3.5.1

# scRNA-seq

## 1.fastp quality control

```bash
files=($(find RNA_fq/paired -maxdepth 1 -name "*_1.fastq.gz" -printf "%f\n"))
for file in "${files[@]}"
do
	#echo ${file}
	name=${file%_*}
	echo ${name}
	fastp --in1 RNA_fq/paired/${name}_1.fastq.gz --in2 RNA_fq/paired/${name}_2.fastq.gz --detect_adapter_for_pe --cut_right --out1 RNA_trim_fq/paired/${name}_R1.clean.fastq.gz --out2 RNA_trim_fq/paired/${name}_R2.clean.fastq.gz --html qc_reports/${name}_fastp.html --json qc_reports/${name}_fastp.json
done

files=($(find RNA_fq/single -maxdepth 1 -name "*.fastq.gz" -printf "%f\n"))
for file in "${files[@]}"
do
	#echo ${file}
	name=${file%%.*}
	echo ${name}
	fastp --in1 RNA_fq/single/${name}.fastq.gz --detect_adapter_for_pe --cut_right --out1 RNA_trim_fq/single/${name}.clean.fastq.gz --html qc_reports/${name}_fastp.html --json qc_reports/${name}_fastp.json
done
```

## 2.STAR alignment

```bash
for file in "${files[@]}"
do
	name=${file%_*}
	echo ${name}
	/work/tools/STAR/bin/Linux_x86_64/STAR --genomeDir ~/genomes/mm10_star --readFilesIn RNA_trim_fq/paired/${name}_R1.clean.fastq.gz RNA_trim_fq/paired/${name}_R2.clean.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix bam/paired/${name}_
done

files=($(find RNA_trim_fq/single -maxdepth 1 -name "*.clean.fastq.gz" -printf "%f\n"))
for file in "${files[@]}"
do
	name=${file%%.*}
	echo ${name}
	/work/tools/STAR/bin/Linux_x86_64/STAR --genomeDir ~/genomes/mm10_star --readFilesIn RNA_trim_fq/single/${name}.clean.fastq.gz --runThreadN 3 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix bam/single/${name}_
done
```

## 3.featureCounts

```bash
/work/tools/featureCounts/subread-2.0.1-Linux-x86_64/bin/featureCounts -T 5 -p -s 0 -M --fraction -t exon -g gene_id -a /work/genomes/Gencode/GRCm38/gencode.vM25.annotation.gtf -o GSE121650_paired_counts_exon2.txt bam/paired/*Aligned.sortedByCoord.out.bam 

/work/tools/featureCounts/subread-2.0.1-Linux-x86_64/bin/featureCounts -T 5 -s 0 -M --fraction -t exon -g gene_id -a /work/genomes/Gencode/GRCm38/gencode.vM25.annotation.gtf -o GSE121650_single_counts_exon2.txt bam/single/*_Aligned.sortedByCoord.out.bam
```



# scNOMe-seq

## 1.fastp quality control

```bash
files=($(find NOMe_fq -maxdepth 1 -name "*_1.fastq.gz" -printf "%f\n"))
for file in "${files[@]}"
do
	#echo ${file}
	name=${file%_*}
	echo ${name}
	fastp --in1 NOMe_fq/${name}_1.fastq.gz --in2 NOMe_fq/${name}_2.fastq.gz -y 0 --cut_right --detect_adapter_for_pe --trim_front1 6 --trim_front2 6 --length_required 20 --out1 NOMe_trim_fq/${name}_R1.clean.fastq.gz --out2 NOMe_trim_fq/${name}_R2.clean.fastq.gz --html qc_reports/${name}_fastp.html --json qc_reports/${name}_fastp.json
done
```

## 2.bismark alignment

```bash
files=($(find NOMe_trim_fq -maxdepth 1 -name "*_R1.clean.fastq.gz" -printf "%f\n"))
for file in "${files[@]}"
do
	name=${file%_*}
	echo ${name}
	mkdir bismark/bam/${name}
	bismark --genome ~/genomes/mm10_bismark/mm10 --bowtie2 --paralle 4 --non_directional -X 1000 -p 4 --bam -1 NOMe_trim_fq/${name}_R1.clean.fastq.gz -2 NOMe_trim_fq/${name}_R2.clean.fastq.gz --output_dir bismark/bam/${name}
done
```

## 3.deduplicate_bismark deduplication

```bash
files=($(find bismark/bam -maxdepth 1 -printf "%f\n"))
for file in "${files[@]}"
do
	echo ${file}
	mkdir bismark/deduplicate_bam/${file}
	deduplicate_bismark --bam -p --output_dir bismark/deduplicate_bam/${file} bismark/bam/${file}/${file}_R1_bismark_bt2_pe.bam
done
```

## 4.methylation extraction

```bash
files=($(find bismark/deduplicate_bam -maxdepth 1 -printf "%f\n"))
for file in "${files[@]}"
do
	echo ${file}
	mkdir bismark/methylation_extractor/${file}
	bismark_methylation_extractor --gzip --bedGraph --CX --buffer_size 10G --genome_folder ~/genomes/mm10_bismark/mm10 -o bismark/methylation_extractor/${file} bismark/deduplicate_bam/${file}/${file}_R1_bismark_bt2_pe.deduplicated.bam
	#rm bismark/methylation_extractor/${file}/CHG_*
	#rm bismark/methylation_extractor/${file}/CHH_*
	#rm bismark/methylation_extractor/${file}/CpG_* 
done
```

## 5.separate CpG and GpC

```bash
files=($(find bismark/deduplicate_bam -maxdepth 1 -printf "%f\n"))
for file in "${files[@]}"
do
	echo ${file}
	mkdir bismark/coverage2cytosine/${file}
	coverage2cytosine --genome_folder ~/genomes/mm10_bismark/mm10 --nome-seq --gc --gzip --dir bismark/coverage2cytosine/${file} -o ${file} bismark/methylation_extractor/${file}/${file}_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz
done
```

