**Author**: Caylyn Railey

**Contact**: cer246@cornell.edu

**Overivew**:Script is going to be interpreted/run in bash, below is pseudo-coded workflow for identifying antisense lncRNAs from RNA-seq data. This workflow is optimized for stranded short-read data and ONT long-read data.

#### **Table of Contents**

- <font color="#7DF9FF">Pull sequencing reads from NCBI</font>
- <font color="#7DF9FF">Map sequencing reads</font>
  - Run Mapping Jobs
    - Stranded short read data
    - ONT long read data
-  <font color="#7DF9FF">Assemble bam files</font>
   - Assemble each bam file from the stranded short read data
   - Assemble each bam file from the ONT long read data
 - <font color="#7DF9FF">Merge all assembly files</font>
   - Merge all assembly files to generate a "new" annotation
   - Identify new genomic elements and antisense transcripts
 - <font color="#7DF9FF">Identify long non-coding RNAs from antisense trancripts</font>
   - Generate transcript nucleotide sequences
   - Run CPC2
   - Identify optimal ORFs
- <font color="#7DF9FF"> Generate genome annotation file and genome fasta file for future uses</font>

## Pull fastq files from <a href="https://www.ncbi.nlm.nih.gov/">NCBI</a> using <a href="https://sra-explorer.info/">SRA Explorer</a> 
  -- Fastest option is to use Aspera to download fastq files

  -- Simplest option is to use the curl command to download fastq files

***Example Aspera command***
```
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR592/007/SRR5920407/SRR5920407_1.fastq.gz . && mv SRR5920407_1.fastq.gz SRR5920407_GSM2736324_Root_without_salt_treatment_2_Camelina_sativa_RNA-Seq_1.fastq.gz 
```

***Example Curl command***
```
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR212/073/SRR21227073/SRR21227073_1.fastq.gz -o SRR21227073_GSM6504938_gi-22-2_Brassica_rapa_RNA-Seq_1.fastq.gz
```
## Run mapping jobs using <a href="https://daehwankimlab.github.io/hisat2/manual/">Hisat2</a> and <a href="https://lh3.github.io/minimap2/minimap2.html">minimap2</a>

Short-read data will be mapped using Hisat2 and long-read data will be mapped using minimap2. If long-read data is multiplexed, refer to this <a href="https://www.protocols.io/view/demultiplexing-nanopore-reads-with-last-14egnxw4zl5d/v7">tutorial</a> for prepping long-read data for mapping.

### Stranded Short-read RNA seq

1) Build Index ***example***

```
/path/to/hisat/hisat2-2.2.1/extract_exons.py your_species_genome_annotation.gff3 > xxxxx.exon
 
/path/to/hisat/hisat2-2.2.1/extract_splice_sites.py your_species_genome_annotation.gff3 > xxxxxx.ss
 
/path/to/hisat/hisat2-2.2.1/hisat2-build --ss xxxxxx.ss --exon xxxxx.exon /path/to/your_species_genome.toplevel.fa xxxxx_hisat -p 32
```

2) Map reads ***example***
```
cd /path/to/downloaded/fastq/

for i in *.fastq.gz;
do
name=$(basename ${i} _.fastq.gz);
/path/to/hisat/hisat2 \
-x /location/of/genome_index \
-1 ${name}_.fastq.gz \
-2 ${name}_.fastq.gz \
--max-intronlen 10000 \
--dta-cufflinks \
--rna-strandness RF \
--summary-file /path/to/mapping_summary/summary.txt \
-p 18 | samtools sort -@ 12 -o /path/to/desired/output/of/bam/${name}.bam -; done
```


### ONT long-read RNA seq

1) Build Index ***example***

```
minimap2 -d <species>_ rep<#>.idx -Q -t 10 -x splice your_species_genome.fa
```

2) Map reads ***example***
```
	mkdir -p mapped
	for bc in $(awk '{print $2}' <species>_ rep<#>_barcode_counts.txt);
	 do echo ${bc};
	 minimap2 -t 10 -a -x splice <species>_ rep<#>.idx oriented/${bc}_reads_dirAdjusted.fastq.gz | \
	 samtools view -b | samtools sort > mapped/mm2_called_$
	{bc}_vs_<species>_ rep<#>.bam;
    done
```

We have bam files now!

## Assemble bam files into gtf using <a href="https://ccb.jhu.edu/software/stringtie/">StringTie</a>

First, assemble each bioproject individually. Then you will merge the assemblies at the final step. Also, assemble stranded short-read data separately from long-read data. 

Assemble stranded short-read bam files ***example***

```
#!/bin/bash

cd /location/of/stranded_short_read_bam/files

for i in *.bam;
do
name=$(basename ${i} .bam);
/location/of/stringtie \
${name}.bam \
-p 14 \
-m 200 \
-f 0.05 \
-j 7 \
-c 10 \
-s 10 \
--rf \
-G /path/to/your_species_genome_annotation.gff3 \
-o /desired/path/to/output/assemblies/${name}.gtf
done
```

Assemble ONT long-read bam files ***example***

```
#!/bin/bash

cd /location/of/ONT_long_read_bam/files

for i in *.bam;
do
name=$(basename ${i} .bam);
/location/of/stringtie \
${name}.bam \
-p 14 \
--fr \
-L \
-m 200 \
-f 0.05 \
-j 5 \
-c 5 \
-s 8 \
-g 25 \
-G /path/to/your_species_genome_annotation.gff3 \
-o desired/path/to/output/assemblies/${name}.gtf
done
```

## Merge all assembly files to generate a "new" annotation
<a href="https://ccb.jhu.edu/software/stringtie/">StringTie</a> will be used to merge the assemblies made from each bam file. <a href="https://ccb.jhu.edu/software/stringtie/gffcompare.shtml">GffCompare</a> will the be used to identify new genomic elements.

***Example merge assemblies command***

```
#!/bin/bash

cd desired/path/to/output/assemblies/

/location/of/stringtie --merge \
-p 16 \
-G /path/to/your_species_genome_annotation.gff3 \
*.gtf \
-m 200 \
-o desired/path/to/output/merged_assembly/xxxxx.gtf
```

***Example GffCompare command to identify new genomic elements***

Run gffcompare for each species against reference annotation (merged.gtf vs mod.gtf) 
```
#!/bin/bash

/location/of/gffcompare \
-r /path/to/your_species_genome_annotation.gff3 \ 
/desired/path/to/output/merged_assembly/xxxxx.gtf \
-o /desired/path/to/output/gffcompare_results
```

<font color="#A03472">Within the gffcompare_results there will be multiple files. Right now, the file that we are most interested in the the *.tmap file. Using the tmap output, filter for antisense transcripts and antisense-intronic transcripts.</font>

Filtering can easily be done in Rstudio using Rscript provided (located on github): antisense_filtering.Rmd

Combine the antisense and antisense intronic transcript id files for further analysis

```
cat antisense_gffcompare_txIDs.txt antisense_intron_gffcompare_txIDs.txt > all_antisense_gffcompare_txIDs.txt
```

## Identify long non-coding RNAs from antisense trancripts using <a href="https://github.com/gao-lab/CPC2_standalone">CPC2</a> and  <a href="https://github.com/TransDecoder/TransDecoder">TransDecoder</a>

CPC2 and TransDecoder need input files at the nucleotide level to determine if a particular transcript has a high probability of being translated or not. In order to generate nucleotide sequences of all the antisense transcripts we can use <a href="https://ccb.jhu.edu/software/stringtie/gff.shtml">Gffread</a> 


1) Get coordinates for antisense transcripts

    Pull genomic information for all the antisense transcripts 

    *working in the following directory :/desired/path/to/output/gffcompare_results

```
grep -wFf  all_antisense_gffcompare_txIDs.txt ../all_merged.gtf > antisense_gffcompare_txIDs.gtf
```
2) Extract the entire transcript sequence using the generated genome annotation file and a genome fasta file

```
gffread -w antisense_gffcompare_txIDs.fa -g /mnt/Kanta/antisense_lncRNAs/Camelina//path/to/your_species_genome.toplevel.fa antisense_gffcompare_txIDs.gtf
```


### Run CPC2 to determine coding capacity

```
/mnt/Kanta/antisense_lncRNAs/CPC2_standalone-1.0.1/bin/CPC2.py \
-i antisense_gffcompare_txIDs.fa \
-o /desired/path/to/output/CPC2_output/

```

#### Test antisense transcripts for homology with known proteins using TransDecoder and Pfam

Now, you might need to download Pfam to your working environment. You can do so with the following command:
```
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.full.gz
```

Inside the Pfam folder follow the directions in the ReadMe to fully install the database to your working environment.


I would suggest making a new folder inside your gffcompare results and labeling it as Pfam

*Working directory: /desired/path/to/output/gffcompare_results/Pfam

1) Translate the longest orf using Transdecoder

```
TransDecoder.LongOrfs -t ../antisense_gffcompare_txIDs.fa -m 50 -S
```

*the output at this step is a new folder that ends with ".transdecoder_dir"*

2) Blast the longest ORF against the Pfam database 


```
cd /*.transdecoder_dir

path/to/PfamScan/pfam_scan.pl -fasta longest_orfs.pep -dir /mnt/Kanta/antisense_lncRNAs/Pfam_files -outfile pfam_antisense.txt -cpu 12
```

### Use the provided Rscript to filter transcripts based on coding capacity for antisense lncRNAs
RScript provided: cpc2_pfam_output_filtering.Rmd


## Generate genome annotation file and genome fasta file for future uses

*Generate genome and genome annotation with HC lncRNAs*

```
grep -wFf hc_anti_txIDs.txt ../gffcompare.annotated.gtf > hc_antisense.gtf

cat  /path/to/your_species_genome_annotation.gtf hc_antisense.gtf > ref_plus_hc_antisense_annotation.gtf

gffread -w ref_plus_hc_antisense_annotation.fa -g /path/to/your_species_genome.toplevel.fa ref_plus_hc_antisense_annotation.gtf
```

*Generate genome and genome annotation with LC lncRNAs*

```
grep -wFf lc_anti_txIDs.txt ../gffcompare.annotated.gtf > lc_antisense.gtf

cat /path/to/your_species_genome_annotation.gtf lc_antisense.gtf > ref_plus_lc_antisense_annotation.gtf

gffread -w ref_plus_lc_antisense_annotation.fa -g /path/to/your_species_genome.toplevel.fa ref_plus_hc_antisense_annotation.gtf
```

*Generate genome and genome annotation with HC+LC lncRNAs*

```
cat hc_antisense.gtf lc_antisense.gtf > hc_lc_antisense.gtf

cat /path/to/your_species_genome_annotation.gtf hc_lc_antisense.gtf > ref_plus_hc_lc_antisense_annotation.gtf

gffread -w  ref_plus_hc_lcantisense_annotation.fa -g /path/to/your_species_genome.toplevel.fa   ref_plus_hc_lc_antisense_annotation.gtf
```

### <font color="#A03472">Now you have a genome annotation and genome file that contain newly identified antisense lncRNAs. These genome and genome annotation files can be used in furhter inquires in expression and transcript architecture!</font> ###












