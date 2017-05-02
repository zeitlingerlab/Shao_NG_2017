# Amazon EC2 AMI setup notes

We start with the Ubuntu 16.04 EBS-backed AMI for the us-east-1 region.

## Boot AMI

AMI was started as type **c4xlarge** with a 250GB EBS boot volume. This instance type has 16 CPU cores and ~30GB of memory.

## Set up data directory

```bash
sudo mkdir /data
sudo chown -R ubuntu:ubuntu /data
```

## Install R 3.3.2 and Bioconductor 3.3

Add the line below to `/etc/apt/sources.list`

```
deb https://cloud.r-project.org/bin/linux/ubuntu xenial/
```

Execute the following to add the appropriate signing key:

```bash
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
```

Update packages and install R:

```bash
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install r-base r-base-dev
sudo apt-get install libxml2-dev libcurl3-dev # for Rsamtools
sudo apt-get install htop samtools pigz sra-toolkit parallel git bowtie bedtools python-pip python-numpy procmail
```

Create ~/.Rprofile: 

```S
update_bioc = function() {
  library(BiocInstaller)
  update.packages(repos=biocinstallRepos(), ask=FALSE)
}

update_all = function() {
  message("Updating R packages...")
  update.packages(ask=F)
  message("Updating Bioconductor packages...")
  update_bioc()
}
````

Install R and Bioconductor packages by sourcing the following script from a root R session:

```S
update.packages(ask=F)
install.packages(c("reshape", "ggplot2", "stringr", "optparse", "dplyr", "knitr", "rmarkdown", "magrittr"))

source("http://bioconductor.org/biocLite.R")
biocLite()

bioc_packages <- c("GenomicRanges", 
                   "BSgenome.Dmelanogaster.UCSC.dm3", 
                   "BSgenome.Hsapiens.UCSC.hg19",
                   "rtracklayer", 
                   "ShortRead",
                   "Rsamtools",
				           "BiocParallel",
				           "chipseq")

biocLite(bioc_packages)
````

## Download analysis code  

```bash
cd /data
git clone https://github.com/zeitlingerlab/Shao_NG_2017 analysis_code
```

## Install cutadapt 

```bash
mkdir /data/software
cd /data/software
git clone https://github.com/marcelm/cutadapt.git
cd cutadapt
 
pip install Cython
sudo python setup.py install
```

## Install SRA tools

```bash
mkdir /data/software/sra
cd /data/software/sra
wget 'https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.7.0/sratoolkit.2.7.0-ubuntu64.tar.gz'
tar -xf sratoolkit.2.7.0-ubuntu64.tar.gz

# Accept defaults
/data/software/sra/sratoolkit.2.4.2-ubuntu64/bin/vdb-config -i
```

## Build UCSC dm3 reference genome for Bowtie

```bash
mkdir /data/genomes
cd /data/genomes
mkdir dm3
cd dm3
wget 'http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz'
tar -xf chromFa.tar.gz
rm chromFa.tar.gz
cat *.fa > dm3.fasta
rm *.fa
bowtie-build -o 1 dm3.fasta dm3
```

## Build UCSC dm3-hg19 reference genome for Bowtie 

```bash
cd /data/genomes
mkdir hg19
cd hg19
wget 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'
tar -xf chromFa.tar.gz
rm chromFa.tar.gz
cat *.fa > hg19.fasta
rm *.fa

cd /data/genomes
mkdir  dm3_hg19
cd dm3_hg19

sed 's/>/>hg19_/' ../hg19/hg19.fasta > hg19_renamed.fasta
sed 's/>/>dm3_/' ../dm3/dm3.fasta > dm3_renamed.fasta

cat  dm3_renamed.fasta hg19_renamed.fasta > dm3_hg19.fasta
rm dm3_renamed.fasta
rm hg19_renamed.fasta

bowtie-build dm3_hg19.fasta dm3_hg19
```

## Preprocess ChIP-nexus FASTQ reads 

```bash
mkdir /data/preprocessed_fastq
cd /data/preprocessed_fastq

### For samples with single fixed barcode and are 50bp sequenced
parallel -uj 2 Rscript /data/analysis_code/pipeline/preprocess_fastq_he_modified.r -f {} -k 22 -b CTGA -r 5 -p 5 -o \`basename {} .fastq.gz\`_processed.fastq.gz ::: /data/raw_fastq/chipnexus_single_barcode/*.fastq.gz   

### For samples with mixed fixed barcode and are 75bp sequenced
parallel -uj 2 Rscript /data/analysis_code/pipeline/preprocess_fastq_he_modified.r -f {} -k 22 -b CTGA,TGAC,GACT,ACTG -t 50 -r 5 -p 5 -o \`basename {} .fastq.gz\`_processed.fastq.gz ::: /data/raw_fastq/chipnexus_mixed_barcode/*.fastq.gz 


mkdir /data/preprocessed_fastq/regular_chipnexus
mkdir /data/preprocessed_fastq/spikein_chipnexus

mv *spikein*.fastq.gz ./spikein_chipnexus/
mv *.fastq.gz  ./regular_chipnexus/
```


## Align ChIP-nexus reads

```bash

mkdir -p /data/bam/regular_chipnexus
cd /data/bam/regular_chipnexus

parallel -uj 2 /data/analysis_code/pipeline/align_chipnexus.sh {} /data/genomes/dm3/dm3 ::: /data/preprocessed_fastq/regular_chipnexus/*.fastq.gz
parallel -uj 2 /data/analysis_code/pipeline/sort_bam {} ::: *.bam


mkdir -p /data/bam/spikein_chipnexus
cd /data/bam/spikein_chipnexus

parallel -uj 2 /data/analysis_code/pipeline/align_chipnexus.sh {} /data/genomes/dm3_hg19/dm3_hg19 ::: /data/preprocessed_fastq/spikein_chipnexus/*.fastq.gz
parallel -uj 2 /data/analysis_code/pipeline/sort_bam {} ::: *.bam
```

## Process aligned ChIP-nexus reads

```bash
mkdir -p /data/rdata/granges/regular_chipnexus
cd /data/rdata/granges/regular_chipnexus

parallel -uj 10 Rscript /data/analysis_code/pipeline/process_chipnexus_bam.r -f {} -n \`basename {} .bam\` ::: /data/bam/regular_chipnexus/*.bam
	
mkdir -p /data/rdata/granges/spikein_chipnexus
cd /data/rdata/granges/spikein_chipnexus	
	
parallel -uj 10 Rscript /data/analysis_code/pipeline/process_chipnexus_bam.r -f {} -n \`basename {} .bam\` ::: /data/bam/spikein_chipnexus/*.bam


```


## Generate normalized ChIP-nexus BigWigs

```bash
mkdir /data/bigwig

cd /data/rdata/granges/regular_chipnexus

### For regular ChIP-nexus samples  

ls -1 $PWD/*_1.granges.rds | sed 's/_1.*//' > nexus_tobe_merged.txt
ls -1 $PWD/*.granges.rds |grep -e _1 -e _2 -v | sed 's/\.granges.*//' > nexus_no_merge.txt
ls $PWD/kc167_dmso_polii_chipnexus_* | sed 's/\.granges.*//' >> nexus_no_merge.txt
ls $PWD/kc167_dmso_tfiia_chipnexus_* | sed 's/\.granges.*//' >> nexus_no_merge.txt

cd /data/bigwig/

parallel -uj 3 Rscript /data/analysis_code/pipeline/generating_bw_from_gr.r -f {}_1.granges.rds -s {}_2.granges.rds -n \`basename {}\`_normalized_and_merged -t chipnexus :::: /data/rdata/granges/regular_chipnexus/nexus_tobe_merged.txt 
parallel -uj 3 Rscript /data/analysis_code/pipeline/generating_bw_from_gr.r -f {}.granges.rds -n \`basename {}\`_normalized -t chipnexus :::: /data/rdata/granges/regular_chipnexus/nexus_no_merge.txt 

Rscript /data/analysis_code/pipeline/generating_bw_from_gr.r -f kc167_dmso_polii_chipnexus_2.granges.rds -n \`basename {}\`_normalized -t chipnexus


### For spikein samples
cd /data/rdata/granges/spikein_chipnexus

ls -1 $PWD/*.granges.rds |grep -e egfp -e nelf -v  > half_life_samples.txt
ls -1 $PWD/*.granges.rds |grep -e egfp -e nelf | grep -e _1 | sed 's/_1.*//'   > knockdown_samples.txt


cd /data/bigwig/

parallel -uj 3 Rscript /data/analysis_code/pipeline/spikein_normalization.r -f {} -n \`basename {} .granges.rds\` :::: /data/rdata/granges/spikein_chipnexus/half_life_samples.txt
parallel -uj 2 Rscript /data/analysis_code/pipeline/spikein_normalization.r -f {}_1.granges.rds -d {}_2.granges.rds -n \`basename {} \`_normalized_and_merged :::: /data/rdata/granges/spikein_chipnexus/knockdown_samples.txt


rename 's/spikein_1/spikein_normalized_1/' *.bw
rename 's/spikein_2/spikein_normalized_2/' *.bw

```

## Align ChIP-seq samples

Align all ChIP-seq samples:

```bash
mkdir -p /data/bam/chipseq
cd /data/bam/chipseq

parallel -uj 3 /data/analysis_code/pipeline/align_chipseq.sh {} /data/genomes/dm3/dm3 ::: /data/raw_fastq/chipseq/*fastq.gz
parallel -uj 2 /data/analysis_code/pipeline/sort_bam {} ::: *chipseq*.bam
````

## R object generation for ChIP-seq samples

Generate R objects for all aligned samples:

```bash
cd /data/bam/chipseq
Rscript /data/analysis_code/pipeline/process_chipseq_bam.r -f kc167_polii_chipseq.bam -e 198 -n kc167_polii_chipseq
Rscript /data/analysis_code/pipeline/process_chipseq_bam.r -f kc167_tfiia_chipseq.bam -e 230 -n kc167_tfiia_chipseq
Rscript /data/analysis_code/pipeline/process_chipseq_bam.r -f kc167_tfiib_chipseq.bam -e 153 -n kc167_tfiib_chipseq

mkdir /data/rdata/granges/chipseq
mv *.RData /data/rdata/granges/chipseq/
```
## Generate normalized ChIP-seq BigWigs
```bash
cd /data/bigwig
parallel -uj 3 Rscript /data/analysis_code/pipeline/generating_bw_from_gr.r -f {}  -n \`basename {} .ranges.RData\`_normalized -t chipseq ::: /data/rdata/granges/chipseq/*.ranges.RData

```

## Download Drosophila PRO-cap data (Lis lab)

```bash
mkdir /data/public_data/
cd /data/public_data/

mkdir ./lis_procap
cd ./lis_procap/

wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1032nnn/GSM1032759/suppl/GSM1032759_PROcap.mn.bedgraph.gz'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1032nnn/GSM1032759/suppl/GSM1032759_PROcap.pl.bedgraph.gz'

```

## Download expression data and annotation from flybase

```bash
cd /data/public_data/

mkdir ./flybase
cd ./flybase/

wget 'ftp://ftp.flybase.net/releases/FB2014_03/precomputed_files/genes/gene_rpkm_report_fb_2014_03.tsv.gz'
wget 'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/gff/dmel-all-no-analysis-r5.57.gff.gz'

zcat ./gene_rpkm_report_fb_2014_03.tsv.gz | grep Kc167 > kc167_rna.txt
awk '{ if ($3 == "mRNA") print }' <(zcat ./dmel-all-no-analysis-r5.57.gff.gz) >flybase_r5.57_mrna.txt
```


## Download TBP data (Zeitlinger lab)

```bash
cd /data/public_data/
mkdir ./zeitlinger_tbp

cd ./zeitlinger_tbp/

wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1333nnn/GSM1333891/suppl/GSM1333891_hsap_k562_tbp_chipnexus_01_negative.bw'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1333nnn/GSM1333891/suppl/GSM1333891_hsap_k562_tbp_chipnexus_01_positive.bw'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1333nnn/GSM1333892/suppl/GSM1333892_hsap_k562_tbp_chipnexus_02_negative.bw'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1333nnn/GSM1333892/suppl/GSM1333892_hsap_k562_tbp_chipnexus_02_positive.bw'
```

## Download human annotation and K562 gro-cap data (Lis lab)

```bash
cd /data/public_data/
mkdir ./human_annotation
cd ./human_annotation

wget 'http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz'

cd /data/public_data/
mkdir ./lis_grocap
cd ./lis_grocap/

wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1480nnn/GSM1480321/suppl/GSM1480321_K562_GROcap_wTAP_minus.bigWig'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1480nnn/GSM1480321/suppl/GSM1480321_K562_GROcap_wTAP_plus.bigWig'
```