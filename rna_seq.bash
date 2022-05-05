#anno ##optional
#projectname=$()
#samplename=$()

#Create envf and Intallation
##wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
##bash Miniconda3-latest-Linux-x86_64.sh
##conda create -n rna-seq
##conda activate rna-seq
##conda config --add channels bioconda
##conda config --add channels conda-forge
##conda install -c bioconda samtools rseqc fastqc htseq bioconductor-deseq2 stringtie 
##conda install trim-galore STAR r

#genome files and build index
##mkdir anno
##cd anno
##mkdir mm10
##mkdir hg19
##cd mm10
##wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
##tar zxvf Mus_musculus_UCSC_mm10.tar.gz
##STAR --runThreadN 30 --runMode genomeGenerate --genomeDir star_index_mm10 --genomeFastaFiles /Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /Annotation/Archives/archive-current/Genes/genes.gtf 
##cd ..
##wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
##tar zxvf hg19.fa.gz
##STAR --runThreadN 30 --runMode genomeGenerate --genomeDir star_index_hg19 --genomeFastaFiles /Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /Annotation/Archives/archive-current/Genes/genes.gtf 



#Activate the environment
conda activate rna-seq
mkdir $projectname
cd $projectname

#Download the public data
mkdir rawdata
for ((i = 12;i<=15;i++)); #
do
fastq-dump --split-3 -O data/ SRR0020$i.sra.sra 
done

#FastQC
cd ../$projectname
mkdir fastQC
cd fastQC
ls ../$projectname/rawdata/*gz |xargs -I [] echo 'nohup fastqc -o ./fastQC/ [] &'>fastqc.sh
bash fastqc.sh

## Adapter Trim[OPTIONAL]
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired ../$projectname/rawdata/$fq1 ../$projectname/rawdata/$fq2  --gzip -o $input_data

#Mapping
STAR --runThreadN 10 --genomeDir ../anno/mm10/star_index_mm10 --readFilesCommand zcat --readFilesIn R1.fastq.gz R2.fastq.gz --outFileNamePrefix /home1/liziyi/EPS/RNA-seq/su_mRNA-seq/mapping/C1-18-TdEPS-1   --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5 &

#Visualization the mapping ratio by R
Rscript ~rcode/mapping.R

#RSeQC
mkdirc rseqc
geneBody_coverage.py -r hg19.housekeeping.bed -i test1.bam,test2.bam,test3.bam  -o /rseqc/geneBody_coverage

#ht-seq and counts2rpkm_tpm
htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m intersection-nonempty yourfile_name.bam ~/reference/hisat2_reference/Homo_sapiens.GRCh38.86.chr_patch_hapl_scaff.gtf > $samplename.counts.txt
Rscript ~rcode/counts2rpkm_tpm.R

#deseq2
Rscript ~rcode/deseq2.R

#downstream
Rscript ~rcode/KEGG.GO.R
Rscript ~rcode/complexHeatmap.R





