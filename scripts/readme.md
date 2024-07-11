# Preprocess for CCLONE

These two scripts show how to run the variant calling with cellSNP-lite, and how to convert the output to anndata format.

## (1) Run cellSNP-lite and annotate output with vep:
run_cellSNP-lite.sh

### Requirements
- ```dir``` containing 'possorted_genome_bam.bam' bam file and 'barcodes.txt' list of barcodes of interest. We recommend that 'barcodes.txt' contains only the subset of barcodes passing the usual QC (excluding very low count barcodes / doublets, etc...), as this will significantly reduce the runtime of cellSNP-lite.
- [cellsnp-lite](https://cellsnp-lite.readthedocs.io/en/latest/) installed and runable via the command line (added to $PATH - or edit path in run_cellSNP-lite.sh)
- [vep](https://www.ensembl.org/info/docs/tools/vep/index.html) installed and runable via the command line
- data:
  - ```hg38.fa.masked``` fasta file of human reference genome file. In this work we used the repeat masked reference genome dowloaded from [ucsc](https://genome.ucsc.edu/index.html). You can download this file via ```wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.masked.gz```.
  - ```00-common_all-with_chr.vcf.gz``` vcf file of common human germline variants for annotation. In this work we downloaded all common variants from dbSNP with ```https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz```. You may need to convert the chromosome name format with: ```awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' 00-common_all.vcf > 00-common_all_chr.vcf``` (with decompression and recompression after conversion). Then create an index with ```bcftools index 00-common_all_chr.vcf```.
  - ```RNAedits.vcf.gz``` vcf file of known RNA edits. The converted file from [REDIportal](http://srv00.recas.ba.infn.it/atlas/) used in this work can be found in the RNAedit/ directory.
    
### Running 
The bash script can be submitted for each chromosome in parallel and is called for one chromosome (chr2) with:
```
bash run_cellSNP-lite.sh -d dir -c chr2 -h hg38.fa.masked -r RNAedits.vcf.gz -s 00-common_all-with_chr.vcf.gz
```
or running for all chromosomes:
```
for i in {1..22};
  do bash run_cellSNP-lite.sh -d path_to_directory_with_bam_file -c chr$i -h hg38.fa.masked -r RNAedits.vcf.gz -s 00-common_all-with_chr.vcf.gz;
done
```

## (2) Convert cellSNP-lite output to anndata
cellSNP-lite_to_anndata.py

### Requirements
- CCLONE package and requirements installed
- run_cellSNP-lite.sh successfully ran for all chromosomes of interest and output saved in ```dir```

### Running 
In it's simplest use and with default parameter, the python script can be run with:
```
python cellSNP-lite_to_anndata.py -d dir
```
To learn more about the run parameters, please run:
```
python cellSNP-lite_to_anndata.py -h
```
Editable run options are:
- running only for a subset of chromosomes
- update min coverage (of 10%) under which variants are excluded from the anndata object
- deactivate filtering of variants in HLA regions
- deactivate merging of highly correlated neighboring variant (to account for variants potentially found on the same reads, and who break the assumption of sampling independance)
