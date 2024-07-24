# Scripts to run CCLONE

The scripts in this directory show how to run CCLONE in three main steps:
(1) preprocessing (variant calling with cellSNP-lite, and conversion to anndata format) 
(2) filtering and running wNMF
(3) [opt] selection of result, annotation and validation based on cell types

A step-by-step tutorial on how to run this analysis starting from step (2) and look at the result is provided in notebooks for two example patients in the ../tutorial/ directory.

## (1) Preprocess:
### (A) Run cellSNP-lite and annotate output with vep:
run_cellSNP-lite.sh

#### Requirements
- ```dir``` containing 'possorted_genome_bam.bam' bam file and 'barcodes.txt' list of barcodes of interest. We recommend that 'barcodes.txt' contains only the subset of barcodes passing the usual QC (excluding very low count barcodes / doublets, etc...), as this will significantly reduce the runtime of cellSNP-lite.
- [cellsnp-lite](https://cellsnp-lite.readthedocs.io/en/latest/) installed and runable via the command line (added to $PATH - or edit path in run_cellSNP-lite.sh)
- [vep](https://www.ensembl.org/info/docs/tools/vep/index.html) installed and runable via the command line
- data:
  - ```hg38.fa.masked``` fasta file of human reference genome file. In this work we used the repeat masked reference genome dowloaded from [ucsc](https://genome.ucsc.edu/index.html). You can download this file via ```wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.masked.gz```.
  - ```00-common_all-with_chr.vcf.gz``` vcf file of common human germline variants for annotation. In this work we downloaded all common variants from dbSNP with ```https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz```. You may need to convert the chromosome name format with: ```awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' 00-common_all.vcf > 00-common_all_chr.vcf``` (with decompression and recompression after conversion). Then create an index with ```bcftools index 00-common_all_chr.vcf```.
  - ```RNAedits.vcf.gz``` vcf file of known RNA edits. The converted file from [REDIportal](http://srv00.recas.ba.infn.it/atlas/) used in this work can be found in the RNAedit/ directory.
    
#### Running 
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

### (B) Convert cellSNP-lite output to anndata
cellSNP-lite_to_anndata.py

#### Requirements
- CCLONE package and requirements installed
- run_cellSNP-lite.sh successfully ran for all chromosomes of interest and output saved in ```dir```

#### Running 
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

## (2) Filter and run wNMF
run_wNMF.py

#### Requirements
- CCLONE package and requirements installed
- cellSNP-lite_to_anndata.py successfully ran and output adata.h5ad file saved in ```dir```

#### Running 
In it's simplest use and with default parameter, the python script can be run with:
```
python run_wNMF.py -d dir -k k
```
With ```dir``` the directory in which the adata.h5ad object is found, and in which the output will be save and ```k``` the number of clones. 

The variants can be filtered based on minimal minor allele frequency (MAF) and on annotation in dbSNP (=known germline variants). In this work we submitted this script for run in parrallel with a range of -k between 2 and 5, --min_MAF of 2 5 and 10 and with and without the -g flag (to filter germline variants), and then selected the result in the next step.

Note if your adata object contains only validated somatic variants, you might want to skip filtering of variants altogether by commenting out the corresponding line in the file, or by setting the thresholds very low, and only run the wNMF for a range of k.

To learn more about the run parameters, please run:
```
python run_wNMF.py -h
```
Editable run options are:
- the input number of factors k
- update min MAF (of 5%) under which variants are excluded from the anndata object
- update min coverage (of 10%) in variant positions under which cells are excluded from the anndata object
- deactivate filtering of known germline variants (with -g flag)
- update number of wNMF EM cycles (default of 100)

## (3) Selection of result, annotation and validation based on cell types
select_and_annotate.py

#### Requirements
- CCLONE package and requirements installed
- run_wNMF.py successfully ran and output file(s) saved in ```dir```

#### Running 
This script performs two tasks:
- first, looping through all output anndata object(s) found in ```dir``` and selecting the best result based on the orthogonality score
- then (optionally, if the necessary information is provided) validating and annotating the factors based on cell types. For that you need cell type labels (either provided as .csv or saved under adata.var.celltypes), as well as the information of either one (or more) healthy and/or one (or more) cancer cell type(s).

This script can be run with (with an example AML dataset in which the T cells are known healthy cells and Blasts are known cancer cells):
```
python select_and_annotate.py -d dir -ct celltypes.csv -H "T cells"
```
or 
```
python select_and_annotate.py -d dir -ct celltypes.csv -C "Blasts"
```
or 
```
python select_and_annotate.py -d dir -ct celltypes.csv -H "T cells" -C "Blasts"
```
With ```dir``` the directory in which the output anndata objects are found, and in which the output will be saved.

To learn more about the run parameters, please run:
```
python run_wNMF.py -h
```
Editable run options are:
- ```-ct```: Provide optional path to the list of cell type labels. We assume the first collumn to be the cell names as saved in the adata.obs_names, and the second column a list of cell types.
- ```-H```: Name(s) of known healthy cell type(s). Comma separated if multiple names provided.
- ```-C```: Name(s) of known cancer cell type(s). Comma separated if multiple names provided.
