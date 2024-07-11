#!/bin/bash


while getopts d:s:r:h:c: flag
do
    case "${flag}" in
        d) dir=${OPTARG};;  # directory containing the possorted_genome_bam.bam and barcodes.txt file - output will be saved here as well
        c) chrom=${OPTARG};;
        s) dbfile=${OPTARG};;  # path to dbSNP vcf file - see readme
        r) rnaEditFile=${OPTARG};;  # path to REDIdb vcf file - downloadable from github
        h) ref=${OPTARG};;  # path to (repeat masked) human reference genome
    esac
done

cells=100 # this is the min number of cells with coverage - opt: set to ~10% of the total number of cells

echo "running var call"
cellsnp-lite -s $dir/possorted_genome_bam.bam -b $dir/barcodes.txt -O $dir/chr${SGE_TASK_ID} -p 22 --minMAF 0.02 --minCOUNT ${cells} --chrom ${chrom} --refseq ${ref}
echo "running vep"
vep -i $dir/chr${SGE_TASK_ID}/cellSNP.base.vcf --cache --dir_cache ~/fast/vep-cache/ --symbol --tab -o $dir/${chrom}/annotations.tsv --force_overwrite --check_existing --exclude_null_alleles --custom  ${rnaEditFile},REDI,vcf,exact --custom ${dbfile},dbSNP-common,vcf,exact --af --fields "Location,Allele,Gene,Feature,Feature_type,Consequence,Amino_acids,Codons,Existing_variation,IMPACT,REDI,AF,dbSNP-common"
