import sys
import pandas as pd 
import anndata
import argparse
from passenger.preprocess.filter import *
from passenger.preprocess.import_data import *
from passenger.cluster.NMF import *

parser = argparse.ArgumentParser(
                    prog='filer_and_run_wNMF',
                    description='Script to run the wNMF for a given anndata object',
                    epilog='The following script takes the an anndata object with necessary REF and ALT matrices in \'layers\', optionally filters the variants and runs the wNMF.')

parser.add_argument('-d', '--dir', required=True, help="Path to the directory containing the 'adata.h5ad' anndata object. Output will also be saved here.")
parser.add_argument('-k', '--n_factors', required=True, type=int, help="Argument defining number of factors.")

parser.add_argument('-g', '--filter_germline', required=False, default=True, type=bool, help="Optional argument defining whether to filter germline variants annotated with True unter 'adata.var.dbSNP'. (Default of True.)")
parser.add_argument('-m', '--min_MAF', required=False, default=5, type=int, help="Optional argument defining minimum minor allele frequency of a variant for it to be kept in the anndata object. (Default of 5.)")
parser.add_argument('-c', '--min_cell_cov', required=False, default=10, type=int, help="Optional argument defining minimum coverage in variant position of a cell for it to be kept in the anndata object. (Default to True.)")
parser.add_argument('-n', '--n_cycles', required=False, default=100, type=int, help="Optional argument defining number of EM cycles of the wNMF to run. (Default of 100.)")

args = parser.parse_args()
path = args.dir
filter_germline, min_MAF, min_cell_cov = args.filter_germline, args.min_MAF, args.min_cell_cov
k, max_cycles = args.n_factors, args.n_cycles

print("Loading adata object.")
adata = anndata.read_h5ad(path+"/adata.h5ad")
print("\tobject has shape: " + str(adata.shape))

print("Filtering anndata object.")
adata = filter_vars(adata, filter_germline=filter_germline, min_MAF=min_MAF, min_cell_cov=min_cell_cov)
print("\tafter filtering object has shape: " + str(adata.shape))

print("Running weighted NMF.")
bootstrap_wNMF(adata, k=k, max_cycles=max_cycles, # n_bootsstrap=2 ## for quick debugging
              )

print("Saving output.")
filename = path + "/adata"
filename += "_filter-" + str(filter_germline)+"-"+str(min_MAF)+"-"+str(min_cell_cov)
filename += "_run-" + str(k)+"-"+str(max_cycles)
filename += ".h5ad"

adata.write_h5ad(filename)