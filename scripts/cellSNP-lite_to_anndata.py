# call with:
# python preproc.py -d path_to_cellsnp-lite-output_dir

import sys
import pandas as pd 
import anndata
import argparse
from passenger.preprocess.filter import *
from passenger.preprocess.import_data import *

parser = argparse.ArgumentParser(
                    prog='cellSNP_to_anndata',
                    description='Script to convert cellsnp-lite output to anndata object',
                    epilog='The following script takes the output of the wrapper to run cellSNP-lite (run_cellSNP-lite.sh) and coverts the output to an anndata object.')

parser.add_argument('-d', '--dir', required=True, help="Path to the directory containing the cellSNP-lite output. The output anndata object will also be saved here.")
parser.add_argument('-c', '--all_chroms', required=False, nargs="*", default="None", help="Optional list of chromosome names to use. (Default of all chromosomes from chr1 to chr23.)")
parser.add_argument('-p', '--perc', required=False, default=10, type=int, help="Optional argument defining minimum percentage of cells a variant needs to be covered in for it to be kept in the anndata object. (Default of 10.)")
parser.add_argument('-f', '--filter_hela', required=False, default=True, type=bool, help="Optional argument defining whether all variants found in the HELA region in chr6 will be filtered out. (Default to True.)")

args = parser.parse_args()
path, all_chroms, perc, filter_hela = args.dir, args.all_chroms, args.perc, args.filter_hela

cell_names = pd.read_csv(path+"/chr1/cellSNP.samples.tsv", header=None)[0].tolist()
all_chroms = ["chr" + str(i) for i in range(1, 23)] # if all_chroms is None else all_chroms

REF, ALT, meta = get_variant_measurement_data(path, all_chroms, cell_names,
                                              perc, filter_hela)
REF, ALT, meta = filter_vars_from_same_read(REF, ALT, meta, dist=np.infty)
meta.pos = meta.pos.astype(str)

adata = anndata.AnnData(X=(REF+ALT).T, layers={"REF":REF.T, "ALT":ALT.T}, var=meta)
adata.write(filename=path+"/adata.h5ad")
