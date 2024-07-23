import sys
import pandas as pd 
import anndata
import numpy as np
import argparse
import glob
import sys
from passenger.plot.plot import get_high_conf_cells


parser = argparse.ArgumentParser(
                    prog='select_and_annotate.py',
                    description='Script to select the best wNMF output of a directory and optionally annotate the clones based on cell types.',
                    epilog="The following script compares the anndata output of all wNMF runs to select the best one based on the orthogonality score. Optionally the factors of the best output is then annotated as healthy / cancer based on cell types. This second step is only done if the necessary cell type labels are either saved in the anndata object (under 'adata.obs.celltype') or provided as a '.csv' file and a reference healthy or cancer cell type is provided.")

parser.add_argument('-d', '--dir', required=True, help="Path to the directory containing all output anndata objects. All objects matching with 'out_adata_*.h5ad' will be compared.")
parser.add_argument('-m', '--min_vars', required=False, default=50, type=int, help="Optional argument defining minimum number of variants needed in the anndata object. In case of noisy data, the wNMF can inadvertedly find patterns that look like co-occurence if the number of variants is too low. Set this parameter to 0 if you only have validated variants in your data. (Default of 50.)")

parser.add_argument('-ct', '--cell_types', required=False, help="Optional argument with path to .csv file containing the cell type labels. We assume the first collumn to be the cell names as saved in the adata.obs_names, and the second column a list of cell types.")
parser.add_argument('-H', '--healthy_ct', required=False, help="Optional argument providing a list of healthy cell type(s).")
parser.add_argument('-C', '--cancer_ct', required=False, help="Optional argument providing a list of cancer cell type(s).")

args = parser.parse_args()
path, min_vars = args.dir, args.min_vars
cell_types, healthy_ct, cancer_ct = args.cell_types, args.healthy_ct, args.cancer_ct
files = glob.glob(path + '/out_adata_*.h5ad')

adata = None
best_s = np.infty
for file in files:
    adata_ = anndata.read_h5ad(file)
    s = adata_.uns["orth_score"]
    if (adata_.shape[1]>min_vars) & (s < best_s):
        adata, best_s = adata_, s
        
if adata is None:
    sys.exit()

if cell_types is not None:
    df_ct = pd.read_csv(cell_types, index_col=0, header=None)
    sub_ct = np.intersect1d(df_ct.index.tolist(), adata.obs_names.tolist())
    if len(sub_ct) == 0:
        print("Indices in "+cell_types+" do not match 'adata.obs_names'. Please double check, fix and rerun.")
        sys.exit()
    if len(sub_ct)<adata.shape[0]:
        print("Indices in "+cell_types+" do not fully match 'adata.obs_names'. We have "+ adata.shape[0]-len(sub_ct) + " cell names missing. Continuing.")
    adata = adata[sub_ct]
    print("Saving cell types to anndata object.")
    adata.obs["celltype"] = df_ct.loc[sub_ct][1]
    
threshold = .05
Cdiff_sub_cells = .3

def get_clone_ct_perc(ct, label, all_types, cell_factor_labels):
    
        ct = ct.split(",") 
        found = np.array([i in all_types for i in ct])
        if ~np.any(found):
            print("None of the provided "+label+" cell type labels are found in 'adata.obs['celltype']'. Skipping.")
            return None
        elif np.any(~found):
            print("The following cell types were not found in 'adata.obs['celltype']':" + ct[~found] + " Continuing with the ones that were found.")
            ct = ct[found]
            
        labeled_cells = np.array([i in ct for i in sub_adata.obs["celltype"]])
        if (np.sum(labeled_cells) / sub_adata.shape[0]) < .05 : 
            print("WARNING: we have less than 5% of "+label+" cell from cell types in the labeled cells. The wNMF and annotation might fail it the expected clonal populations are too small. Continuing.")
            
        n_cells_labeled = [np.sum(labeled_cells[(cell_factor_labels==i)]) for i in np.unique(cell_factor_labels)]
        n_cells_labeled /= np.sum(labeled_cells)
        return n_cells_labeled

if "celltype" in adata.obs.columns and ((healthy_ct is not None) or (cancer_ct is not None)):
    all_types = np.unique(adata.obs["celltype"])
    
    Cdiff = get_high_conf_cells(adata, Cdiff_sub_cells=Cdiff_sub_cells)
    sub_adata = adata[Cdiff]
    
    cell_factor_labels = np.argmax(sub_adata.obsm["C"], axis=1)
   
    if healthy_ct is not None:
        n_cells_healthy = get_clone_ct_perc(healthy_ct, "healthy", all_types, cell_factor_labels)
    if cancer_ct is not None:
        n_cells_cancer = get_clone_ct_perc(cancer_ct, "cancer", all_types, cell_factor_labels)
    
    labels = np.repeat("ambiguous", adata.uns["k"])
    
    if n_cells_cancer is not None and n_cells_healthy is not None:
        cancer_factors = (n_cells_cancer>threshold)
        healthy_factors = (n_cells_healthy>threshold)
        idx = cancer_factors & healthy_factors
        if np.any(idx):
                  print("Factors "+str(np.where(idx)[0])+ " contain both heahlty and cancer cell types. We assume the wNMF has failed and do not label the factors.")
        else:
            labels[cancer_factors]="cancer"
            labels[healthy_factors]="healthy"
            print("Labeling factors "+ str(np.where(cancer_factors)[0]) + " as cancer.")
            print("Labeling factors "+ str(np.where(healthy_factors)[0]) + " as healthy.")
        idx = np.any((~cancer_factors) & (~healthy_factors))
        if np.any(idx):
            print("Factor(s) "+ str(np.where(idx)[0]) + " contain neither healthy nor cancer cell types. Labeling as ambiguous.")
    elif n_cells_cancer is not None:
        cancer_factors = (n_cells_cancer>threshold)
        labels[cancer_factors]="cancer"
        labels[~cancer_factors]="healthy"
        print("Labeling factors "+ str(np.where(cancer_factors)[0]) + " as cancer.")
        print("Labeling factors "+ str(np.where(~cancer_factors)[0]) + " as healthy.")
    elif n_cells_healthy is not None:
        healthy_factors = (n_cells_healthy>threshold)
        labels[healthy_factors]="healthy"
        labels[~healthy_factors]="cancer"
        print("Labeling factors "+ str(np.where(~healthy_factors)[0]) + " as cancer.")
        print("Labeling factors "+ str(np.where(healthy_factors)[0]) + " as healthy.")

    if np.any(labels=="healthy") & np.any(labels=="cancer"): # SUCCESS
        print("wNMF validation and factor labeling SUCCESS. Saving the labels to adata object.")
        adata.uns["factor_labels"]=labels
        C = adata.obsm["C"]
        healthy = np.sum(C[:,labels == "healthy"], axis=1)
        cancer = np.sum(C[:,labels == "cancer"], axis=1)
        conf = np.abs(healthy-cancer)
        cell_labels = np.repeat("undetermined", adata.shape[0])
        cell_labels[(conf>Cdiff_sub_cells) & (cancer>healthy)] = "cancer"
        cell_labels[(conf>Cdiff_sub_cells) & (cancer<healthy)] = "healthy"
        print(cell_labels)
        adata.obs["cell_labels"]=cell_labels
        adata.obs["cancer_weight"]=cancer
        adata.obs["healthy_weight"]=healthy
    else:
        print("wNMF validation and factor labeling FAILED. The found factors likely do not correspond to genetic clones.")

print("Saving adata object to '"+path+"final.h5ad'.")
adata.write_h5ad(path+"final.h5ad")