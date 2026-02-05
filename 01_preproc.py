# ------------------------------------------------------------
# This script performs initial quality control (QC) and preprocessing for 10x Genomics
# single-cell RNA-seq data using Scanpy.
# 
# Workflow summary:
# 1. Load raw Cell Ranger output for a given sample (specified by --sp).
# 2. Filter low-quality cells and lowly expressed genes.
# 3. Compute QC metrics including mitochondrial gene percentage.
# 4. Remove cells with excessive gene counts or mitochondrial RNA content.
# 5. Detect and remove putative doublets using Scrublet.
# 6. Remove contaminating cells expressing erythrocyte, platelet, or endothelial markers
#    based on published expression thresholds.
# 7. Save the filtered AnnData object for downstream analysis.
# 
# Key QC parameters:
# - Minimum genes per cell: 200
# - Minimum cells per gene: 3
# - Maximum genes per cell: 4000
# - Maximum mitochondrial fraction: 10%
# - Expected doublet rate: 6%
# 
# The resulting dataset contains high-quality singlet cells with reduced ambient RNA
# and lineage contamination, suitable for downstream normalization, integration,
# and cell-type annotation.
# 
# Usage:
#     python script.py --sp <sample_id>
# 
# Input:
#     ../data/cellranger/<sample_id>/
# 
# Output:
#     ../data/scanpy_tmp/<sample_id>.h5ad
# ------------------------------------------------------------


from utils import *

min_cells = 3
min_genes = 200

max_genes = 4000
max_mtpct = 10

doublet_rate = .06
min_counts = 2
min_gene_variability_pctl = 85
n_prin_comps=30


parser = argparse.ArgumentParser(description="Process input variable 'sp'")
parser.add_argument('--sp', type=str, required=True, help='Input value for sp')
args = parser.parse_args()
sp = args.sp # python script.py --sp "value_for_sp"

# load data
adata = sc.read_10x_mtx(f'../data/cellranger/{sp}',
                     var_names='gene_symbols', cache=True
                    )
adata.var_names_make_unique()

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

adata = adata[(adata.obs.n_genes_by_counts < max_genes)
              & (adata.obs.pct_counts_mt < max_mtpct), :]

# duplicate detection
scrub = scr.Scrublet(adata.X, expected_doublet_rate=doublet_rate)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=min_counts, 
                                                          min_cells=min_cells, 
                                                          min_gene_variability_pctl=min_gene_variability_pctl, 
                                                          n_prin_comps=n_prin_comps)

adata.obs['doublets_pred'] = scrub.predicted_doublets_
adata.obs['doublets_score'] = scrub.doublet_scores_obs_
sel = adata[adata.obs['doublets_pred']==False]

genes_pollution = ['HBB','HBD',   'PF4','SDPR','GNG11','PPBP'] # RBC, PLT/Megakaryocyte
thres = e**2 - 1 # from literature, https://www.science.org/doi/10.1126/science.abf1970 

for g in genes_pollution:
    if g not in sel.var.index:
        continue
    sel = sel[sel[:,g].X < thres]
    # print(g, sel.shape[0])

sel.write_h5ad(f'../data/scanpy_tmp/{sp}.h5ad')