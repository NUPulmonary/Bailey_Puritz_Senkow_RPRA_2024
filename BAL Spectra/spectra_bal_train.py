import scanpy as sc
from spectra import spectra as spc
import pickle

data_dir = "/projects/b1038/Pulmonary/cpuritz/PASC/data/01BAL/spectra"
adata = sc.read_h5ad(f"{data_dir}/adata.h5ad")
with open(f"{data_dir}/input_gene_sets.pkl", 'rb') as f:
    input_gene_sets = pickle.load(f)
with open(f"{data_dir}/L.pkl", 'rb') as f:
	L = pickle.load(f)
model = spc.est_spectra(adata = adata,
			gene_set_dictionary = input_gene_sets,
			L = L,
			use_highly_variable = True,
			cell_type_key = "cell_type_spectra",
			use_weights = True,
			lam = 0.01,
			delta = 0.001,
			use_cell_types = True,
			n_top_vals = 50,
			num_epochs = 10000,
			verbose = True,
			rho = None,
			kappa = None)
with open(f"{data_dir}/model.pkl", 'wb') as f:
	pickle.dump(model, f)
adata.write(f"{data_dir}/adata.h5ad")
