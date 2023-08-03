import scvi

data_dir = "/projects/b1038/Pulmonary/cpuritz/PASC/data/NEP_exon_only/models"
model_adata = scvi.model.SCVI.load(data_dir)
model_adata.train(max_epochs = 500, use_gpu = True, check_val_every_n_epoch = 2, early_stopping = True)
model_adata.save(data_dir, overwrite = True, prefix = "trained")
