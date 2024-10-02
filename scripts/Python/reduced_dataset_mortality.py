from clock_utils import run_kfold_testing_mortality

with open('GTEx_4x_FC_genes_longitudinal.json', 'r') as file:
    GTEx_4x_FC_genes = json.load(file)
del GTEx_4x_FC_genes['Bladder']
del GTEx_4x_FC_genes['Thyroid']
        
feature_subsets = { "GTEx_4x_FC_longitudinal" : GTEx_4x_FC_genes}

paral_n = sys.argv[1]
mortality_results = run_kfold_testing_mortality("mortality_reduced/", ["Status", "Time"], "Age", feature_subsets, n_split=paral_n,
                                validation_folds=2,
                                save_coefs=True, file_prefix="mortality_", save_dir="mortality_reduced/")
