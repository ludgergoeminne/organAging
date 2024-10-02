from clock_utils import run_kfold_testing_mortality

with open('GTEx_4x_FC_genes.json', 'r') as file:
    GTEx_4x_FC_genes = json.load(file)
del GTEx_4x_FC_genes['Bladder']
del GTEx_4x_FC_genes['Organismal']
del GTEx_4x_FC_genes['Conventional']
del GTEx_4x_FC_genes['Multi-organ']

        
feature_subsets = { "GTEx_4x_FC" : GTEx_4x_FC_genes}

paral_n = sys.argv[1]
mortality_results = run_kfold_testing_mortality("kfolds/test_kf/", ["Status", "Time"], "Age", feature_subsets, n_split=paral_n,
                                save_coefs=True, file_prefix="_aft_", save_dir="output/", model="aft")
