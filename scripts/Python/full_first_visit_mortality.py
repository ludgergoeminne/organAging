from clock_utils import run_kfold_testing_mortality

with open('feature_subsets/GTEx_4x_FC_genes.json', 'r') as file:
    GTEx_4x_FC_genes = json.load(file)
del GTEx_4x_FC_genes['Bladder']

feature_subsets = { "GTEx_4x_FC" : GTEx_4x_FC_genes}

paral_n = sys.argv[1]

mortality_results = run_kfold_testing_mortality("kfolds/test_kf/", ["Status", "Time"], "Age", feature_subsets, n_split=paral_n,
                                validation_folds=2,
                                save_coefs=True, file_prefix="_mortality_", save_dir="output/")