from clock_utils import run_kfold_testing_chronological

with open('feature_subsets/GTEx_4x_FC_genes.json', 'r') as file:
    GTEx_4x_FC_genes = json.load(file)
del GTEx_4x_FC_genes['Bladder']

feature_subsets = { "GTEx_4x_FC" : GTEx_4x_FC_genes}

paral_n = sys.argv[1]

results = run_kfold_testing_chronological("kfolds/test_kf/", "Age", feature_subsets, n_split=paral_n, save_dir="output/", validation_folds =5)