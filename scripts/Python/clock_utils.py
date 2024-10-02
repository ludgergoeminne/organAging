import os
import sys
import numpy as np
import pandas as pd
import pickle
import json
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import StratifiedKFold, KFold, GridSearchCV
from sksurv.linear_model import CoxnetSurvivalAnalysis
from lifelines import WeibullAFTFitter


def aft_model(data, time, status):
    model = WeibullAFTFitter(penalizer=0.01)
    model.fit(data, duration_col=time, event_col=status)
    return model


def elastic_grid_search(X, y, validation_folds=5, alpha_search=np.logspace(-4, -1, 4),
                        l1_ratio_search=np.linspace(0.01, 0.4, 3), n_bins=5, n_jobs=-1, random_state = 1716326954):

    binned_y = pd.cut(y, bins=n_bins, labels=False)
    kf = StratifiedKFold(n_splits=validation_folds)
    split = kf.split(binned_y, binned_y)

    param_grid = {'alpha': alpha_search,
                  'l1_ratio': l1_ratio_search,
                  }

    grid_model = GridSearchCV(ElasticNet(random_state=random_state, fit_intercept=True), param_grid, scoring="neg_mean_absolute_error",
                              cv=split, n_jobs=n_jobs)
    grid_model.fit(X, y)
    
    print("Best l1_ratio " + str(grid_model.best_params_['l1_ratio']), flush=True)
    return grid_model


def mortality_grid_search(X, y, validation_folds=3, l1_ratio_search=[0.5, 1.0], alpha_min_ratio="auto",
                          random_state=1716326956, n_jobs=-1):
    kf = KFold(n_splits=validation_folds, shuffle=True, random_state=random_state)

    y = y.to_records(index=False)
    coxnet = CoxnetSurvivalAnalysis(alpha_min_ratio=alpha_min_ratio, fit_baseline_model=True, n_alphas=50).fit(X, y)
    gcv = GridSearchCV(
        coxnet,
        {"alphas": [[v] for v in coxnet.alphas_],
         "l1_ratio": l1_ratio_search},
        cv=kf, error_score=0, n_jobs=n_jobs).fit(X, y)
    print("Best alpha " + str(gcv.best_params_['alphas']))
    print("Best l1_ratio " + str(gcv.best_params_['l1_ratio']))

    return (gcv)


def save_coefficients(feature_results, feature_subsets, file_prefix, save_dir, intercept=True):
    for feature_subset_name in feature_subsets.keys():
        for organ in feature_subsets[feature_subset_name].keys():
            colnames = list(feature_subsets[feature_subset_name][organ])
            if intercept:
                df = pd.concat([pd.DataFrame(feature_results[feature_subset_name][organ]['intercept']),
                                pd.DataFrame(feature_results[feature_subset_name][organ]['features_coefficients'])], axis=1)

                colnames.insert(0, "Intercept")
                df.columns = colnames
                df.to_csv(save_dir+organ+file_prefix + "_coefs_" + feature_subset_name + ".csv",
                          index=False)
            else:
                df = pd.DataFrame(feature_results[feature_subset_name][organ]['features_coefficients'])
                df.columns = colnames
                df.to_csv(save_dir + organ + file_prefix + "_coefs_" + feature_subset_name + ".csv",
                          index=False)


def create_empty_feature_results(feature_subsets):
    feature_results = dict()
    for feature_subset_name in feature_subsets.keys():
        feature_results[feature_subset_name] = dict()
        for organ in feature_subsets[feature_subset_name].keys():
            feature_results[feature_subset_name][organ] = dict()
            feature_results[feature_subset_name][organ]['features_coefficients'] = list()
            feature_results[feature_subset_name][organ]['intercept'] = list()
            feature_results[feature_subset_name][organ]['correlation'] = list()
    return feature_results


def check_features(data, feature_subsets):
    for feature_subset_name in feature_subsets.keys():
        for organ in feature_subsets[feature_subset_name].keys():
            # subset by availiable columns (in case wasn't done before and sort alphabetically
            feature_subsets[feature_subset_name][organ] = sorted(
                set(feature_subsets[feature_subset_name][organ]).intersection(data.columns))
    return feature_subsets


def run_kfold_testing_chronological(kfold_directory, age_column_name, feature_subsets, n_split,
                                validation_folds=5,
                                save_coefs=True, file_prefix="", save_dir="", n_jobs=-1):
    train_data = pd.read_csv(kfold_directory + "train_imputed_" + str(n_split) + ".csv", index_col=0)
    train_data.rename(columns=lambda x: x.replace('.', '-'), inplace=True)
    feature_subsets = check_features(train_data, feature_subsets)

    feature_results = create_empty_feature_results(feature_subsets)

    print("fold "+str(n_split), flush=True)
    

    y_train = train_data[age_column_name]

    for feature_subset_name in feature_subsets.keys():
        print(feature_subset_name, flush=True)
        for organ in feature_subsets[feature_subset_name].keys():
            print(organ, flush=True)
            current_feature_subsets = feature_subsets[feature_subset_name][organ]
            X_subset = train_data[current_feature_subsets]

            
            trained_model = elastic_grid_search(X_subset, y_train, validation_folds=validation_folds,
                                                      n_jobs=n_jobs, l1_ratio_search=[0.0, 0.2, 0.4, 0.8, 1.0])
            if save_coefs:
                feature_results[feature_subset_name][organ]['features_coefficients'].append(trained_model.best_estimator_.coef_)
                feature_results[feature_subset_name][organ]['intercept'].append(trained_model.best_estimator_.intercept_)

    if save_coefs:
        save_coefficients(feature_results, feature_subsets, file_prefix+str(n_split)+"_", save_dir, intercept=True)
    return feature_results



def run_kfold_testing_mortality(kfold_directory, event_column_name, age_column_name, feature_subsets, n_split,
                                validation_folds=2,
                                save_coefs=True, file_prefix="mortality_", save_dir="", model="aft", n_jobs=-1):
    train_data = pd.read_csv(kfold_directory + "train_imputed_" + str(n_split) + ".csv", index_col=0)
    train_data['Status'] = train_data['Status'].astype(bool)
    
    train_data.rename(columns=lambda x: x.replace('.', '-'), inplace=True)
    feature_subsets = check_features(train_data, feature_subsets)

    feature_results = create_empty_feature_results(feature_subsets)

    print("fold "+str(n_split), flush=True)
    y_train = train_data[event_column_name]

    for feature_subset_name in feature_subsets.keys():
        print(feature_subset_name, flush=True)
        for organ in feature_subsets[feature_subset_name].keys():
            print(organ, flush=True)
            current_feature_subsets = feature_subsets[feature_subset_name][organ]
            X_subset = train_data[current_feature_subsets]

            if model == "elastic":
                trained_model = mortality_grid_search(X_subset, y_train, validation_folds=validation_folds,
                                                      n_jobs=n_jobs)
                feature_results[feature_subset_name][organ]['features_coefficients'].append(trained_model.best_estimator_.coef_.flatten())
            elif model == "aft":
                data = pd.concat([X_subset, y_train], axis=1)
                trained_model = aft_model(data, "Time", "Status")
                feature_subsets[feature_subset_name][organ] = trained_model.summary["coef"].index.get_level_values("covariate")[:-1]
                feature_results[feature_subset_name][organ]['features_coefficients'].append(trained_model.summary["coef"].values[:-1])
            else:
                raise AttributeError("Model must be elastic. Will be extended to include svm later.")
                

    if save_coefs:
        save_coefficients(feature_results, feature_subsets, file_prefix+str(n_split)+"_", save_dir, intercept=False)
    return feature_results