import numpy as np
import scipy.io as sio
import pandas as pd
import h5py
import joblib
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.inspection import permutation_importance

import skexplain
from skexplain import ExplainToolkit
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix,
)
from sklearn.inspection import permutation_importance
import shap

from joblib import Parallel, delayed

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
print("Working directory:", os.getcwd())

# from skexplain.main import ExplainToolkit

# 
print("Hello, Python is working!")

# SHAP analysis for Xe data (stage 3)

# ============================================================
# 1-1. Load dataset from .mat file
# ============================================================

# Load the .mat file
# mat = sio.loadmat("/Volumes/LH - drive/research/Xe-code/Xe/out_Xe_lhs_new_3e7_test_fixK_U_MC_lowerT_alot_Rac1100_new_lowf.mat")

mat_path = "/Volumes/LH - drive/research/Xe-code/Xe/out_Xe_lhs_new_3e7_test_fixK_U_MC_lowerT_alot_Rac1100_new_lowf.mat"

f = h5py.File(mat_path, 'r')
print(list(f.keys()))

# Option A (recommended): if the .mat file directly contains a numeric variable named 'results'
# results = mat["results"]

# Option B: if the .mat file contains a struct named 'results_Xe_struct'
# and the numeric array is stored in the field 'results',
# then you can access it like this:
# results_struct = mat["results_Xe_struct"]          # this is a numpy structured array
# results = results_struct["results"][0, 0]          # field access; shape = (size_Xe, nparam_Xe)

# Access struct field: results_Xe_struct/results

results = f["results"][:]

print("results shape:", results.shape)

n_params,n_samples = results.shape


# ======================================================================
# 1-2. Define col_names exactly as in MATLAB script (order matters!)
# ======================================================================

col_names = [
    'ts_cc','kappa_gcc','kappa_rcc','Rs','Rp',
    'Ubse','Qc_pd', 'dQc', 'Ti', 'eta_ref',
    'Xm_init','Ms_init','Fr_w',
    'Xe_init','Fr_Xe',
    'Fd_mor','Fd_p',
    'Er_Pmor_pd',
    'Er_xec','Er_xe128v130','Er_xe128v132','Er_xe130v132',
    'Er_xe131v132','Er_xe134v132','Er_xe136v132'
]

# Ensure consistency between col_names and results matrix size
assert len(col_names) == n_params, (
    f"col_names length ({len(col_names)}) does not match results columns ({n_params})"
)

# ======================================================================
# 1-3. Select 12 RF input parameters (X_rf) matching MATLAB indexing
# ======================================================================

rf_feature_indices = list(range(0, 7)) + list(range(8, 10)) + list(range(14, 17))
# X_rf = results[:, rf_feature_indices]
X_rf = results[rf_feature_indices,:]

print("X_rf shape:", X_rf.shape)
print("Selected RF feature names:", [col_names[i] for i in rf_feature_indices])


# ======================================================================
# 1-4. Define parameter names used for labeling 
# ======================================================================

param_names = [
    'ts','kappa_g','kappa_r','Rs','Rp',
    'Ubse','Qc_pd', 'Ti', 'eta_ref',
    'Fr_Xe',
    'Fd_mor','Fd_p'
]

# ======================================================================
# 1-5. Identify indices of residual error columns (6 outputs)
#    Direct Python equivalent of MATLAB find(strcmp(...))
# ======================================================================

def find_col_index(name_list, target):
    """Return index of a column name (equivalent to MATLAB find(strcmp()))"""
    return name_list.index(target)

Erxec_idx       = find_col_index(col_names, 'Er_xec')
Erxe128v130_idx = find_col_index(col_names, 'Er_xe128v130')
Erxe128v132_idx = find_col_index(col_names, 'Er_xe128v132')
Erxe130v132_idx = find_col_index(col_names, 'Er_xe130v132')
Erxe131v132_idx = find_col_index(col_names, 'Er_xe131v132')
Erxe134v132_idx = find_col_index(col_names, 'Er_xe134v132')
Erxe136v132_idx = find_col_index(col_names, 'Er_xe136v132')


resid_col_ids = [
    Erxec_idx, 
    Erxe128v130_idx, 
    Erxe130v132_idx, 
    Erxe131v132_idx, 
    Erxe134v132_idx, 
    Erxe136v132_idx
]

resid_col_names = [
    '130Xe', '128Xe/130Xe', '130Xe/132Xe',
    '131Xe/132Xe', '134Xe/132Xe', '136Xe/132Xe'
]

print("Residual output column indices:", resid_col_ids)
print("Residual output names:", resid_col_names)


# ======================================================================
# 1-6. Construct Y_resid matrix of outputs
# ======================================================================

Y_resid = results[resid_col_ids,:]
print("Y_resid shape:", Y_resid.shape)

# # Define success based on RMSE threshold
# threshold = 1
# Y = Y_resid < threshold;  # success = 1, failure = 0

# ======================================================================
# 1-7. Convert X & Y to fit the format of sklearn
# ======================================================================

X_rf = X_rf.T  # 
Y_resid= Y_resid.T

# ubse_from_results = results[ find_col_index(col_names, 'Ubse'), : ]
# ubse_in_X = X_rf[:, 5]
# print("max abs diff:", np.max(np.abs((ubse_from_results - ubse_in_X)/ubse_from_results)))

from sklearn.preprocessing import StandardScaler

sc = StandardScaler()
X_scaled = sc.fit_transform(X_rf)

# # ============================================================
# # 2. Train–test split
# # ============================================================

# X_train, X_test, y_train, y_test = train_test_split(
#     X_rf, Y_resid, test_size=0.2, random_state=42
# )


# # ============================================================
# # 3. Train baseline Random Forest Regressor
# # ============================================================

# # RandomForestRegressor supports multi-output directly.
# # rf = RandomForestRegressor(
# #     n_estimators=200,     # number of trees
# #     max_depth=None,
# #     n_jobs=-1,            # use all CPU cores
# #     max_samples=500_000, # max smaples used by a tree
# #     random_state=42,
# # )
# # rf.fit(X_train, y_train)

# models = []
# for i in range(y_train.shape[1]):  # 6 output
#     print(f"training target {i}")
#     rf = RandomForestRegressor(
#         n_estimators=150,
#         max_depth=18,
#         min_samples_leaf=10,
#         max_samples=200_000,
#         n_jobs=-1,
#         random_state=0
#     )
#     rf.fit(X_train, y_train[:, i])
#     models.append(rf)

# # save
# joblib.dump(models, "rf_models_all_outputs.joblib")
# # # models = joblib.load("rf_models_all_outputs.joblib")


# # ============================================================
# # 4. Evaluate baseline model (MSE / R²)
# # ============================================================

# # y_pred = rf.predict(X_test)

# y_pred_list = []

# for rf in models:
#     y_pred_list.append(rf.predict(X_test))

# # 
# y_pred = np.column_stack(y_pred_list)

# # Raw values = one MSE per output dimension
# mse_per_output = mean_squared_error(y_test, y_pred, multioutput="raw_values")
# r2_per_output  = r2_score(y_test, y_pred, multioutput="raw_values")

# # Uniform average gives a single value summarizing performance
# mse_mean = mean_squared_error(y_test, y_pred, multioutput="uniform_average")
# r2_mean  = r2_score(y_test, y_pred, multioutput="uniform_average")

# print("\n=== Baseline RF Performance ===")
# for i, col in enumerate(resid_col_names):
#     print(f"{col}: MSE={mse_per_output[i]:.4f}, R2={r2_per_output[i]:.4f}")
# print(f"Overall (mean across 6 outputs): MSE={mse_mean:.4f}, R2={r2_mean:.4f}")

# # # save
# # metrics_df = pd.DataFrame({
# #     "output_name": resid_col_names,   # e.g., ["resid1", ..., "resid6"]
# #     "mse": mse_per_output,
# #     "r2": r2_per_output
# # })
# # metrics_df.to_csv("rf_metrics_per_output.csv", index=False)
# # print("Saved: rf_metrics_per_output.csv")

# # overall_df = pd.DataFrame({
# #     "metric": ["mse_mean", "r2_mean"],
# #     "value": [mse_mean, r2_mean]
# # })
# # overall_df.to_csv("rf_metrics_overall.csv", index=False)
# # print("Saved: rf_metrics_overall.csv")



# # ============================================================
# # 5. Feature importance and top-k feature retraining
# # ============================================================

# # importances = rf.feature_importances_             # shape (12,)
# # indices = np.argsort(importances)[::-1]           # sorted descending

# # print("\n=== Feature Importances (Gini Importance) ===")
# # for rank, idx in enumerate(indices):
# #     print(f"Rank {rank+1}: {param_names[idx]} importance={importances[idx]:.4f}")



# for out_i, rf in enumerate(models):
#     importances = rf.feature_importances_            # shape (12,)
#     indices = np.argsort(importances)[::-1]          # ranking

#     print(f"\n=== Feature Importances (Impurity) for output {out_i} ===")
#     for rank, idx in enumerate(indices):
#         print(f"Rank {rank+1}: {param_names[idx]} importance={importances[idx]:.4f}")

#     # select out limited samples for permutation importance
#     n_subsample = min(50_000, X_test.shape[0])
#     rng = np.random.default_rng(42)  # 
#     sub_idx = rng.choice(X_test.shape[0], size=n_subsample, replace=False)

#     X_pi = pd.DataFrame(X_test[sub_idx], columns=param_names)
#     y_pi = y_test[sub_idx, out_i]   #

#     perm_importance = permutation_importance(
#     rf,         # e.g. model for output 0
#     X_pi,
#     y_pi,
#     n_repeats=10,
#     random_state=42
#     )
#     indices = np.argsort(perm_importance.importances_mean)[::-1]          # ranking 

#     print(f"\n=== Feature Importances (Permutation) for output {out_i} ===")
#     for rank, idx in enumerate(indices):
#         print(f"Rank {rank+1}: {param_names[idx]} importance={perm_importance.importances_mean[idx]:.4f}")


# # ============================================================
# # NNN. redo with classify tree
# # ============================================================

X_rf = X_scaled

threshold = 1.0
Y_bin = (np.abs(Y_resid) < threshold).astype(int)
print("Y_bin shape:", Y_bin.shape)

# train/test split
X_train_c, X_test_c, y_train_c, y_test_c = train_test_split(
    X_rf,
    Y_bin,
    test_size=0.2,
    random_state=42
)



# clf_models_list = []      # list of classifiers, one per output
# perm_results_list = []    # store permutation_importance results (optional)

# rng = np.random.default_rng(42)

# for i in range(Y_resid.shape[1]):
#     print(f"\n==============================")
#     print(f"=== Output {i}: {resid_col_names[i]} ===")

#     # Select binary labels for this output
#     y_train_i = y_train_c[:, i]    # shape: (n_train,)
#     y_test_i  = y_test_c[:, i]     # shape: (n_test,)

#     # Train classifier for this single output
#     clf = RandomForestClassifier(
#         n_estimators=200,
#         max_depth=None,
#         n_jobs=-1,
#         random_state=42,
#     )
#     clf.fit(X_train_c, y_train_i)
#     clf_models_list.append(clf)

#     # ----------------------------------------------------
#     # 4. Classification metrics for this output
#     # ----------------------------------------------------
#     y_pred_cls = clf.predict(X_test_c)

#     acc  = accuracy_score(y_test_i, y_pred_cls)
#     # zero_division=0 avoids warnings if one class is missing in predictions
#     prec = precision_score(y_test_i, y_pred_cls, zero_division=0)
#     rec  = recall_score(y_test_i, y_pred_cls, zero_division=0)
#     f1   = f1_score(y_test_i, y_pred_cls, zero_division=0)

#     print("\n=== Classification metrics ===")
#     print(f"Accuracy : {acc:.4f}")
#     print(f"Precision: {prec:.4f}")
#     print(f"Recall   : {rec:.4f}")
#     print(f"F1-score : {f1:.4f}")

#     cm = confusion_matrix(y_test_i, y_pred_cls)
#     print("Confusion matrix:")
#     print(cm)

#     # 5. Permutation importance for this output
#     #    (Use a subsample to avoid being too slow / OOM)
#     # ----------------------------------------------------
#     n_sub = min(10_000, X_test_c.shape[0])  # subsample size
#     sub_idx = rng.choice(X_test_c.shape[0], size=n_sub, replace=False)
#     X_pi = X_test_c[sub_idx]
#     y_pi = y_test_i[sub_idx]

#     perm = permutation_importance(
#         clf,
#         X_pi,
#         y_pi,
#         n_repeats=5,       # you can increase if it is fast enough
#         random_state=42,
#         n_jobs=1,          # safer on macOS; avoids kernel crashes
#         scoring="f1",      # or "accuracy"
#     )
#     perm_results_list.append(perm)

#     indices = np.argsort(perm.importances_mean)[::-1]

#     print("\n=== Permutation feature importance (F1-based) ===")
#     for rank, idx_f in enumerate(indices):
#         print(
#             f"Rank {rank+1}: {param_names[idx_f]} "
#             f"importance={perm.importances_mean[idx_f]:.6f}"
#         )


# # ----- save
# # ------------------------------------------------------------
# # 1. Collect all metrics into a structured list
# # ------------------------------------------------------------
# metrics_list = []

# for i in range(Y_resid.shape[1]):
#     metrics_list.append({
#         "output_index": i,
#         "output_name": resid_col_names[i],
#         "accuracy": float(accuracy_score(y_test_c[:, i],
#                                         clf_models_list[i].predict(X_test_c))),
#         "precision": float(precision_score(y_test_c[:, i],
#                                           clf_models_list[i].predict(X_test_c),
#                                           zero_division=0)),
#         "recall": float(recall_score(y_test_c[:, i],
#                                      clf_models_list[i].predict(X_test_c),
#                                      zero_division=0)),
#         "f1": float(f1_score(y_test_c[:, i],
#                              clf_models_list[i].predict(X_test_c),
#                              zero_division=0)),
#         "confusion_matrix": confusion_matrix(y_test_c[:, i],
#                                              clf_models_list[i].predict(X_test_c))
#     })

# # ------------------------------------------------------------
# # 2. Extract permutation importance results into a dict
# # ------------------------------------------------------------
# perm_list = []
# for i, perm in enumerate(perm_results_list):
#     perm_list.append({
#         "output_index": i,
#         "output_name": resid_col_names[i],
#         "importances_mean": perm.importances_mean,
#         "importances_std": perm.importances_std,
#         "sorted_indices": np.argsort(perm.importances_mean)[::-1]
#     })

# # ------------------------------------------------------------
# # 3. Create final dictionary bundle
# # ------------------------------------------------------------
# save_bundle = {
#     "mat_path": mat_path, 
#     "models": clf_models_list,        # list of 6 classifiers
#     "metrics": metrics_list,          # list of dicts
#     "permutation_importance": perm_list,

#     "param_names": param_names,       # feature names
#     "output_names": resid_col_names,  # 6 outputs

#     "threshold": threshold,           # classification threshold used
# #    "scaler": scaler if "scaler" in globals() else None,  # store scaler if used
# }

# # ------------------------------------------------------------
# # 4. Save with joblib
# # ------------------------------------------------------------
# joblib.dump(save_bundle, "rf_classification_results.pkl")

# print("\nSaved all results to rf_classification_results.pkl")

# ------------------------------------------------------------
# reload the models
# -------------------------------------------------------------
# ------------------------------------------------------------
# 1. Load the saved bundle
# ------------------------------------------------------------
save_bundle = joblib.load("rf_classification_results.pkl")
print("Loaded rf_classification_results.pkl")

# ------------------------------------------------------------
# 2. Unpack the contents into variables (same names as before)
# ------------------------------------------------------------
mat_path        = save_bundle["mat_path"]
clf_models_list = save_bundle["models"]                  # list of 6 classifiers
metrics_list    = save_bundle["metrics"]                 # list of dicts
perm_list       = save_bundle["permutation_importance"]  # permutation importance results

param_names     = save_bundle["param_names"]             # feature names
resid_col_names = save_bundle["output_names"]            # 6 outputs

threshold       = save_bundle["threshold"]

# 3. Prepare data for MATLAB .mat file
mat_data = {}


# Basic info
mat_data['mat_path'] = mat_path
mat_data['param_names'] = np.array(param_names, dtype=object)
mat_data['output_names'] = np.array(resid_col_names, dtype=object)
mat_data['threshold'] = threshold

# Metrics (convert list of dicts to structured format)
for i, metrics in enumerate(metrics_list):
    prefix = f'model_{i+1}_'
    for key, value in metrics.items():
        if isinstance(value, (int, float, np.number)):
            mat_data[prefix + key] = value
        elif isinstance(value, np.ndarray):
            mat_data[prefix + key] = value
        elif isinstance(value, dict):
            # For nested dicts like classification_report
            for subkey, subvalue in value.items():
                if isinstance(subvalue, (int, float, np.number)):
                    mat_data[f'{prefix}{key}_{subkey}'] = subvalue
                elif isinstance(subvalue, dict):
                    for subsubkey, subsubvalue in subvalue.items():
                        mat_data[f'{prefix}{key}_{subkey}_{subsubkey}'] = subsubvalue

# Permutation importance
# if perm_list is not None:
#     for i, perm in enumerate(perm_list):
#         if hasattr(perm, 'importances_mean'):
#             mat_data[f'perm_importance_mean_{i+1}'] = perm.importances_mean
#             mat_data[f'perm_importance_std_{i+1}'] = perm.importances_std

n_targets = len(perm_list)
n_features = len(param_names)

RF_importance = np.zeros((n_targets, n_features))

for i, perm_result in enumerate(perm_list):
    # perm_result.importances_mean
        RF_importance[i, :] = perm_result["importances_mean"]

mat_data['RF_importance'] = RF_importance

# Model parameters (save important model info, not the model itself)
for i, clf in enumerate(clf_models_list):
    prefix = f'model_{i+1}_'
    # Save feature importances
    if hasattr(clf, 'feature_importances_'):
        mat_data[prefix + 'feature_importances'] = clf.feature_importances_
    
    # Save model parameters
    params = clf.get_params()
    for key, value in params.items():
        if isinstance(value, (int, float, bool)):
            mat_data[f'{prefix}param_{key}'] = value
        elif isinstance(value, str):
            mat_data[f'{prefix}param_{key}'] = value
    
    # Save number of estimators info
    if hasattr(clf, 'n_estimators'):
        mat_data[prefix + 'n_estimators'] = clf.n_estimators
    if hasattr(clf, 'n_features_in_'):
        mat_data[prefix + 'n_features'] = clf.n_features_in_

# 4. Save to .mat file
sio.savemat('rf_classification_results.mat', mat_data, oned_as='column')
print("Saved to rf_classification_results.mat")

# Print what was saved
print("\nSaved variables:")
for key in sorted(mat_data.keys()):
    if isinstance(mat_data[key], np.ndarray):
        print(f"  {key}: array shape {mat_data[key].shape}")
    else:
        print(f"  {key}: {type(mat_data[key]).__name__}")




# ============================================================
# 6. SHAP interaction values for all outputs (one clf per output)
#    NOTE: this is expensive; we use a small subset for each clf.
# ============================================================


print("\n=== SHAP interaction values for each output (class = 1) ===")

shap_interactions_list = []   # optional: store mean |interaction| per output
rng_shap = np.random.default_rng(42)

# choose how many samples to use for SHAP --> cluster irginal to n_shap seeds
n_shap = min(100, X_train_c.shape[0])

X_summary = shap.kmeans(X_train_c, n_shap).data
X_shap = X_summary

def compute_inter_for_one_model(out_i, clf_i, X_shap, n_shap):
    # # random seed
    # rng = np.random.default_rng(42 + out_i) 
    # idx_shap = rng.choice(X_train_c.shape[0], size=n_shap, replace=False)
    # X_shap = X_train_c[idx_shap]
    
    # 
    explainer = shap.TreeExplainer(clf_i)
    # 2.
    shap_inter = explainer.shap_interaction_values(X_shap)
    
    # select out Class 1 (match observations)
    if isinstance(shap_inter, list):
        # take interaction values for the positive class (label 1)
        shap_inter_pos = np.array(shap_inter[1])
    else:
        if shap_inter.ndim == 4:
            shap_inter_pos = shap_inter[..., 1]
        else:
            shap_inter_pos = shap_inter

    return np.abs(shap_inter_pos).mean(axis=0)

shap_interactions_list = Parallel(n_jobs=4)(
    delayed(compute_inter_for_one_model)(i, clf, X_shap, n_shap)
    for i, clf in enumerate(clf_models_list)
)


# ============ display ==================

for out_i, shap_inter_pos_mean_i in enumerate(shap_interactions_list):
    mean_abs_inter = shap_inter_pos_mean_i  # (n_features, n_features)

    # Optionally print top few strongest pairwise interactions
    # (flatten upper triangle and sort)
    n_features = mean_abs_inter.shape[0]
    upper_idx = np.triu_indices(n_features, k=1)
    inter_values = mean_abs_inter[upper_idx]
    sorted_idx = np.argsort(inter_values)[::-1]

    print("Top 5 pairwise interactions (by mean | SHAP interaction|):")
    for k in range(5):
        if k >= len(sorted_idx):
            break
        flat_idx = sorted_idx[k]
        i_feat = upper_idx[0][flat_idx]
        j_feat = upper_idx[1][flat_idx]
        # i_feat = int(sorted_idx[k,0])
        # j_feat = int(sorted_idx[k,1])
        # val=float(mean_abs_inter[i_feat,j_feat,1])
        print(
            f"  {param_names[i_feat]} x {param_names[j_feat]}: "
            f"{inter_values[flat_idx]:.4e}"
        ) 

results_dict = {}
for i, mean_abs_inter in enumerate(shap_interactions_list):
    # 
    results_dict[resid_col_names[i]] = mean_abs_inter

output_filename = "all_mean_abs_interactions.npz"
np.savez_compressed(output_filename, **results_dict)

print(f"\n all {len(shap_interactions_list)} shaps are saved to {output_filename}")


### ==================
# reload
npz_filename = "all_mean_abs_interactions.npz"
loaded_data = np.load(npz_filename)

print(f"Loaded {npz_filename}")
print(f"Keys in .npz file: {loaded_data.files}")

# 
param_names = [
    'ts','kappa_g','kappa_r','Rs','Rp',
    'Ubse','Qc_pd', 'Ti', 'eta_ref',
    'Fr_Xe',
    'Fd_mor','Fd_p'
]

resid_col_names = [
    '130Xe', '128Xe/130Xe', '130Xe/132Xe',
    '131Xe/132Xe', '134Xe/132Xe', '136Xe/132Xe'
]

n_targets = len(resid_col_names)
n_features = len(param_names)

print(f"Loaded {n_targets} targets, {n_features} features")


#
SHAP_interactions = np.zeros((n_targets, n_features, n_features))
for i, target_name in enumerate(resid_col_names):
    SHAP_interactions[i, :, :] = loaded_data[target_name]

#  Prepare data for MATLAB .mat file
mat_data = {
    # (n_targets, n_features, n_features)
    'SHAP_interactions': SHAP_interactions,
    
    # 
    'feature_names': param_names,
    'target_names': resid_col_names,
    
    
    'n_targets': n_targets,
    'n_features': n_features,
}

sio.savemat("SHAP_interactions.mat", mat_data)
print(f"\n✓ Saved to SHAP_interactions.mat")

loaded_data.close()



# for out_i, clf_i in enumerate(clf_models_list):
#     print("\n------------------------------")
#     print(f"Output {out_i}: {resid_col_names[out_i]}")

#     # subset of training data for SHAP
#     idx_shap = rng_shap.choice(X_train_c.shape[0], size=n_shap, replace=False)
#     X_shap = X_train_c[idx_shap]

#     # build TreeExplainer for this classifier
#     explainer_i = shap.TreeExplainer(clf_i)

#     # for binary classifiers, shap_inter_i is a [12*12] list of length 2 (class 0, class 1)
#     shap_inter_i = explainer_i.shap_interaction_values(X_shap)
#     # shap_inter_i = explainer_i.shap_interaction_values(X_shap, check_additivity=False)
   
#     if isinstance(shap_inter_i, list):
#         # take interaction values for the positive class (label 1)
#         shap_inter_pos = np.array(shap_inter_i[1])
#     else:
#         # some versions return a single array already
#        if shap_inter.ndim == 4:
        #     shap_inter_pos = shap_inter[..., 1]
        # else:
        #     shap_inter_pos = shap_inter

#     print("SHAP interaction array shape (class 1):", shap_inter_pos.shape)
#     # shap_inter_pos shape: (n_shap, n_features, n_features)

#     # Example: mean absolute interaction matrix for this output
#     mean_abs_inter = np.abs(shap_inter_pos).mean(axis=0)  # (n_features, n_features)
#     shap_interactions_list.append(mean_abs_inter)

#     # Optionally print top few strongest pairwise interactions
#     # (flatten upper triangle and sort)
#     n_features = mean_abs_inter.shape[0]
#     upper_idx = np.triu_indices(n_features, k=1)
#     inter_values = mean_abs_inter[upper_idx]
#     sorted_idx = np.argsort(inter_values)[::-1]

#     print("Top 5 pairwise interactions (by mean | SHAP interaction|):")
#     for k in range(5):
#         if k >= len(sorted_idx):
#             break
#         flat_idx = sorted_idx[k]
#         i_feat = upper_idx[0][flat_idx]
#         j_feat = upper_idx[1][flat_idx]
#         # i_feat = int(sorted_idx[k,0])
#         # j_feat = int(sorted_idx[k,1])
#         # val=float(mean_abs_inter[i_feat,j_feat,1])
#         print(
#             f"  {param_names[i_feat]} x {param_names[j_feat]}: "
#             f"{inter_values[flat_idx]:.4e}"
#         )



# ------------------------------------------
# # Select top-k most important features
# k = 6     # TODO: choose the number of features you want to keep
# topk_idx = indices[:k]
# topk_features = [param_names[i] for i in topk_idx]

# print(f"\nUsing top-{k} features:", topk_features)

# X_train_topk = X_train[:, topk_idx]
# X_test_topk  = X_test[:, topk_idx]

# rf_topk = RandomForestRegressor(
#     n_estimators=200,
#     max_depth=None,
#     n_jobs=-1,
#     random_state=42,
# )

# rf_topk.fit(X_train_topk, y_train)
# y_pred_topk = rf_topk.predict(X_test_topk)

# mse_per_output_topk = mean_squared_error(y_test, y_pred_topk, multioutput="raw_values")
# r2_per_output_topk  = r2_score(y_test, y_pred_topk, multioutput="raw_values")
# mse_mean_topk = mean_squared_error(y_test, y_pred_topk, multioutput="uniform_average")
# r2_mean_topk  = r2_score(y_test, y_pred_topk, multioutput="uniform_average")

# print("\n=== RF with Top-k Features Performance ===")
# for i, col in enumerate(target_cols):
#     print(f"{col}: MSE={mse_per_output_topk[i]:.4f}, R2={r2_per_output_topk[i]:.4f}")
# print(f"Overall (mean across 6 outputs): MSE={mse_mean_topk:.4f}, R2={r2_mean_topk:.4f}")

# ============================================================
# 6. Compute Friedman H-statistic using scikit-explain
# ============================================================

# IMPORTANT:
# Computing H-statistic on all 3e6 samples is infeasible.
# We subsample ~20k to ~50k samples for interpretability analysis.


# n_outputs = y_train.shape[1]

# all_hstats = {}   

# for out_i, rf in enumerate(models):
#     print(f"\n=== Computing H-statistics for output {out_i} ===")

#     # 1. select out limited samples
#     n_subsample = min(50_000, X_train.shape[0])
#     rng = np.random.default_rng(42)  # 
#     sub_idx = rng.choice(X_train.shape[0], size=n_subsample, replace=False)

#     X_sub = pd.DataFrame(X_train[sub_idx], columns=param_names)
#     y_sub = y_train[sub_idx, out_i]   #

#     # 2. Initialize ExplainToolkit
#     explainer = ExplainToolkit(
#         estimators=[(f"rf_y{out_i}", rf)],  # (name, estimator)
#         X=X_sub,
#         y=y_sub,
#         estimator_output="raw",
#     )

#     print("  Computing 1D and 2D partial dependence...")
#     pd_1d = explainer.pd(features="all",   n_bins=25, subsample=1.0)
#     pd_2d = explainer.pd(features="all_2d", n_bins=20, subsample=1.0)

#     # 3. Friedman H-statistic
#     pairs = list(pd_2d.coords["features"].values)
#     hstat_ds = explainer.friedman_h_stat(pd_1d, pd_2d, features=pairs)

#     print("  Done. Example:")
#     print(hstat_ds.isel(features=0).values)

#     # 4. convert to DataFrame
#     h_df = hstat_ds.to_dataframe().reset_index()
#     # 
#     h_df["output_index"] = out_i

#     all_hstats[out_i] = h_df

# # check
# print("\n=== H-statistics for output 0 (head) ===")
# print(all_hstats[0].head())

# # ============================================================
# # 7. Package and save all RF results into a single file
# # ============================================================

# # 7-1. Collect feature importances for each output
# importances_list = [rf.feature_importances_ for rf in models]    # list of arrays, shape (n_features,)
# importances_array = np.vstack(importances_list)                  # shape (n_outputs, n_features)
# mean_importance = np.mean(importances_array, axis=0)             # averaged importance across outputs



# # 7-2. Combine all H-statistic DataFrames into a single table
# import pandas as pd
# h_all_df = pd.concat(all_hstats.values(), ignore_index=True)     # shape = (pairs * n_outputs, ...)

# # 7-3. Pack everything into one Python dictionary
# results = {
#     "models": models,   # list of RandomForestRegressor models (one per output)

#     "meta": {
#         "mat_path": mat_path,
#         "rf_feature_indices": rf_feature_indices,
#         "param_names": param_names,
#         "resid_col_ids": resid_col_ids,
#         "resid_col_names": resid_col_names,
#     },

#     "metrics": {
#         "mse_per_output": mse_per_output,   # array of shape (n_outputs,)
#         "r2_per_output":  r2_per_output,    # array of shape (n_outputs,)
#         "mse_mean": mse_mean,
#         "r2_mean":  r2_mean,
#     },

#     "feature_importance": {
#         "per_output": importances_array,    # shape (n_outputs, n_features)
#         "mean":       mean_importance,      # shape (n_features,)
#     },

#     "hstats": {
#         "per_output_df": all_hstats,        # dict: {output_index: DataFrame}
#         "all_df":        h_all_df,          # merged H-statistic DataFrame
#     },
# }

# # 7-4. Save everything into a single .joblib file
# joblib.dump(results, "rf_all_results.joblib")
# print("Saved packed RF results to rf_all_results.joblib")


# print("end")



