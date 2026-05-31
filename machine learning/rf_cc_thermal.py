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

# random forest training and feature importance analysis for cc & thermal data (stage 1,2)

# ============================================================
# 1-1. Load dataset from .mat file
# ============================================================

# Load the .mat file
# mat = sio.loadmat("/Volumes/LH - drive/research/Xe-code/Xe/out_Xe_lhs_new_3e7_test_fixK_U_MC_lowerT_alot_Rac1100_new_lowf.mat")

mat_path = "/Volumes/LH - drive/research/Xe-code/cc/out_cc_growth_lhs_Krw_lowf.mat"

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
    'ts','\kappa_g','\kappa_r','Rs','Rp','frw',
    'rmse_formation_age','rmse_surface_age'
]

# Ensure consistency between col_names and results matrix size
assert len(col_names) == n_params, (
    f"col_names length ({len(col_names)}) does not match results columns ({n_params})"
)

# ======================================================================
# 1-3. Select 6 RF input parameters (X_rf) matching MATLAB indexing
# ======================================================================

rf_feature_indices = list(range(0, 6)) 
# X_rf = results[:, rf_feature_indices]
X_rf = results[rf_feature_indices,:]

print("X_rf shape:", X_rf.shape)
print("Selected RF feature names:", [col_names[i] for i in rf_feature_indices])


# ======================================================================
# 1-4. Define parameter names used for labeling 
# ======================================================================

param_names = [
   'ts','\kappa_g','\kappa_r','Rs','Rp','frw'
]

# ======================================================================
# 1-5. Identify indices of residual error columns (6 outputs)
#    Direct Python equivalent of MATLAB find(strcmp(...))
# ======================================================================

def find_col_index(name_list, target):
    """Return index of a column name (equivalent to MATLAB find(strcmp()))"""
    return name_list.index(target)

rmsef_idx       = find_col_index(col_names, 'rmse_formation_age')
rmses_idx = find_col_index(col_names, 'rmse_surface_age')


resid_col_ids = [
    rmsef_idx, 
    rmses_idx
]

resid_col_names = [
    'rmse_cc_formation_age', 'rmse_cc_surface_age'
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

# ============================================================
# 2. Train–test split
# ============================================================

X_train, X_test, y_train, y_test = train_test_split(
    X_rf, Y_resid, test_size=0.2, random_state=42
)



# # ============================================================
# # NNN. with classify tree
# # ============================================================

X_rf = X_scaled

threshold = np.array([0.1, 0.2])
Y_bin = (np.abs(Y_resid) < threshold).astype(int)
print("Y_bin shape:", Y_bin.shape)

# train/test split
X_train_c, X_test_c, y_train_c, y_test_c = train_test_split(
    X_rf,
    Y_bin,
    test_size=0.2,
    random_state=42
)



clf_models_list = []      # list of classifiers, one per output
perm_results_list = []    # store permutation_importance results (optional)

rng = np.random.default_rng(42)

for i in range(Y_resid.shape[1]):
    print(f"\n==============================")
    print(f"=== Output {i}: {resid_col_names[i]} ===")

    # Select binary labels for this output
    y_train_i = y_train_c[:, i]    # shape: (n_train,)
    y_test_i  = y_test_c[:, i]     # shape: (n_test,)

    # Train classifier for this single output
    clf = RandomForestClassifier(
        n_estimators=200,
        max_depth=None,
        n_jobs=-1,
        random_state=42,
    )
    clf.fit(X_train_c, y_train_i)
    clf_models_list.append(clf)

    # ----------------------------------------------------
    # 4. Classification metrics for this output
    # ----------------------------------------------------
    y_pred_cls = clf.predict(X_test_c)

    acc  = accuracy_score(y_test_i, y_pred_cls)
    # zero_division=0 avoids warnings if one class is missing in predictions
    prec = precision_score(y_test_i, y_pred_cls, zero_division=0)
    rec  = recall_score(y_test_i, y_pred_cls, zero_division=0)
    f1   = f1_score(y_test_i, y_pred_cls, zero_division=0)

    print("\n=== Classification metrics ===")
    print(f"Accuracy : {acc:.4f}")
    print(f"Precision: {prec:.4f}")
    print(f"Recall   : {rec:.4f}")
    print(f"F1-score : {f1:.4f}")

    cm = confusion_matrix(y_test_i, y_pred_cls)
    print("Confusion matrix:")
    print(cm)

    # 5. Permutation importance for this output
    #    (Use a subsample to avoid being too slow / OOM)
    # ----------------------------------------------------
    n_sub = min(10_000, X_test_c.shape[0])  # subsample size
    sub_idx = rng.choice(X_test_c.shape[0], size=n_sub, replace=False)
    X_pi = X_test_c[sub_idx]
    y_pi = y_test_i[sub_idx]

    perm = permutation_importance(
        clf,
        X_pi,
        y_pi,
        n_repeats=5,       # 
        random_state=42,
        n_jobs=1,          # safer on macOS; avoids kernel crashes
        scoring="f1",      # or "accuracy"
    )
    perm_results_list.append(perm)

    indices = np.argsort(perm.importances_mean)[::-1]

    print("\n=== Permutation feature importance (F1-based) ===")
    for rank, idx_f in enumerate(indices):
        print(
            f"Rank {rank+1}: {param_names[idx_f]} "
            f"importance={perm.importances_mean[idx_f]:.6f}"
        )


# ----- save
# ------------------------------------------------------------
# 1. Collect all metrics into a structured list
# ------------------------------------------------------------
metrics_list = []

for i in range(Y_resid.shape[1]):
    metrics_list.append({
        "output_index": i,
        "output_name": resid_col_names[i],
        "accuracy": float(accuracy_score(y_test_c[:, i],
                                        clf_models_list[i].predict(X_test_c))),
        "precision": float(precision_score(y_test_c[:, i],
                                          clf_models_list[i].predict(X_test_c),
                                          zero_division=0)),
        "recall": float(recall_score(y_test_c[:, i],
                                     clf_models_list[i].predict(X_test_c),
                                     zero_division=0)),
        "f1": float(f1_score(y_test_c[:, i],
                             clf_models_list[i].predict(X_test_c),
                             zero_division=0)),
        "confusion_matrix": confusion_matrix(y_test_c[:, i],
                                             clf_models_list[i].predict(X_test_c))
    })

# ------------------------------------------------------------
# 2. Extract permutation importance results into a dict
# ------------------------------------------------------------
perm_list = []
for i, perm in enumerate(perm_results_list):
    perm_list.append({
        "output_index": i,
        "output_name": resid_col_names[i],
        "importances_mean": perm.importances_mean,
        "importances_std": perm.importances_std,
        "sorted_indices": np.argsort(perm.importances_mean)[::-1]
    })

# ------------------------------------------------------------
# 3. Create final dictionary bundle
# ------------------------------------------------------------
save_bundle = {
    "mat_path": mat_path, 
    "models": clf_models_list,        # list of 6 classifiers
    "metrics": metrics_list,          # list of dicts
    "permutation_importance": perm_list,

    "param_names": param_names,       # feature names
    "output_names": resid_col_names,  # 6 outputs

    "threshold": threshold,           # classification threshold used
#    "scaler": scaler if "scaler" in globals() else None,  # store scaler if used
}

# ------------------------------------------------------------
# 4. Save with joblib
# ------------------------------------------------------------
joblib.dump(save_bundle, "rf_classification_results_cc.pkl")

print("\nSaved all results to rf_classification_results_cc.pkl")


# ------------------------------------------------------------
# reload the models
# -------------------------------------------------------------
# ------------------------------------------------------------
# 1. Load the saved bundle
# ------------------------------------------------------------
save_bundle = joblib.load("rf_classification_results_cc.pkl")
print("Loaded rf_classification_results_cc.pkl")

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
sio.savemat('rf_classification_results_cc.mat', mat_data, oned_as='column')
print("Saved to rf_classification_results_cc.mat")

# Print what was saved
print("\nSaved variables:")
for key in sorted(mat_data.keys()):
    if isinstance(mat_data[key], np.ndarray):
        print(f"  {key}: array shape {mat_data[key].shape}")
    else:
        print(f"  {key}: {type(mat_data[key]).__name__}")






# ============================================================
# 1-1. Load dataset from .mat file
# ============================================================

# Load the .mat file
# mat = sio.loadmat("/Volumes/LH - drive/research/Xe-code/Xe/out_Xe_lhs_new_3e7_test_fixK_U_MC_lowerT_alot_Rac1100_new_lowf.mat")

mat_path = "/Volumes/LH - drive/research/Xe-code/thermal/out_thermal_lhs_test_MC_lowerT_Rac1100_lowf.mat"

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
    'RMSE_pd','RMSE_t','RMSE_Q'
]

# Ensure consistency between col_names and results matrix size
assert len(col_names) == n_params, (
    f"col_names length ({len(col_names)}) does not match results columns ({n_params})"
)

# ======================================================================
# 1-3. Select 6 RF input parameters (X_rf) matching MATLAB indexing
# ======================================================================

rf_feature_indices = list(range(0, 7))+list(range(8, 10)) 
# X_rf = results[:, rf_feature_indices]
X_rf = results[rf_feature_indices,:]

print("X_rf shape:", X_rf.shape)
print("Selected RF feature names:", [col_names[i] for i in rf_feature_indices])


# ======================================================================
# 1-4. Define parameter names used for labeling 
# ======================================================================

param_names = [
  'ts_cc','kappa_gcc','kappa_rcc','Rs','Rp',
    'Ubse','Qc_pd', 'Ti', 'eta_ref'
]

# ======================================================================
# 1-5. Identify indices of residual error columns (6 outputs)
#    Direct Python equivalent of MATLAB find(strcmp(...))
# ======================================================================

def find_col_index(name_list, target):
    """Return index of a column name (equivalent to MATLAB find(strcmp()))"""
    return name_list.index(target)


rmseTpd_idx       = find_col_index(col_names, 'RMSE_pd')
rmseTt_idx = find_col_index(col_names, 'RMSE_t')
rmseQ_idx = find_col_index(col_names, 'RMSE_Q') 

resid_col_ids = [
    rmseTpd_idx, 
    rmseTt_idx,
    rmseQ_idx
]

resid_col_names = [
    'rmse_Tpd', 'rmse_Tt','rmse_Q'
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

# ============================================================
# 2. Train–test split
# ============================================================

X_train, X_test, y_train, y_test = train_test_split(
    X_rf, Y_resid, test_size=0.2, random_state=42
)



# # ============================================================
# # NNN. with classify tree
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



clf_models_list = []      # list of classifiers, one per output
perm_results_list = []    # store permutation_importance results (optional)

rng = np.random.default_rng(42)

for i in range(Y_resid.shape[1]):
    print(f"\n==============================")
    print(f"=== Output {i}: {resid_col_names[i]} ===")

    # Select binary labels for this output
    y_train_i = y_train_c[:, i]    # shape: (n_train,)
    y_test_i  = y_test_c[:, i]     # shape: (n_test,)

    # Train classifier for this single output
    clf = RandomForestClassifier(
        n_estimators=200,
        max_depth=None,
        n_jobs=-1,
        random_state=42,
    )
    clf.fit(X_train_c, y_train_i)
    clf_models_list.append(clf)

    # ----------------------------------------------------
    # 4. Classification metrics for this output
    # ----------------------------------------------------
    y_pred_cls = clf.predict(X_test_c)

    acc  = accuracy_score(y_test_i, y_pred_cls)
    # zero_division=0 avoids warnings if one class is missing in predictions
    prec = precision_score(y_test_i, y_pred_cls, zero_division=0)
    rec  = recall_score(y_test_i, y_pred_cls, zero_division=0)
    f1   = f1_score(y_test_i, y_pred_cls, zero_division=0)

    print("\n=== Classification metrics ===")
    print(f"Accuracy : {acc:.4f}")
    print(f"Precision: {prec:.4f}")
    print(f"Recall   : {rec:.4f}")
    print(f"F1-score : {f1:.4f}")

    cm = confusion_matrix(y_test_i, y_pred_cls)
    print("Confusion matrix:")
    print(cm)

    # 5. Permutation importance for this output
    #    (Use a subsample to avoid being too slow / OOM)
    # ----------------------------------------------------
    n_sub = min(10_000, X_test_c.shape[0])  # subsample size
    sub_idx = rng.choice(X_test_c.shape[0], size=n_sub, replace=False)
    X_pi = X_test_c[sub_idx]
    y_pi = y_test_i[sub_idx]

    perm = permutation_importance(
        clf,
        X_pi,
        y_pi,
        n_repeats=5,       # you can increase if it is fast enough
        random_state=42, 
        n_jobs=1,          # safer on macOS; avoids kernel crashes
        scoring="f1",      # or "accuracy"
    )
    perm_results_list.append(perm)

    indices = np.argsort(perm.importances_mean)[::-1]

    print("\n=== Permutation feature importance (F1-based) ===")
    for rank, idx_f in enumerate(indices):
        print(
            f"Rank {rank+1}: {param_names[idx_f]} "
            f"importance={perm.importances_mean[idx_f]:.6f}"
        )


# ----- save
# ------------------------------------------------------------
# 1. Collect all metrics into a structured list
# ------------------------------------------------------------
metrics_list = []

for i in range(Y_resid.shape[1]):
    metrics_list.append({
        "output_index": i,
        "output_name": resid_col_names[i],
        "accuracy": float(accuracy_score(y_test_c[:, i],
                                        clf_models_list[i].predict(X_test_c))),
        "precision": float(precision_score(y_test_c[:, i],
                                          clf_models_list[i].predict(X_test_c),
                                          zero_division=0)),
        "recall": float(recall_score(y_test_c[:, i],
                                     clf_models_list[i].predict(X_test_c),
                                     zero_division=0)),
        "f1": float(f1_score(y_test_c[:, i],
                             clf_models_list[i].predict(X_test_c),
                             zero_division=0)),
        "confusion_matrix": confusion_matrix(y_test_c[:, i],
                                             clf_models_list[i].predict(X_test_c))
    })

# ------------------------------------------------------------
# 2. Extract permutation importance results into a dict
# ------------------------------------------------------------
perm_list = []
for i, perm in enumerate(perm_results_list):
    perm_list.append({
        "output_index": i,
        "output_name": resid_col_names[i],
        "importances_mean": perm.importances_mean,
        "importances_std": perm.importances_std,
        "sorted_indices": np.argsort(perm.importances_mean)[::-1]
    })

# ------------------------------------------------------------
# 3. Create final dictionary bundle
# ------------------------------------------------------------
save_bundle = {
    "mat_path": mat_path, 
    "models": clf_models_list,        # list of 6 classifiers
    "metrics": metrics_list,          # list of dicts
    "permutation_importance": perm_list,

    "param_names": param_names,       # feature names
    "output_names": resid_col_names,  # 6 outputs

    "threshold": threshold,           # classification threshold used
#    "scaler": scaler if "scaler" in globals() else None,  # store scaler if used
}

# ------------------------------------------------------------
# 4. Save with joblib
# ------------------------------------------------------------
joblib.dump(save_bundle, "rf_classification_results_T.pkl")

print("\nSaved all results to rf_classification_results_T.pkl")


# ------------------------------------------------------------
# reload the models
# -------------------------------------------------------------
# ------------------------------------------------------------
# 1. Load the saved bundle
# ------------------------------------------------------------
save_bundle = joblib.load("rf_classification_results_T.pkl")
print("Loaded rf_classification_results_T.pkl")

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
sio.savemat('rf_classification_results_T.mat', mat_data, oned_as='column')
print("Saved to rf_classification_results_T.mat")

# Print what was saved
print("\nSaved variables:")
for key in sorted(mat_data.keys()):
    if isinstance(mat_data[key], np.ndarray):
        print(f"  {key}: array shape {mat_data[key].shape}")
    else:
        print(f"  {key}: {type(mat_data[key]).__name__}")


