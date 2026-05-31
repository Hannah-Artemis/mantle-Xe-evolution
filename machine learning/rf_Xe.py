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

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
print("Working directory:", os.getcwd())

# from skexplain.main import ExplainToolkit

# 
print("Hello, Python is working!")

# random forest training and feature importance analysis for Xe data (stage 3)

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

X_rf = X_rf.astype("float32")
Y_resid = Y_resid.astype("float32")

# ============================================================
# 2. Train–test split
# ============================================================

X_train, X_test, y_train, y_test = train_test_split(
    X_rf, Y_resid, test_size=0.2, random_state=42
)


# ============================================================
# 3. Train baseline Random Forest Regressor
# ============================================================

# RandomForestRegressor supports multi-output directly.
# rf = RandomForestRegressor(
#     n_estimators=200,     # number of trees
#     max_depth=None,
#     n_jobs=-1,            # use all CPU cores
#     max_samples=500_000, # max smaples used by a tree
#     random_state=42,
# )
# rf.fit(X_train, y_train)

models = []
for i in range(y_train.shape[1]):  # 6 output
    print(f"training target {i}")
    rf = RandomForestRegressor(
        n_estimators=150,
        max_depth=18,
        min_samples_leaf=10,
        max_samples=200_000,
        n_jobs=-1,
        random_state=0
    )
    rf.fit(X_train, y_train[:, i])
    models.append(rf)

# save
joblib.dump(models, "rf_models_all_outputs.joblib")
# # models = joblib.load("rf_models_all_outputs.joblib")


# ============================================================
# 4. Evaluate baseline model (MSE / R²)
# ============================================================

# y_pred = rf.predict(X_test)

y_pred_list = []

for rf in models:
    y_pred_list.append(rf.predict(X_test))

# 
y_pred = np.column_stack(y_pred_list)

# Raw values = one MSE per output dimension
mse_per_output = mean_squared_error(y_test, y_pred, multioutput="raw_values")
r2_per_output  = r2_score(y_test, y_pred, multioutput="raw_values")

# Uniform average gives a single value summarizing performance
mse_mean = mean_squared_error(y_test, y_pred, multioutput="uniform_average")
r2_mean  = r2_score(y_test, y_pred, multioutput="uniform_average")

print("\n=== Baseline RF Performance ===")
for i, col in enumerate(resid_col_names):
    print(f"{col}: MSE={mse_per_output[i]:.4f}, R2={r2_per_output[i]:.4f}")
print(f"Overall (mean across 6 outputs): MSE={mse_mean:.4f}, R2={r2_mean:.4f}")

# # save
# metrics_df = pd.DataFrame({
#     "output_name": resid_col_names,   # e.g., ["resid1", ..., "resid6"]
#     "mse": mse_per_output,
#     "r2": r2_per_output
# })
# metrics_df.to_csv("rf_metrics_per_output.csv", index=False)
# print("Saved: rf_metrics_per_output.csv")

# overall_df = pd.DataFrame({
#     "metric": ["mse_mean", "r2_mean"],
#     "value": [mse_mean, r2_mean]
# })
# overall_df.to_csv("rf_metrics_overall.csv", index=False)
# print("Saved: rf_metrics_overall.csv")



# ============================================================
# 5. Feature importance and top-k feature retraining
# ============================================================

# importances = rf.feature_importances_             # shape (12,)
# indices = np.argsort(importances)[::-1]           # sorted descending

# print("\n=== Feature Importances (Gini Importance) ===")
# for rank, idx in enumerate(indices):
#     print(f"Rank {rank+1}: {param_names[idx]} importance={importances[idx]:.4f}")



for out_i, rf in enumerate(models):
    importances = rf.feature_importances_            # shape (12,)
    indices = np.argsort(importances)[::-1]          # ranking

    print(f"\n=== Feature Importances (Impurity) for output {out_i} ===")
    for rank, idx in enumerate(indices):
        print(f"Rank {rank+1}: {param_names[idx]} importance={importances[idx]:.4f}")

    # select out limited samples for permutation importance
    n_subsample = min(50_000, X_test.shape[0])
    rng = np.random.default_rng(42)  # 
    sub_idx = rng.choice(X_test.shape[0], size=n_subsample, replace=False)

    X_pi = pd.DataFrame(X_test[sub_idx], columns=param_names)
    y_pi = y_test[sub_idx, out_i]   #

    perm_importance = permutation_importance(
    rf,         # e.g. model for output 0
    X_pi,
    y_pi,
    n_repeats=10,
    random_state=42
    )
    indices = np.argsort(perm_importance.importances_mean)[::-1]          # ranking 

    print(f"\n=== Feature Importances (Permutation) for output {out_i} ===")
    for rank, idx in enumerate(indices):
        print(f"Rank {rank+1}: {param_names[idx]} importance={perm_importance.importances_mean[idx]:.4f}")


# ============================================================
# NNN. redo with classify tree
# ============================================================

threshold = 1.0
Y_bin = (np.abs(y_pred[:, i]) < threshold).astype(int)
from sklearn.ensemble import RandomForestClassifier

clf = RandomForestClassifier(
    n_estimators=200,
    max_depth=None,
    n_jobs=-1,
    random_state=42
)
clf.fit(X_rf, Y_bin)

perm = permutation_importance(
    clf,
    X_rf,
    Y_bin,
    n_repeats=10,
    random_state=42
)

indices = np.argsort(perm.importances_mean)[::-1]
for idx in indices:
    print(param_names[idx], perm.importances_mean[idx])




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


n_outputs = y_train.shape[1]

all_hstats = {}   

for out_i, rf in enumerate(models):
    print(f"\n=== Computing H-statistics for output {out_i} ===")

    # 1. select out limited samples
    n_subsample = min(50_000, X_train.shape[0])
    rng = np.random.default_rng(42)  # 
    sub_idx = rng.choice(X_train.shape[0], size=n_subsample, replace=False)

    X_sub = pd.DataFrame(X_train[sub_idx], columns=param_names)
    y_sub = y_train[sub_idx, out_i]   #

    # 2. Initialize ExplainToolkit
    explainer = ExplainToolkit(
        estimators=[(f"rf_y{out_i}", rf)],  # (name, estimator)
        X=X_sub,
        y=y_sub,
        estimator_output="raw",
    )

    print("  Computing 1D and 2D partial dependence...")
    pd_1d = explainer.pd(features="all",   n_bins=25, subsample=1.0)
    pd_2d = explainer.pd(features="all_2d", n_bins=20, subsample=1.0)

    # 3. Friedman H-statistic
    pairs = list(pd_2d.coords["features"].values)
    hstat_ds = explainer.friedman_h_stat(pd_1d, pd_2d, features=pairs)

    print("  Done. Example:")
    print(hstat_ds.isel(features=0).values)

    # 4. convert to DataFrame
    h_df = hstat_ds.to_dataframe().reset_index()
    # 
    h_df["output_index"] = out_i

    all_hstats[out_i] = h_df

# check
print("\n=== H-statistics for output 0 (head) ===")
print(all_hstats[0].head())

# ============================================================
# 7. Package and save all RF results into a single file
# ============================================================

# 7-1. Collect feature importances for each output
importances_list = [rf.feature_importances_ for rf in models]    # list of arrays, shape (n_features,)
importances_array = np.vstack(importances_list)                  # shape (n_outputs, n_features)
mean_importance = np.mean(importances_array, axis=0)             # averaged importance across outputs



# 7-2. Combine all H-statistic DataFrames into a single table
import pandas as pd
h_all_df = pd.concat(all_hstats.values(), ignore_index=True)     # shape = (pairs * n_outputs, ...)

# 7-3. Pack everything into one Python dictionary
results = {
    "models": models,   # list of RandomForestRegressor models (one per output)

    "meta": {
        "mat_path": mat_path,
        "rf_feature_indices": rf_feature_indices,
        "param_names": param_names,
        "resid_col_ids": resid_col_ids,
        "resid_col_names": resid_col_names,
    },

    "metrics": {
        "mse_per_output": mse_per_output,   # array of shape (n_outputs,)
        "r2_per_output":  r2_per_output,    # array of shape (n_outputs,)
        "mse_mean": mse_mean,
        "r2_mean":  r2_mean,
    },

    "feature_importance": {
        "per_output": importances_array,    # shape (n_outputs, n_features)
        "mean":       mean_importance,      # shape (n_features,)
    },

    "hstats": {
        "per_output_df": all_hstats,        # dict: {output_index: DataFrame}
        "all_df":        h_all_df,          # merged H-statistic DataFrame
    },
}

# 7-4. Save everything into a single .joblib file
joblib.dump(results, "rf_all_results.joblib")
print("Saved packed RF results to rf_all_results.joblib")


print("end")



