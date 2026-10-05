import warnings
warnings.filterwarnings("ignore")

from sklearn.metrics import roc_auc_score, balanced_accuracy_score, precision_score, confusion_matrix, f1_score
from joblib import Parallel, delayed
from sklearn.utils import shuffle
from xgboost import XGBClassifier
from sklearn import svm
import lightgbm as lgb
import anndata as ad
import scanpy as sc
import optuna


cancer_type_map = {
"gse_261194": 0,
"gse_268201": 0,
"gse_109761": 0,
"gse_75367": 0,
"gse_51827": 0,
"gse_55807": 0,
"gse_67939": 0,
"gse_86978": 0,
"gse_111065": 0,
"gse_111842": 0,
"gse_67980": 1,
"gse_74639": 2,
"gse_60407": 3,
"gse_117623": 4,
"gse_176078": 0,
"our_spikes": 0,
"gse_125449": 4,
"gse_131907": 2,
"gse_181294": 1,
"gse_194247": 3,
}


def add_dataset_and_cancer_type(
    adata,
    dataset_column="dataset",
    cancer_column="label"
):
    adata.obs[dataset_column] = [
        "_".join(str(cell_name).split("_")[:2])
        for cell_name in adata.obs_names
    ]

    adata.obs[cancer_column] = (
        adata.obs[dataset_column]
        .map(cancer_type_map)
        .fillna("Unknown")
    )

    return adata


def load_set(set_path):
    adata = sc.read_h5ad(set_path)
    adata = add_dataset_and_cancer_type(adata)

    X = adata.X
    y = adata.obs['label'].values
    X, y = shuffle(X, y, random_state=42)
    return X, y

test_data, test_classes = load_set(set_path="./ready_data/test_tumor_multiclass.h5ad")

train_data, train_classes = load_set(set_path="./ready_data/train_tumor_multiclass.h5ad")

val_data, val_classes = load_set(set_path="./ready_data/val_tumor_multiclass.h5ad")

test_data_ctc, test_classes_ctc = load_set(set_path="./ready_data/ctc_normalized_data.h5ad")


def objective(trial,
              model_name,
              train_data,
              train_classes,
              val_data,
              val_classes,
              test_data,
              test_classes,
              test_data_ctc,
              test_classes_ctc):

    if model_name == 'lightgbm':
        num_leaves = trial.suggest_int('num_leaves', 2, 128)
        learning_rate = trial.suggest_float('learning_rate', 1e-5, 1.0, log=True)
        n_estimators = trial.suggest_int('n_estimators', 50, 1000)
        max_depth = trial.suggest_int('max_depth', 2, 64)
        reg_alpha = trial.suggest_float('reg_alpha', 1e-8, 1.0, log=True)
        reg_lambda = trial.suggest_float('reg_lambda', 1e-8, 1.0, log=True)
        min_child_samples = trial.suggest_int('min_child_samples', 2, 64)
        is_unbalance = trial.suggest_categorical('is_unbalance', [True, False])
        model = lgb.LGBMClassifier(num_leaves=num_leaves, max_depth=max_depth, is_unbalance=is_unbalance,
                                learning_rate=learning_rate, n_estimators=n_estimators, random_state=42, verbosity=-1,
                                reg_alpha=reg_alpha, reg_lambda=reg_lambda, min_child_samples=min_child_samples)
    
    elif model_name == 'xgboost':
        max_depth = trial.suggest_int('max_depth', 2, 32)
        learning_rate = trial.suggest_float('learning_rate', 1e-5, 1.0, log=True)
        n_estimators = trial.suggest_int('n_estimators', 50, 500)
        reg_alpha = trial.suggest_float('reg_alpha', 1e-8, 10.0, log=True)
        reg_lambda = trial.suggest_float('reg_lambda', 1e-8, 10.0, log=True)
        min_child_weight = trial.suggest_float('min_child_weight', 1e-3, 10.0, log=True)
        model = XGBClassifier(max_depth=max_depth, learning_rate=learning_rate, n_estimators=n_estimators, 
                            reg_alpha=reg_alpha, reg_lambda=reg_lambda, min_child_weight=min_child_weight,
                            seed=42)
    
    model.fit(train_data, train_classes)
    val_pred = model.predict(val_data)
    test_pred = model.predict(test_data)
    test_pred_ctc = model.predict(test_data_ctc)

    val_balanced_accuracy = balanced_accuracy_score(val_classes, val_pred)
    test_balanced_accuracy = balanced_accuracy_score(test_classes, test_pred)
    test_ctc_balanced_accuracy = balanced_accuracy_score(test_classes_ctc, test_pred_ctc)

    val_macro_f1 = f1_score(val_classes, val_pred, average="macro", zero_division=0)
    test_macro_f1 = f1_score(test_classes, test_pred, average="macro", zero_division=0)
    test_ctc_macro_f1 = f1_score(test_classes_ctc, test_pred_ctc, average="macro", zero_division=0)

    val_weighted_f1 = f1_score(val_classes, val_pred, average="weighted", zero_division=0)
    test_weighted_f1 = f1_score(test_classes, test_pred, average="weighted", zero_division=0)
    test_ctc_weighted_f1 = f1_score(test_classes_ctc, test_pred_ctc, average="weighted", zero_division=0)

    val_cm = confusion_matrix(val_classes, val_pred).tolist()
    test_cm = confusion_matrix(test_classes, test_pred).tolist()
    test_ctc_cm = confusion_matrix(test_classes_ctc, test_pred_ctc).tolist()

    trial.set_user_attr('test_balanced_accuracy', test_balanced_accuracy)
    trial.set_user_attr('test_macro_f1', test_macro_f1)
    trial.set_user_attr('test_weighted_f1', test_weighted_f1)
    trial.set_user_attr('test_confusion_matrix', test_cm)

    trial.set_user_attr('val_balanced_accuracy', val_balanced_accuracy)
    trial.set_user_attr('val_macro_f1', val_macro_f1)
    trial.set_user_attr('val_weighted_f1', val_weighted_f1)
    trial.set_user_attr('val_confusion_matrix', val_cm)

    trial.set_user_attr('test_ctc_balanced_accuracy', test_ctc_balanced_accuracy)
    trial.set_user_attr('test_ctc_macro_f1', test_ctc_macro_f1)
    trial.set_user_attr('test_ctc_weighted_f1', test_ctc_weighted_f1)
    trial.set_user_attr('test_ctc_confusion_matrix', test_ctc_cm)

    return val_macro_f1


def run_study(model_name, n_trials=100):
    study = optuna.create_study(direction="maximize", study_name=f"study_{model_name}")

    study.optimize(lambda trial: objective(trial, model_name,
                                           train_data=train_data,
                                           train_classes=train_classes,
                                           val_data=val_data,
                                           val_classes=val_classes,
                                           test_data=test_data,
                                           test_classes=test_classes,
                                           test_data_ctc=test_data_ctc,
                                           test_classes_ctc=test_classes_ctc),
                                                                        n_trials=n_trials)
    print(f"Best parameters for {model_name}:", study.best_params)
    print(f"Best score for {model_name}:", study.best_value)
    
    df = study.trials_dataframe()
    df.to_csv(f"./optuna_results/model_{model_name}_multiclass_study_results.csv")
    
    return study 


model_names = ["xgboost", "lightgbm"]

studies = Parallel(n_jobs=len(model_names))(
    delayed(run_study)(model_name) for model_name in model_names
)
best_scores = {study.study_name: study.best_value for study in studies}
best_model = max(best_scores, key=best_scores.get)
print("Best overall model:", best_model)
print("Best score:", best_scores[best_model])