# ml_classification_pipeline.py
# Purpose:
#   Reproducible binary classification workflow using a public
#   person-level analytic table.
#
# Inputs:
#   output/nhanes_ckd_analytic.csv
#
# Outputs:
#   output/ml_model_metrics.csv

from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


DATA_PATH = Path("output/nhanes_ckd_analytic.csv")
OUT_PATH = Path("output/ml_model_metrics.csv")


def make_stage(row: pd.Series) -> float:
    """Construct a simple ordinal CKD severity score from eGFR and urine ACR."""
    egfr = row.get("egfr", np.nan)
    acr = row.get("urine_acr", np.nan)

    if pd.isna(egfr) or pd.isna(acr):
        return np.nan
    if egfr >= 60 and acr < 30:
        return 0
    if egfr >= 90 and acr >= 30:
        return 1
    if 60 <= egfr < 90 and acr >= 30:
        return 2
    if 45 <= egfr < 60:
        return 3
    if 30 <= egfr < 45:
        return 4
    if 15 <= egfr < 30:
        return 5
    return 6


def evaluate_classifier(name: str, model: Pipeline, X_test: pd.DataFrame, y_test: pd.Series) -> dict:
    """Compute standard binary-classification metrics."""
    y_pred = model.predict(X_test)
    y_score = model.predict_proba(X_test)[:, 1]

    tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0, 1]).ravel()

    return {
        "model": name,
        "n_test": int(len(y_test)),
        "tp": int(tp),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "accuracy": accuracy_score(y_test, y_pred),
        "sensitivity": tp / (tp + fn) if (tp + fn) else np.nan,
        "specificity": tn / (tn + fp) if (tn + fp) else np.nan,
        "precision": precision_score(y_test, y_pred, zero_division=0),
        "f1": f1_score(y_test, y_pred, zero_division=0),
        "auc": roc_auc_score(y_test, y_score),
    }


def main() -> None:
    df = pd.read_csv(DATA_PATH)

    df["ckd_stage"] = df.apply(make_stage, axis=1)
    df["ckd_any"] = np.where(df["ckd_stage"] > 0, 1, 0)
    df.loc[df["ckd_stage"].isna(), "ckd_any"] = np.nan

    features = [
        "age",
        "BMXBMI",
        "diabetes",
        "hypertension",
    ]
    features = [c for c in features if c in df.columns]

    model_df = df[features + ["ckd_any"]].dropna(subset=["ckd_any"]).copy()
    X = model_df[features]
    y = model_df["ckd_any"].astype(int)

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.30,
        random_state=123,
        stratify=y,
    )

    numeric_preprocess = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )

    preprocess = ColumnTransformer(
        transformers=[
            ("numeric", numeric_preprocess, features),
        ],
        remainder="drop",
    )

    logistic = Pipeline(
        steps=[
            ("preprocess", preprocess),
            ("classifier", LogisticRegression(max_iter=2000)),
        ]
    )

    random_forest = Pipeline(
        steps=[
            ("preprocess", preprocess),
            (
                "classifier",
                RandomForestClassifier(
                    n_estimators=500,
                    min_samples_leaf=10,
                    class_weight="balanced",
                    random_state=123,
                    n_jobs=-1,
                ),
            ),
        ]
    )

    models = {
        "logistic_regression": logistic,
        "random_forest": random_forest,
    }

    results = []
    for name, model in models.items():
        model.fit(X_train, y_train)
        results.append(evaluate_classifier(name, model, X_test, y_test))

    metrics = pd.DataFrame(results)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(OUT_PATH, index=False)

    print(metrics.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
