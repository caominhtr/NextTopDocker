import pandas as pd
import numpy as np
import pickle
import sys
import os

if len(sys.argv) < 4:
    print("Usage: python logreg.py <model_nums...> <merge_file> <output_file>")
    sys.exit(1)

num_models = [int(x) for x in sys.argv[1:-2]]
merge_file = sys.argv[-2]
output_file = sys.argv[-1]

base_dir = os.path.dirname(__file__)


if not all(10 <= x <= 100 and x % 10 == 0 for x in num_models):
    raise ValueError("Please provide a valid model")
else:
    with open(os.path.join(base_dir, "model", "threshold.pkl"), "rb") as f:
        threshold = pickle.load(f)

    df = pd.read_csv(merge_file, header=None)
    df.columns = ["ID", "Pose", "Smina", "Gnina"]
    df = df.dropna()

    for num_model in num_models:

        with open(os.path.join(base_dir, "model", f"LogReg{num_model}.pkl"), "rb") as f:
            model = pickle.load(f)
        with open(os.path.join(base_dir, "scale", f"scaler_{num_model}.pkl"), "rb") as f:
            scale = pickle.load(f)

        thres = threshold[int(num_model / 10 - 1)]

        X_test = df[["Smina", "Gnina"]]
        X_test_scaled = scale.transform(X_test)
        y_pred_test = model.predict_proba(X_test_scaled)[:, 1]

        df[f"Pred_{num_model}"] = y_pred_test
        df[f"Near-native_{num_model}"] = df[f"Pred_{num_model}"].apply(lambda x: "Yes" if x >= thres else "No")

    df.to_csv(output_file, header=True, index=False)
