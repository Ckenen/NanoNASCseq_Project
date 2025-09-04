#!/usr/bin/env python
import sys
import os
import pickle
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn import metrics
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr


def plot(m, out_pdf):
    xs = m["TC_ratio"]
    ys = m["Pe"]
    plt.figure(figsize=(3, 3))
    plt.scatter(xs, ys, marker="o", edgecolor="black", alpha=1)
    lim = max(max(xs), max(ys)) * 1.2
    plt.plot([0, lim], [0, lim], lw=1, ls="--", color="grey")
    plt.xlim(0, lim)
    plt.ylim(0, lim)
    plt.xticks([])
    plt.yticks([])
    plt.xlabel("Observed T-to-C mismatch ratio")
    plt.ylabel("Predicted Pe")
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()
    
    
def train_pe_model(infile, outfile, p=False):
    
    

def main():
    """
    
    
    """
    infile, outfile = sys.argv[1:]
    assert outfile.endswith(".pkl")
    prefix = os.path.splitext(outfile)[0]
    out_pdf = prefix + ".pdf"
    out_pkl = prefix + ".pkl"
    out_tsv = prefix + ".tsv"

    features = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TG']
    features = ["%s_ratio" % x for x in features]

    m = pd.read_csv(infile, sep="\t")
    x = m[features]
    y = m["TC_ratio"]
    x_train, x_test, y_train, y_test = x, x, y, y
    print("Shape of train X:", x_train.shape)
    print("Shape of train Y:", y_train.shape)
    print("Shape of test X:", x_test.shape)
    print("Shape of test Y:", y_test.shape)
    
    lr = LinearRegression()
    lr.fit(x_train, y_train)
    y_pred = lr.predict(x_test)
    print(pearsonr(y_test, y_pred))
    print(spearmanr(y_test, y_pred))
    MSE = metrics.mean_squared_error(y_test, y_pred)
    RMSE = np.sqrt(metrics.mean_squared_error(y_test, y_pred))
    print('MSE:', MSE)
    print('RMSE:', RMSE)
    m["Pe"] = y_pred
    
    pickle.dump(lr, open(out_pkl, "wb"))
    m.to_csv(out_tsv)
    plot(m, out_pdf)
    
    
if __name__ == "__main__":
    main()