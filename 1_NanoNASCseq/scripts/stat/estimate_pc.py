#!/usr/bin/env python
import pickle
import optparse
import numpy as np
import pandas as pd
from scipy.stats import binom


def createMkn(Akn, p_e): # Left out from Akn
    M = np.zeros(Akn.shape)
    for n in range(Akn.shape[1]):
        for k in range(Akn.shape[0]):
            Ekn = np.sum(Akn[(k+1):, n]) * binom.pmf(k, n, p_e)
            if Ekn > 0.01 * Akn[k, n]:
                M[k, n] = 1
    return M


def EstepAkn(Akn, Mkn, p_c): # Alters Akn in place - modify initial step
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            if Mkn[k, n] == 1:
                num = 0
                denom = 0
                for kp in range(Mkn.shape[0]):
                    if Mkn[kp, n] == 1:
                        num = num + binom.pmf(k, n, p_c) * Akn[kp, n]
                        denom = denom + binom.pmf(kp, n, p_c)
                Akn[k, n] = num / denom
    return Akn


def MstepP_c(Akn):
    num = 0
    denom = 0
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            num = num + k * Akn[k, n]
            denom = denom + n * Akn[k, n]
    p_c = num / denom
    return p_c


def load_mismatch_ratio(path):
    return pd.read_csv(path, sep="\t")


def load_t_tc_count(path):
    array = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t, tc, count = map(int, line.strip("\n").split("\t"))
            array.append([t, tc, count])
    max_t = 0
    max_tc = 0
    if len(array) > 0:
        max_t = max([x[0] for x in array])
        max_tc = max([x[1] for x in array])
    # dtype can not be int, because some time the value wil be NaN in the downstream code.
    Akn = np.zeros((max_tc + 1, max_t + 1))
    for t, tc, count in array:
        Akn[tc, t] = count
    return Akn


def load_pe_model(path):
    return pickle.load(open(path, "rb"))


def predict_pc(p_e, Akn):
    l = 0
    r = 1
    p_c0 = (l + r) / 2
    p_c = p_c0
    Mkn = createMkn(Akn, p_e)
    while r - l >= 10e-8:
        Akn = EstepAkn(Akn, Mkn, p_c)
        p_c_old = p_c
        p_c = MstepP_c(Akn)
        if p_c < p_c_old:
            r = p_c
        else:
            l = p_c
    return p_c


def main():
    parser = optparse.OptionParser(usage="%prog -r ratios.tsv -e matrix.tsv -m model.pkl -o out.tsv")
    parser.add_option("-r", "--ratio", dest="ratio", metavar="PATH", help="Mismatch ratio of all mismatch types at 1 row.")
    parser.add_option("-e", "--event", dest="event", metavar="PATH", help="Matrix of the number of T-C mismatch event under certain number of T.")
    parser.add_option("-m", "--model", dest="model", metavar="PATH", help="The model for predicting the Pe value from misamtch ratios.")
    parser.add_option("-o", "--outfile", dest="outfile", metavar="PATH", help="Output file.")
    options, args = parser.parse_args()
    outfile = options.outfile
    if outfile is None:
        print("Must provide outfile!")
        exit(1)
    ratios = load_mismatch_ratio(options.ratio)
    Akn = load_t_tc_count(options.event)
    model = load_pe_model(options.model)
    features = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TG']
    features = ["%s_ratio" % x for x in features]
    try:
        p_e = model.predict(ratios[features])[0]
        p_c = predict_pc(p_e, Akn)
    except ValueError as e:
        print(e)
        p_e = np.NaN
        p_c = np.NaN
    ratios["Pe"] = p_e
    ratios["Pc"] = p_c
    ratios["SNR"] = ratios["Pc"] / ratios["Pe"]
    ratios.to_csv(outfile, sep="\t", index=False)
    

if __name__ == "__main__":
    main()