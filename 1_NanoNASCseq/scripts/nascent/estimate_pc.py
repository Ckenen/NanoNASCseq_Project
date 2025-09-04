#!/usr/bin/env python
import sys
import pickle
import numpy as np
from scipy.stats import binom


def createMkn(Akn, p_e): #Left out from Akn
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


def main():
    f_pkl, f_tsv = sys.argv[1:]

    with open(f_pkl, "rb") as f:
        umis = pickle.load(f)
    umis = list(filter(lambda umi: umi[1] >= 2, umis))

    p_e = 0.00015401167351422416

    if len(umis) > 0:
        max_length = max([umi[2]["T"] for umi in umis])
        
        Akn = np.zeros((max_length + 1, max_length + 1))
        for umi in umis:
            n = umi[2]["T"]
            k = umi[3]["T-C"]
            Akn[k, n] += 1
        Akn0 = Akn

        l = 0
        r = 1
        p_c0 = (l + r) / 2

        p_c = p_c0
        Mkn = createMkn(Akn0, p_e)
        Akn = Akn0

        while r - l >= 10e-8:
            # print(r-l)
            Akn = EstepAkn(Akn, Mkn, p_c)
            p_c_old = p_c
            p_c = MstepP_c(Akn)
            if p_c < p_c_old:
                r = p_c
            else:
                l = p_c
    else:
        p_c = np.nan

    # print(p_c, p_e, p_c / p_e, sep="\t")
    with open(f_tsv, "w+") as fw:
        fw.write("Pe\tPc\tPc/Pe\n")
        fw.write("%f\t%f\t%f\n" % (p_e, p_c, p_c/p_e))


if __name__ == "__main__":
    main()