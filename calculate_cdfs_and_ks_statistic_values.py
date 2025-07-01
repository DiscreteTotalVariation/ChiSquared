# -*- coding: utf-8 -*-

import numpy as np
import cv2
from math import comb
from scipy.stats import chi2

import matplotlib
import matplotlib.pyplot as plt

from os import listdir, makedirs
from os.path import isdir, isfile, join
from tqdm import tqdm

ROOT_DIR="chi_data"

EXACT_DISTRIBUTIONS_DIR=join(ROOT_DIR, "chi_exact_distributions")
CDF_EXACT_VS_APPROXIMATED_DIR=join(ROOT_DIR, "chi_cdf_exact_vs_approximated")
KOLMOGOROV_SMIRNOV_STATISTIC_DIR=join(ROOT_DIR, "chi_kss")
HISTOGRAM_CONFIGURATION_PATTERN="N_%d_n_%d.txt"

def get_true_chi_statistic_value(integer_value, N, n):
    return integer_value*n/N-N

def get_and_cache_kolmogorov_smirnov_statistic(N, n, skip_possible_values=False):
    exact_distributions_dir=EXACT_DISTRIBUTIONS_DIR
    exact_vs_approximated_dir=CDF_EXACT_VS_APPROXIMATED_DIR
    kss_dir=KOLMOGOROV_SMIRNOV_STATISTIC_DIR

    makedirs(exact_vs_approximated_dir, exist_ok=True)
    makedirs(kss_dir, exist_ok=True)

    # check whther the comparison exists
    exact_vs_approximated_path=join(exact_vs_approximated_dir, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
    if isfile(exact_vs_approximated_path) is False:
        
        # load the distribution
        distribution=dict()
        input_path=join(exact_distributions_dir, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
        with open(input_path, "r", encoding="utf8") as f:
            lines=[l.strip("\r\n") for l in f.readlines()]
        for l in lines:
            a, b=list(map(int, l.split()))
            distribution[a]=b
        
        exact_cdf_values=[]
        approximated_cdf_values=[]
        chi_statistic_values=[]

        # calculate the exact CDF and the approximated CDF
        df=n-1
        numerator=0
        denominator=n**N
        for value, occurrence in distribution.items():
            if occurrence==0:
                continue
            chi_statistic=get_true_chi_statistic_value(integer_value=value, N=N, n=n)
            numerator+=occurrence
            exact_cdf_value=numerator/denominator
            approximated_cdf_value=chi2.cdf(chi_statistic, df)

            chi_statistic_values.append(chi_statistic)
            exact_cdf_values.append(exact_cdf_value)
            approximated_cdf_values.append(approximated_cdf_value)

        # save the CDF values
        with open(exact_vs_approximated_path, "w", encoding="utf8") as fo:
            for chi_statistic, exact_cdf_value, approximated_cdf_value in zip(chi_statistic_values, exact_cdf_values, approximated_cdf_values):
                fo.write(str(chi_statistic)+" "+str(exact_cdf_value)+" "+str(approximated_cdf_value)+"\n")

    # check the Kolmogorov-Smirnov statistic
    kss_path=join(kss_dir, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
    if isfile(kss_path) is False:
        d=np.loadtxt(exact_vs_approximated_path)
        if np.prod(d.shape)==3 and len(d.shape)==1:
            d=d[np.newaxis, :]
        if len(d.shape)==1:
            raise("A problem with the shape of the data.")

        kss=np.max(np.abs(d[:, 1]-d[:, 2]))
        with open(kss_path, "w", encoding="utf8") as fo:
            fo.write(str(kss))

    # finally, load the Kolmogorov-Smirnov statistic value
    with open(kss_path, "r", encoding="utf8") as f:
        ks=float(f.read())

    return ks

def main():

    input_dir=EXACT_DISTRIBUTIONS_DIR

    if isdir(input_dir) is False:
        return

    names=sorted(listdir(input_dir))

    show_status=True

    taker=names
    if show_status is True:
        taker=tqdm(names)
    for name in taker:
        parts=name[:name.rfind(".")].split("_")
        N=int(parts[1])
        n=int(parts[3])

        cdf_output_path=join(CDF_EXACT_VS_APPROXIMATED_DIR, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
        kss_output_path=join(KOLMOGOROV_SMIRNOV_STATISTIC_DIR, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))

        if isfile(cdf_output_path) is True and isfile(kss_output_path) is True:
            continue

        get_and_cache_kolmogorov_smirnov_statistic(N=N, n=n, skip_possible_values=True)

    pass
    
if __name__=="__main__":
    main()
