# -*- coding: utf-8 -*-

import argparse

import numpy as np

from os import listdir, makedirs
from os.path import isfile, join
from tqdm import tqdm

ROOT_DIR="chi_data"

POSSIBLE_VALUES_DIR=join(ROOT_DIR, "chi_possible_values")
EXACT_DISTRIBUTIONS_DIR=join(ROOT_DIR, "chi_exact_distributions")

HISTOGRAM_CONFIGURATION_PATTERN="N_%d_n_%d.txt"

def precompute_combinations(n_max):
    combs=[[0]*(n_max+1) for _ in range(n_max+1)]
    for n in range(n_max+1):
        combs[n][0]=combs[n][n]=1
        for k in range(1, n):
            combs[n][k]=combs[n-1][k-1]+combs[n-1][k]
    return combs

def calculate_possible_values(N, n, only_final=False, use_cached=True, show_status=False):
    if use_cached is True:
        input_path=join(POSSIBLE_VALUES_DIR, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
        if isfile(input_path) is True:
            with open(input_path, "r", encoding="utf8") as f:
                lines=[l.strip("\r\n") for l in f.readlines()]
            return list(map(int, lines))
    
    # the states of the histograms of the currently observed size in terms of number fo bins will be stored as a combination of the subsample size and the sum of squared bin values
    # these two values, e.g., a and b will be stored as a single integer i=a*factor+b and b will always be less than N+1, since N is the maximum sample size
    factor=N+1
    
    previous=[0]
    # the maximum possible value of the sum of squared bin values is N*N
    previous_flags=np.zeros(factor*(N*N+1))
    
    all_values=set()
    taker=list(range(n))
    if show_status is True:
        print("Calculating possible values...")
        taker=tqdm(taker)
    
    # the following is a modified breadth-first search algorithm
    for i in taker:
        next=[]
        for previous_combination in previous:
            # separating the previous combination into the sample size and the sum of squared bin values
            previous_N=previous_combination%factor
            previous_value=previous_combination//factor
            # trying to add a new bin to the histogram with all possible values based on what remaing after subtracting the previously taken sample size
            for current_N in range(0, N-previous_N+1):
                next_N=previous_N+current_N
                next_value=previous_value+current_N**2
                all_values.add(next_value)
                # combining the new sample size and the new sum of squared bin values into a single integer to describe the next state
                next_combination=next_N+factor*next_value
                if previous_flags[next_combination]==False:
                    next.append(next_combination)
                    previous_flags[next_combination]=True
        previous=next
    
    # selecting the required possbile sums of squared bin values
    values=set()
    if only_final is True:
        for current_combination in next:
            current_N=current_combination%factor
            value=current_combination//factor
            if current_N==N:
                values.add(value)
    else:
        for value in all_values:
            values.add(value)
    
    return sorted(values)

def cache_exact_distribution_up_to_both_ns(upper_N, upper_n, show_status=False, use_cached=True, skip_zeros=True):
    if (upper_n==0):
        return

    # calculating the possible values of the sums of squared bin values to speed up the calculations
    possible_values=calculate_possible_values(N=upper_N, n=upper_n, only_final=False, use_cached=use_cached, show_status=show_status)
    value_idx={x:i for i, x in enumerate(possible_values)}
    values_count=len(possible_values)
    
    # precomputing the binomial coefficients for all possible sample sizes - taking combs[a][b] will have the same results as calling math.comb(a, b)
    combs=precompute_combinations(upper_N)

    values_dir=POSSIBLE_VALUES_DIR
    values_path=join(values_dir, HISTOGRAM_CONFIGURATION_PATTERN%(upper_N, upper_n))
    if isfile(values_path) is False:
        makedirs(values_dir, exist_ok=True)
        possible_numerators=sorted(possible_values)
        with open(values_path, "w", encoding="utf8") as fo:
            for numerator in possible_numerators:
                fo.write(str(numerator)+"\n")

    value_idx={x:i for i, x in enumerate(possible_values)}
    values_count=len(possible_values)
    
    # preparing the previous and current counts of ways to form histograms with a specific sample size and a sum of squared bin values
    previous_counts=[]
    current_counts=[]
    for i in range(upper_N+1):
        current_counts.append([0]*values_count)
        previous_counts.append([0]*values_count)
    
    for i in range(upper_N+1):
        current_counts[i][value_idx[i**2]]=combs[upper_N][i]
    
    exact_distributions_dir=EXACT_DISTRIBUTIONS_DIR
    makedirs(exact_distributions_dir, exist_ok=True)

    Ns=list(range(1, upper_N+1))
    ns=list(range(2, upper_n+1))
    missing_count=0
    for N in Ns:
        for n in ns:
            exact_distribution_path=join(exact_distributions_dir, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
            if isfile(exact_distribution_path) is False:
                missing_count+=1

    if missing_count==0:
        return
    
    taker=list(range(1, n))
    if show_status is True:
        print("Calculating the exact distribution(s)...")
        taker=tqdm(taker)

    for idx in taker:
        # to save memory, only two counts variables are used and here their role is swapped
        previous_counts, current_counts=current_counts, previous_counts
        # clearing the current counts
        for i in range(len(current_counts)):
            for j in range(len(current_counts[i])):
                current_counts[i][j]=0
        # going over all previously calculated counts for specific sample sizes (the sample size is in the variable previous_sum)
        for previous_sum in range(0, N+1):
            remaining=N-previous_sum
            # going over all sums of squared bin values that could have appeared for the previous samples
            for previous_value_idx in range(0, values_count):
                # the previous sum of squared bin values
                previous_value=possible_values[previous_value_idx]
                previous_count=previous_counts[previous_sum][previous_value_idx]
                # this is done only if there are some previous counts, which speeds things up
                if (previous_count>0):
                    for new_last in range(0, remaining+1):
                        # the new sample size is the previous sample size plus the currently taken part of the remaining sample
                        new_sum=previous_sum+new_last
                        # increasing the sum of the squared bins by the square of the newly added bin
                        new_value=previous_value+new_last**2
                        # adding to the number of ways in which a histogram with new_sum bins and a sum of squared bins equal to new_value can be formed
                        current_counts[new_sum][value_idx[new_value]]+=previous_count*combs[remaining][new_last]

        # writing the exact distribution for histograms with idx+1 bins and various sample sizes between 1 and upper_N
        for N in range(1, upper_N+1):
            distribution=dict()
            for current_value_idx in range(0, values_count):
                current=current_counts[N][current_value_idx]
                if skip_zeros is False or current>0:
                    distribution[possible_values[current_value_idx]]=current
            exact_distribution_path=join(exact_distributions_dir, HISTOGRAM_CONFIGURATION_PATTERN%(N, idx+1))
            if isfile(exact_distribution_path) is False:
                values_counts=distribution
                scaling=combs[upper_N][N]
                with open(exact_distribution_path, "w", encoding="utf8") as fo:
                    for k, v in sorted(values_counts.items()):
                        if v!=0:
                            fo.write(str(k)+" "+str(v//scaling)+"\n")

def main():

    parser=argparse.ArgumentParser()
    parser.add_argument("-N", type=int, help="The upper sample size.", default=100)
    parser.add_argument("-n", type=int, help="The upper bins count.", default=100)

    args=parser.parse_args()

    upper_N=args.N
    upper_n=args.n

    cache_exact_distribution_up_to_both_ns(upper_N=upper_N, upper_n=upper_n, show_status=True, use_cached=True, skip_zeros=True)
    
if __name__=="__main__":
    main()
