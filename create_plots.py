# -*- coding: utf-8 -*-    

import numpy as np
import cv2

from os import makedirs
from os.path import isfile, join
from tqdm import tqdm

import matplotlib
import matplotlib.pyplot as plt

from scipy.stats import chi2

ROOT_DIR="chi_data"

EXACT_DISTRIBUTIONS_DIR=join(ROOT_DIR, "chi_exact_distributions")

CDF_EXACT_VS_APPROXIMATED_DIR=join(ROOT_DIR, "chi_cdf_exact_vs_approximated")

KOLMOGOROV_SMIRNOV_STATISTIC_DIR=join(ROOT_DIR, "chi_kss")

HISTOGRAM_CONFIGURATION_PATTERN="N_%d_n_%d.txt"

HISTOGRAM_PLOT_CONFIGURATION_PATTERN="N_%d_n_%d.png"

LOG_P_TEMPLATE="log_p_N_%d_n_%d.png"

PLOTS_OUTPUT_DIR="img"

KSS_RATIOS_PLOT_PATH="ratios.png"

GOOD_APPROXIMATION_BORDER_PLOT="good_approximation_border.png"

BORDER_P_TEMPLATE="border_p_N_%d_n_%d.png"

def get_true_chi_statistic_value(integer_value, N, n):
    return integer_value*n/N-N

def crop_image(img):
    gray=cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mask=gray<255
    coords=np.argwhere(mask)
    y0, x0=coords.min(axis=0)
    y1, x1=coords.max(axis=0)+1
    cropped=img[y0:y1, x0:x1, :]
    return cropped

def create_cdf_comparison_plots(show_status=True):
    input_dir=CDF_EXACT_VS_APPROXIMATED_DIR
    output_dir=PLOTS_OUTPUT_DIR

    figsize=(8, 5)
    tick_fontsize=20
    label_fontsize=30
    linewidth=2
        
    dpi=None
    format="png"

    makedirs(output_dir, exist_ok=True)

    N_and_n_and_bounds_combinations=[(4, 4, (None, None)), (20, 4, (0.001, 0.999)), (10, 20, (0.001, 0.999)), (55, 10, (0.001, 0.999)), (200, 100, (0.001, 0.999)), (100, 200, (0.001, 0.999))]

    taker=N_and_n_and_bounds_combinations
    if show_status:
        taker=tqdm(taker)
    for N, n, bounds in taker:
    
        lower_bound, upper_bound=bounds

        output_path=join(output_dir, HISTOGRAM_PLOT_CONFIGURATION_PATTERN%(N, n))

        if isfile(output_path) is True:
            continue

        input_path=join(input_dir, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))

        d=np.loadtxt(input_path)
        if np.prod(d.shape)==3 and len(d.shape)==1:
            d=d[np.newaxis, :]
        if len(d.shape)==1:
            continue
        
        i=0
        start=None
        for x, cdf_exact, cdf_approx in d:
            if start is None and lower_bound is not None and cdf_exact>=lower_bound:
                start=i
            if upper_bound is not None and cdf_exact>=upper_bound and i>0:
                if start is None:
                    start=0
                d=d[start:i, :]
                break
            i+=1
        
        x=d[:, 0]
        cdf_exact=d[:, 1]
        cdf_approx=d[:, 2]

        error=np.abs(cdf_exact-cdf_approx)
        max_error_idx=np.argmax(error)
        x_max_error=x[max_error_idx]
        y_exact=cdf_exact[max_error_idx]
        y_approx=cdf_approx[max_error_idx]

        matplotlib.rc("font", size=tick_fontsize)
        fig=plt.figure(figsize=figsize)
        fig.add_subplot(1, 1, 1)

        plt.figure(figsize=figsize)
        plt.plot(x, cdf_exact, "o-", label="Exact CDF", color="blue", linewidth=linewidth)
        plt.plot(x, cdf_approx, "o--", label="Approximated CDF", color="red", linewidth=linewidth)

        plt.annotate(
            "", 
            xy=(x_max_error, y_exact), 
            xytext=(x_max_error, y_approx), 
            arrowprops=dict(arrowstyle="<->", color="green", lw=2)
        )

        plt.text(
            x_max_error+0.2, (y_exact + y_approx)/2,
            f"Maximum error = {error[max_error_idx]:.3f}",
            color="green",
            va="center"
        )

        plt.xlabel("$\\chi^2$ statistic", fontsize=label_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        plt.ylabel("CDF", fontsize=label_fontsize)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output_path, dpi=dpi, format=format)
        plt.close()

        img=cv2.imread(output_path, 6)
        
        cropped=crop_image(img)

        cv2.imwrite(output_path, cropped)

def create_logarithm_of_p_plots(show_status=True):

    N_and_n_combinations=[(10, 20), (55, 10), (200, 100), (100, 200)]

    output_dir=PLOTS_OUTPUT_DIR

    makedirs(output_dir, exist_ok=True)

    taker=N_and_n_combinations
    if show_status is True:
        taker=tqdm(taker)
    for N, n in taker:
        output_path=join(output_dir, LOG_P_TEMPLATE%(N, n))

        figsize=(13, 8)
        tick_fontsize=20
        label_fontsize=30

        dpi=300
        format="png"

        distribution=dict()
        input_path=join(EXACT_DISTRIBUTIONS_DIR, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
        with open(input_path, "r", encoding="utf8") as f:
            lines=[l.strip("\r\n") for l in f.readlines()]
        for l in lines:
            a, b=list(map(int, l.split()))
            distribution[a]=b

        p_values=[]
        chi_statistic_values=[]

        denominator=n**N
        for value, occurrence in distribution.items():
            if occurrence==0:
                continue

            chi_statistic=get_true_chi_statistic_value(integer_value=value, N=N, n=n)
            chi_statistic_values.append(chi_statistic)
            p_values.append(occurrence/denominator)

        matplotlib.rc("font", size=tick_fontsize)
        fig=plt.figure(figsize=figsize)
        fig.add_subplot(1, 1, 1)

        plt.figure(figsize=figsize)
        plt.plot(chi_statistic_values, p_values, "-", linewidth=2, label="p")
        plt.xlabel("$\\chi^2$ statistic", fontsize=label_fontsize)
        plt.ylabel("Probability", fontsize=label_fontsize)
        plt.yscale("log")
        plt.savefig(output_path, dpi=dpi, format=format)
        plt.close()

        img=cv2.imread(output_path, 6)
        cropped=crop_image(img)
        cv2.imwrite(output_path, cropped)

def create_kss_plots(show_status=True):
    input_dir=KOLMOGOROV_SMIRNOV_STATISTIC_DIR

    output_dir=PLOTS_OUTPUT_DIR

    makedirs(output_dir, exist_ok=True)

    combinations=[(list(range(2, 400+1)), list(range(1, 400+1)), join(output_dir, "kss_1_2.png")), (list(range(50, 400+1)), list(range(50, 400+1)), join(output_dir, "kss_50_50.png"))]

    taker=combinations
    if show_status is True:
        taker=tqdm(taker)

    for ns, Ns, output_path in taker:

        orientation="vertical"
        transpose=True

        reverse_y_axis=True

        figsize=(13, 8)
        tick_fontsize=20
        label_fontsize=30

        dpi=300
        format="png"

        errors=np.zeros((len(Ns), len(ns)))

        for Ni in range(len(Ns)):
            N=Ns[Ni]
            for ni in range(len(ns)):
                n=ns[ni]
                input_path=join(input_dir, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
                with open(input_path, "r", encoding="utf8") as f:
                    errors[Ni, ni]=float(f.read())
        if transpose is True:
            errors=errors.T

        matplotlib.rc("font", size=tick_fontsize)
        fig=plt.figure(figsize=figsize)
        fig.add_subplot(1, 1, 1)

        first_source=Ns
        second_source=ns
        if transpose is True:
            first_source, second_source=Ns, ns
        else:
            first_source, second_source=ns, Ns

        if reverse_y_axis is True:
            errors=errors[::-1, :]
            second_source=second_source[::-1]

        xticks=[i for i in range(len(first_source)) if first_source[i]%50==0]
        yticks=[i for i in range(len(second_source)) if second_source[i]%50==0]

        if xticks[0]!=0:
            xticks=[0]+xticks
        if reverse_y_axis is False:
            if yticks[0]!=0:
                yticks=[0]+yticks
        else:
            if yticks[-1]!=len(second_source)-1:
                yticks=[len(second_source)-1]+yticks

        xtick_labels=[first_source[i] for i in xticks]
        ytick_labels=[second_source[i] for i in yticks]

        plt.figure(figsize=figsize)
        plt.imshow(errors, cmap="turbo", interpolation="nearest")
        plt.colorbar(label="Kolmogorov-Smirnov statistic", orientation=orientation)
        if transpose is True:
            plt.xlabel("N", fontsize=label_fontsize)
            plt.ylabel("n", fontsize=label_fontsize)
        else:
            plt.xlabel("n", fontsize=label_fontsize)
            plt.ylabel("N", fontsize=label_fontsize)
        plt.xticks(xticks, xtick_labels)
        plt.yticks(yticks, ytick_labels)
        plt.savefig(output_path, dpi=dpi, format=format)
        plt.close()

        img=cv2.imread(output_path, 6)

        cropped=crop_image(img)

        cv2.imwrite(output_path, cropped)

def create_ratios_plot():
    input_dir=KOLMOGOROV_SMIRNOV_STATISTIC_DIR

    output_dir=PLOTS_OUTPUT_DIR

    makedirs(output_dir, exist_ok=True)

    output_path=join(output_dir, KSS_RATIOS_PLOT_PATH)

    ns=list(range(2, 100+1))
    ratios=[-2, 1, 2, 3, 4, 5]

    figsize=(13, 8)
    tick_fontsize=20
    label_fontsize=30

    threshold_plot_value=0.01

    dpi=300
    format="png"

    matplotlib.rc("font", size=tick_fontsize)
    fig=plt.figure(figsize=figsize)
    fig.add_subplot(1, 1, 1)

    xticks=[n for n in ns if n%10==0]

    if ns[0]%5!=0:
        xticks=[ns[0]]+xticks

    xtick_labels=xticks

    plt.figure(figsize=figsize)

    line_style="-"
    linewidth=3

    for ratio in ratios:
        x=[]
        y=[]
        for n in ns:
            if ratio>0:
                N=ratio*n
            else:
                N=n//-ratio
            if N==0:
                continue
            input_path=join(input_dir, "N_%d_n_%d.txt"%(N, n))
            if isfile(input_path) is False:
                continue
            x.append(n)
            with open(input_path, "r", encoding="utf8") as f:
                y.append(float(f.read()))

        if ratio>0:
            if N==n:
                plt.plot(x, y, line_style, label="$N=n$", linewidth=linewidth)
            else:
                plt.plot(x, y, line_style, label="$N="+str(ratio)+"\\cdot n$", linewidth=linewidth)
        else:
            plt.plot(x, y, line_style, label="$N=\\left\\lfloor \\frac{n}{"+str(-ratio)+"} \\right\\rfloor$", linewidth=linewidth)

    plt.plot(ns, [threshold_plot_value]*len(ns), "--", label=str(threshold_plot_value))

    plt.xlabel("n", fontsize=label_fontsize)
    plt.ylabel("Kolmogorov-Smirnov statistic", fontsize=label_fontsize)
    plt.xticks(xticks, xtick_labels)
    plt.legend()
    plt.savefig(output_path, dpi=dpi, format=format)
    plt.close()

    img=cv2.imread(output_path, 6)

    cropped=crop_image(img)

    cv2.imwrite(output_path, cropped)

def create_good_approximation_border_plot():
    output_dir=PLOTS_OUTPUT_DIR

    makedirs(output_dir, exist_ok=True)

    thresholds=[0.02, 0.05, 0.1]
    output_path=join(output_dir, GOOD_APPROXIMATION_BORDER_PLOT)

    ns=list(range(5, 250+1))

    upper_N=10000

    figsize=(13, 8)
    tick_fontsize=20
    label_fontsize=30
    line_style="-"
    linewidth=3

    dpi=300
    format="png"

    matplotlib.rc("font", size=tick_fontsize)
    fig=plt.figure(figsize=figsize)
    fig.add_subplot(1, 1, 1)

    for threshold in thresholds:
        plot_ns=[]
        plot_Ns=[]
        for n in ns:
            N=1
            while(N<=upper_N):
                if isfile(join(KOLMOGOROV_SMIRNOV_STATISTIC_DIR, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))) is True:
                    with open(join(KOLMOGOROV_SMIRNOV_STATISTIC_DIR, HISTOGRAM_CONFIGURATION_PATTERN%(N, n)), "r", encoding="utf8") as f:
                        ks=float(f.read())
                if ks<threshold:
                    break
                N+=1

            if N<=upper_N:
                plot_ns.append(n)
                plot_Ns.append(N)

        plt.plot(plot_ns, plot_Ns, line_style, label="K-S statistic below "+str(threshold), linewidth=linewidth)

    plt.xlabel("$n$", fontsize=label_fontsize)
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.ylabel("$N$", fontsize=label_fontsize)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, format=format)
    plt.close()

    img=cv2.imread(output_path, 6)

    cropped=crop_image(img)

    cv2.imwrite(output_path, cropped)

def create_border_probability_plot():
    N=55
    n=10

    df=n-1

    output_dir=PLOTS_OUTPUT_DIR

    threshold=10**-4

    offset=7
    offset2=3

    makedirs(output_dir, exist_ok=True)

    output_path=join(output_dir, BORDER_P_TEMPLATE%(N, n))

    figsize=(13, 8)
    tick_fontsize=20
    label_fontsize=30

    dpi=300
    format="png"

    distribution=dict()
    input_path=join(EXACT_DISTRIBUTIONS_DIR, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
    with open(input_path, "r", encoding="utf8") as f:
        lines=[l.strip("\r\n") for l in f.readlines()]
    for l in lines:
        a, b=list(map(int, l.split()))
        distribution[a]=b

    chi_statistic_values=[]
    exact_p_values=[]
    approximated_p_values=[]

    numerator=0
    denominator=n**N
    for value, occurrence in distribution.items():
        if occurrence==0:
            continue

        chi_statistic=get_true_chi_statistic_value(integer_value=value, N=N, n=n)

        exact_p_value=1-numerator/denominator
        numerator+=occurrence
        approximated_p_value=chi2.sf(chi_statistic, df)

        exact_p_values.append(exact_p_value)
        approximated_p_values.append(approximated_p_value)

        chi_statistic_values.append(chi_statistic)

    for i in range(len(exact_p_values)):
        if exact_p_values[i]<threshold:
            break

    i1=i-offset
    i2=i+offset

    ii1=None
    ii2=None
    for i in range(i1, i2+1):
        if (exact_p_values[i]<threshold)!=(approximated_p_values[i]<threshold):
            if ii1 is None:
                ii1=i
        else:
            if ii1 is not None and ii2 is None:
                ii2=i-1

    iii1=ii1-offset2
    iii2=ii2+offset2

    matplotlib.rc("font", size=tick_fontsize)
    fig=plt.figure(figsize=figsize)
    fig.add_subplot(1, 1, 1)

    plt.figure(figsize=figsize)
    plt.plot(chi_statistic_values[iii1:iii2+1], exact_p_values[iii1:iii2+1], "o-", linewidth=2, label="Exact $p$-value")
    plt.plot(chi_statistic_values[iii1:iii2+1], approximated_p_values[iii1:iii2+1], "o-", linewidth=2, color="green", label="Approximated $p$-value")
    plt.plot(chi_statistic_values[ii1:ii2+1], approximated_p_values[ii1:ii2+1], "o-", linewidth=2, color="red", label="Problematic approximated $p$-value")
    plt.plot(chi_statistic_values[iii1:iii2+1], [threshold]*(iii2-iii1+1), "--", linewidth=2, color="black", label=str(threshold))
    plt.xlabel("$\\chi^2$ statistic", fontsize=label_fontsize)
    plt.ylabel("$p$-value", fontsize=label_fontsize)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, format=format)
    plt.close()

    img=cv2.imread(output_path, 6)

    cropped=crop_image(img)

    cv2.imwrite(output_path, cropped)

def create_input_data_for_type_i_error_plots(show_status=True):

    thresholds=[0.05, 0.001, 0.0001, 0.00001]

    for ti, threshold in enumerate(thresholds):
        if show_status is True:
            print("%d / %d"%(ti+1, len(thresholds)))

        Ns=list(range(1, 750+1))
        n=10

        output_path="type_i_error_for_approximation_n_%d_a_%f.txt"%(n, threshold)

        previous_data=dict()

        if isfile(output_path) is True:
            with open(output_path, "r", encoding="utf8") as f:
                lines=[l.strip("\r\n") for l in f.readlines()]
            for l in lines:
                parts=l.split()
                N=int(parts[0])
                a=float(parts[1])
                b=float(parts[2])
                previous_data[N]=(a, b)

        input_dir=CDF_EXACT_VS_APPROXIMATED_DIR

        taker=Ns
        if show_status is True:
            taker=tqdm(taker)
        with open(output_path, "w", encoding="utf8") as fo:
            for N in taker:
                exact_type_i_error=None
                approximated_type_i_error=None
                if N in previous_data.keys():
                    exact_type_i_error, approximated_type_i_error=previous_data[N]
                else:

                    distribution=dict()
                    input_path=join(EXACT_DISTRIBUTIONS_DIR, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))
                    with open(input_path, "r", encoding="utf8") as f:
                        lines=[l.strip("\r\n") for l in f.readlines()]
                    for l in lines:
                        a, b=list(map(int, l.split()))
                        distribution[a]=b

                    distribution_integer_values=[k for k, v in sorted(distribution.items())]
                    distribution_counts=[v for k, v in sorted(distribution.items())]

                    input_path=join(input_dir, HISTOGRAM_CONFIGURATION_PATTERN%(N, n))

                    exact_p_values=[]
                    approximated_p_values=[]

                    d=np.loadtxt(input_path)
                    if np.prod(d.shape)==3 and len(d.shape)==1:
                        d=d[np.newaxis, :]
                    if len(d.shape)==1:
                        return

                    for i in range(d.shape[0]):
                        exact_p_value=1 if i==0 else 1-d[i-1, 1]

                        approximated_p_value=1-d[i, 2]

                        exact_p_values.append(exact_p_value)
                        approximated_p_values.append(approximated_p_value)

                    problematic_count=0
                    after_problematic_count=0
                    problematic_statistic_count=0
                    after_problematic_statistic_count=0
                    interesting_statistic_values=[]
                    interesting_probabilities=[]

                    exact_rejected_count=0
                    approximated_rejected_count=0

                    for integer_value, current_count, exact_p_value, approximated_p_value in zip(distribution_integer_values, distribution_counts, exact_p_values, approximated_p_values):
                        if exact_p_value<threshold:
                            exact_rejected_count+=current_count
                        if approximated_p_value<=threshold:
                            approximated_rejected_count+=current_count
                        if (exact_p_value<threshold)!=(approximated_p_value<=threshold):
                            problematic_count+=current_count
                            problematic_statistic_count+=1
                            interesting_statistic_values.append(get_true_chi_statistic_value(integer_value=integer_value, N=N, n=n))
                            interesting_probabilities.append(current_count/(n**N))
                        elif problematic_count>0:
                            after_problematic_count+=current_count
                            after_problematic_statistic_count+=1
                            interesting_statistic_values.append(get_true_chi_statistic_value(integer_value=integer_value, N=N, n=n))
                            interesting_probabilities.append(current_count/(n**N))

                    exact_type_i_error=exact_rejected_count/(n**N)
                    approximated_type_i_error=approximated_rejected_count/(n**N)

                fo.write("%d %.16f %.16f\n"%(N, exact_type_i_error, approximated_type_i_error))
                fo.flush()

def create_type_i_error_plots(show_status=True):
    output_dir=PLOTS_OUTPUT_DIR

    makedirs(output_dir, exist_ok=True)

    thresholds=[0.05, 0.001, 0.0001, 0.00001]

    taker=thresholds
    if show_status is True:
        taker=tqdm(thresholds)
    for threshold in taker:
        n=10

        input_path="type_i_error_for_approximation_n_%d_a_%f.txt"%(n, threshold)
        output_path=join(output_dir, "type_i_error_for_approximation_n_%d_a_%f.png"%(n, threshold))

        lower_N=50

        Ns=[]
        exact_type_i_errors=[]
        approximated_type_i_errors=[]

        with open(input_path, "r", encoding="utf8") as f:
            lines=[l.strip("\r\n") for l in f.readlines()]
        for l in lines:
            parts=l.split()
            N=int(parts[0])
            a=float(parts[1])
            b=float(parts[2])
            if lower_N is not None and N<lower_N:
                continue
            Ns.append(N)
            exact_type_i_errors.append(a)
            approximated_type_i_errors.append(b)

        figsize=(13, 8)
        tick_fontsize=20
        label_fontsize=30
        linewidth=2

        dpi=300
        format="png"

        matplotlib.rc("font", size=tick_fontsize)
        fig=plt.figure(figsize=figsize)
        fig.add_subplot(1, 1, 1)

        plt.plot(Ns, approximated_type_i_errors, "-", linewidth=linewidth, label="Approximated distribution")
        plt.plot(Ns, exact_type_i_errors, "-", linewidth=linewidth, label="Exact distribution")

        plt.xlabel("$N$", fontsize=label_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        plt.ylabel("Type I error rate", fontsize=label_fontsize)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output_path, dpi=dpi, format=format)
        plt.close()

        img=cv2.imread(output_path, 6)

        cropped=crop_image(img)

        cv2.imwrite(output_path, cropped)

def main():
    
    print("Creating the CDF comparison plots...")
    create_cdf_comparison_plots()

    print("Creating the plots for the logarithm of chi-squared statistic probability...")
    create_logarithm_of_p_plots()

    print("Creating the Kolmogorov-Smirnov statistic value plots...")
    create_kss_plots()

    print("Creating the Kolmogorov-Smirnov statistic values for points with various ratios between N and n...")
    create_ratios_plot()

    print("Creating the good approximation border plot...")
    create_good_approximation_border_plot()

    print("Creating the border probability plot...")
    create_border_probability_plot()

    print("Creating input data for type I error plots...")
    create_input_data_for_type_i_error_plots()

    print("Creating type I error plots...")
    create_type_i_error_plots()

    pass
    
if __name__=="__main__":
    main()
