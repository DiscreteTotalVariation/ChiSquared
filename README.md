# Zero-disparity Distribution Synthesis: Fast Exact Calculation of Chi-Squared Statistic Distribution for Discrete Uniform Histograms

This is the repository with the source code and experiments that were conducted for the calculation of the exact distribution of the chi-squared statistic for histograms with uniform distribution proposed by [Nikola Banić](https://scholar.google.com/citations?user=QSH8C_QAAAAJ&hl=en) and [Neven Elezović](https://scholar.google.com/citations?user=MlXwbFIAAAAJ&hl=en) in the paper [Zero-disparity Distribution Synthesis: Fast Exact Calculation of Chi-Squared Statistic Distribution for Discrete Uniform Histograms](https://arxiv.org/abs/2506.23416).

## Data

Some precalculated exact distributions are still being compressed at the moment of writing this text.

Some precalculated exact and approximated CDFs are available [here](https://www.dropbox.com/scl/fi/64n5yhm1zir1yq1tmj1iz/chi_cdf_exact_vs_approximated.7z?rlkey=bxug2xn5whrtaoaz6q7561ngb&dl=1).

Some precalculated values of the Kolmogorov-Smirnov statistic are available [here](https://www.dropbox.com/scl/fi/f619rkeesr1tmb6s5wfnu/chi_kss.7z?rlkey=vcgvwsip53vxcvlhx8j5oy61s&dl=1).

Once the [7z](https://www.7-zip.org/) archives with the precalculated data are downloaded, it is required to extract them inside of the `chi_data` directory for the code mentioned below to be able to use the data.

## Code

The code is written for Python 3 and it can be used to recreate the main results presented in the paper. It can also be changed to check the behavior in conditions not mentioned in the paper.

### Exact distribution

The file [calculate_exact_distributions.py](calculate_exact_distributions.py) contains the script for the calculating the exact distribution for a given range of values of `N` and `n`. By default, the script will store the results in the directory `chi_data/chi_exact_distributions`. In addition to that, it will also store the possible sums of squared bin values for certain values of `N` and `n`.

### Exact and approximated CDFs and the Kolmogorov-Smirnov statistic values

The file [calculate_cdfs_and_ks_statistic_values.py](calculate_cdfs_and_ks_statistic_values.py) contains the script for the calculating the exact and the approximated values of CDF as well as calculating the values of the Kolmogorov-Smirnov statistic values. The calculations are done for the values of `N` and `n` that have a corresponding exact distribution in the directory `chi_data/chi_exact_distributions`. The results will be stored in the directories `chi_data/chi_cdf_exact_vs_approximated` and `chi_data/chi_kss`.

## Authors

* **[Nikola Banić](https://scholar.google.com/citations?user=QSH8C_QAAAAJ&hl=en)**
* **[Neven Elezović](https://scholar.google.com/citations?user=MlXwbFIAAAAJ&hl=en)**

See also the list of [contributors](https://github.com/DiscreteTotalVariation/ChiSquared/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
