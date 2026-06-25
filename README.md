# Zero-disparity Distribution Synthesis: Fast Exact Calculation of Chi-Squared Statistic Distribution for Discrete Uniform Histograms

This is the repository with the source code and experiments that were conducted for the calculation of the exact distribution of the chi-squared statistic for histograms with uniform distribution proposed by [Nikola Banić](https://scholar.google.com/citations?user=QSH8C_QAAAAJ&hl=en) and [Neven Elezović](https://scholar.google.com/citations?user=MlXwbFIAAAAJ&hl=en) in the paper [Zero-disparity Distribution Synthesis: Fast Exact Calculation of Chi-Squared Statistic Distribution for Discrete Uniform Histograms](https://arxiv.org/abs/2506.23416).

## Data

Some precalculated exact distributions can be downloaded by using [this script](download_exact_distributions.sh). **Warning.** This notice must be read in full before the script is executed. The script retrieves a prohibitively large volume of data: all but the last of the archive parts are 100 GB each, so the download approaches one terabyte in aggregate, and upon extraction the full set of exact distributions stores the raw counts, compounding the storage demand still further. Execution therefore carries a substantial risk of exhausting all available disk space, monopolising network bandwidth for a protracted period, and rendering the system inoperable. For these reasons its use is emphatically discouraged and should be regarded as a measure of last resort. The script ought to be executed only when the complete dataset is unequivocally required, when its consequences are fully understood, and when ample storage and bandwidth have been provisioned in advance. In all but these exceptional circumstances, the two substantially smaller archives described below must be used instead.

Instead, it is recommended to use the following two smaller archives, which cover only a grid of values of `N` and `n` rather than all of them. The grid is the full combination of `N` taking the multiples of 10 from 10 to 500 with `n` taking every integer from 2 to 20 and then the multiples of 10 from 30 to 500, which amounts to 3350 `(N, n)` pairs in total:

* The exact distributions for that grid (storing counts) are available [here](https://www.dropbox.com/scl/fi/1tbyqrrbl53ivooii8fw1/exact_distributions.7z?rlkey=3wtvdslwnvywo4bq3v4b5de2v&st=tes9gauv&dl=1).
* The corresponding probabilities for a chi-squared value are available [here](https://www.dropbox.com/scl/fi/7p8wvy4neps2vqe26qxup/exact_distributions_fp.7z?rlkey=yqvsc3dy1oclsh7dnrwaa9mtk&st=3khh8w2z&dl=1). In this archive each line holds three numbers in the form `integer_form value exponent_for_the_base_10`, where `integer_form` is the integer form `s` of the chi-squared value (see Eq. (4) in the paper) and the corresponding probability is given as `value` × 10^`exponent_for_the_base_10`. The same grid of values of `N` and `n` is covered as in the archive above. This archive is provided only for convenience and inspection; it is not used by any of the scripts.

Some precalculated exact and approximated CDFs are available [here](https://www.dropbox.com/scl/fi/64n5yhm1zir1yq1tmj1iz/chi_cdf_exact_vs_approximated.7z?rlkey=bxug2xn5whrtaoaz6q7561ngb&dl=1).

Some precalculated values of the Kolmogorov-Smirnov statistic are available [here](https://www.dropbox.com/scl/fi/f619rkeesr1tmb6s5wfnu/chi_kss.7z?rlkey=vcgvwsip53vxcvlhx8j5oy61s&dl=1).

Once the [7z](https://www.7-zip.org/) archives with the precalculated data are downloaded, it is required to extract them inside of the `chi_data` directory for the code mentioned below to be able to use the data.

## Code

The code is written for Python 3 and it can be used to recreate the main results presented in the paper. It can also be changed to check the behavior in conditions not mentioned in the paper.

### Exact distribution

The file [calculate_exact_distributions.py](calculate_exact_distributions.py) contains the script for the calculating the exact distribution for a given range of values of `N` and `n`. By default, the script will store the results in the directory `chi_data/chi_exact_distributions`. In addition to that, it will also store the possible sums of squared bin values for certain values of `N` and `n`.

### Exact and approximated CDFs and the Kolmogorov-Smirnov statistic values

The file [calculate_cdfs_and_ks_statistic_values.py](calculate_cdfs_and_ks_statistic_values.py) contains the script for calculating the exact and the approximated values of CDF as well as calculating the values of the Kolmogorov-Smirnov statistic values. The calculations are done for the values of `N` and `n` that have a corresponding exact distribution in the directory `chi_data/chi_exact_distributions`. The results will be stored in the directories `chi_data/chi_cdf_exact_vs_approximated` and `chi_data/chi_kss`.

### Plots

The file [create_plots.py](create_plots.py) contains the scripts for generating the plots used in the paper. It relies on the input data in the `chi_data` directory obtained either by downloading the precalcualted data or by generating it from scratch with the provided code, and it stores the plots in the `img` directory.

## Authors

* **[Nikola Banić](https://scholar.google.com/citations?user=QSH8C_QAAAAJ&hl=en)**
* **[Neven Elezović](https://scholar.google.com/citations?user=MlXwbFIAAAAJ&hl=en)**

See also the list of [contributors](https://github.com/DiscreteTotalVariation/ChiSquared/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
