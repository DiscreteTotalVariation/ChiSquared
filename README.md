# Zero-disparity Distribution Synthesis: Fast Exact Calculation of Chi-Squared Statistic Distribution for Discrete Uniform Histograms

This is the repository with the source code and experiments that were conducted for the calculation of the exact distribution of the chi-squared statistic for histograms with uniform distribution proposed by [Nikola Banić](https://scholar.google.com/citations?user=QSH8C_QAAAAJ&hl=en) and [Neven Elezović](https://scholar.google.com/citations?user=MlXwbFIAAAAJ&hl=en) in the paper [Zero-disparity Distribution Synthesis: Fast Exact Calculation of Chi-Squared Statistic Distribution for Discrete Uniform Histograms](https://arxiv.org/abs/2506.23416).

## Data

All of the download links in this section point directly to files; **selecting a link starts the download immediately** rather than opening a description page. The approximate download sizes are stated next to each link so that the cost of a download is known in advance. Note that all of these archives are compressed, so the extracted data occupies considerably more space than the figures below.

### Sample files

To make it possible to inspect the format and the content of the data without downloading anything, a few small example files are included directly in this repository, under [`samples`](samples). They correspond to two of the cases shown in the figures of the paper, namely `N = 20`, `n = 4` and `N = 10`, `n = 20` (see Fig. 1), which are small enough to be read at a glance:

* [`samples/exact_distributions`](samples/exact_distributions) contains examples of the exact distributions in the form of counts. Each line holds two numbers, `integer_form count`, where `integer_form` is the integer form `s` of the chi-squared value (see Eq. (4) in the paper) and `count` is the number of histograms attaining it. These counts are exact arbitrary-precision integers, which is why the full archives of counts are so large.
* [`samples/exact_distributions_fp`](samples/exact_distributions_fp) contains the same distributions as floating-point probabilities. Each line holds three numbers, `integer_form value exponent_for_the_base_10`, where the probability of `integer_form` is given as `value` × 10^`exponent_for_the_base_10`. This is the form that is convenient for numerics, plotting, and any approximation. Probabilities far below the resolution of ordinary double-precision arithmetic (roughly 10^-16) are negligible and can be treated as zero for all practical purposes.
* [`samples/exact_distributions_real_fp`](samples/exact_distributions_real_fp) contains the same probabilities written as ordinary floating-point numbers. Each line holds two numbers, `integer_form realfp`, where `realfp` is the probability of `integer_form` as a single floating-point number (for example, `100 0.01067086943658069`) rather than split into a value and a base-10 exponent. This is the most convenient form for direct numerical use.

The files named `N_10_n_20.txt` in the three sample directories describe the same case (`N = 10`, `n = 20`) in the different forms and can be compared directly: for example, the count 670442572800 out of the 20^10 = 10240000000000 possible assignments corresponds to the probability 6.54729075 × 10^-2, which is written as `6.54729075 -2` in the `exact_distributions_fp` form and as `0.0654729075` in the `exact_distributions_real_fp` form.

### Full downloads

The recommended download is the archive of floating-point probabilities, which is by far the smallest of the bulk archives and is sufficient for numerics, plotting, and approximation:

* The probabilities for a grid of values of `N` and `n` are available [here](https://www.dropbox.com/scl/fi/7p8wvy4neps2vqe26qxup/exact_distributions_fp.7z?rlkey=yqvsc3dy1oclsh7dnrwaa9mtk&st=3khh8w2z&dl=1) (about 1.1 GB). The format is the one described for `samples/exact_distributions_fp` above. This archive is not used by any of the scripts; it is provided for direct use in numerical work and inspection.
* The same probabilities written as ordinary floating-point numbers are available [here](https://www.dropbox.com/scl/fi/ki9tcipu05kwz2i9qudly/exact_distributions_real_fp.7z?rlkey=712j0oo0n6k9fpv0i7m6wyhj3&st=i0wc8zs5&dl=1) (about 490 MB). This archive covers exactly the same grid as the one above, but each line holds two numbers, `integer_form realfp`, where `realfp` is the probability of `integer_form` written directly as a single floating-point number (for example, `100 0.01067086943658069`) rather than split into a value and a base-10 exponent. This is the most convenient form for direct numerical use. Like the archive above, it is not used by any of the scripts.

The grid covered by the bulk archives below is the full combination of `N` taking the multiples of 10 from 10 to 500 with `n` taking every integer from 2 to 20 and then the multiples of 10 from 30 to 500, which amounts to 3350 `(N, n)` pairs in total.

If the exact integer counts themselves are required, they are available as a larger archive:

* The exact distributions for the grid (storing counts) are available [here](https://www.dropbox.com/scl/fi/1tbyqrrbl53ivooii8fw1/exact_distributions.7z?rlkey=3wtvdslwnvywo4bq3v4b5de2v&st=tes9gauv&dl=1) (about 18 GB). The format is the one described for `samples/exact_distributions` above. This is the archive that the scripts use once it is extracted into `chi_data`.

The complete exact distributions for the entire range of `N` and `n`, rather than just the grid, can be downloaded by using [this script](download_exact_distributions.sh). **Warning.** This notice must be read in full before the script is executed. The script retrieves a prohibitively large volume of data: the archive is split into eight parts, of which the first seven are about 102 GB each and the last is about 95 GB, so the download totals roughly 810 GB, and upon extraction the full set of exact distributions stores the raw counts, compounding the storage demand still further. Execution therefore carries a substantial risk of exhausting all available disk space, monopolising network bandwidth for a protracted period, and rendering the system inoperable. For these reasons its use is emphatically discouraged and should be regarded as a measure of last resort. The script ought to be executed only when the complete dataset is unequivocally required, when its consequences are fully understood, and when ample storage and bandwidth have been provisioned in advance. In all but these exceptional circumstances, the grid archives above must be used instead. As a safeguard, the script does not download anything until the warning has been acknowledged interactively by typing `I UNDERSTAND`; this prompt can be bypassed for non-interactive use by passing the `--yes` flag.

Some precalculated exact and approximated CDFs are available [here](https://www.dropbox.com/scl/fi/64n5yhm1zir1yq1tmj1iz/chi_cdf_exact_vs_approximated.7z?rlkey=bxug2xn5whrtaoaz6q7561ngb&dl=1) (about 7.5 GB).

Some precalculated values of the Kolmogorov-Smirnov statistic are available [here](https://www.dropbox.com/scl/fi/f619rkeesr1tmb6s5wfnu/chi_kss.7z?rlkey=vcgvwsip53vxcvlhx8j5oy61s&dl=1) (about 2.7 MB).

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
