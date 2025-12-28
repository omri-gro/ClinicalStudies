import numpy as np
import scipy.stats as stats
# from scipy.sparse import coo_matrix
from joblib import Parallel, delayed


def col_totals(df, replace_columns=False):
    tots_dit = {
        'Total Imm': ['Promyelocyte', 'Myelocyte', 'Metamyelocyte'],
        'Total Neutrophil': ['Segmented Neutrophil', 'Band Neutrophil'],
        'Total Var Lym': ['Large Granular Lymphocyte', 'Atypical Lymphocyte', 'Aberrant Lymphocyte'],
        'Total Lym': ['Lymphocyte', 'Large Granular Lymphocyte', 'Atypical Lymphocyte', 'Aberrant Lymphocyte'],
        'LymRegLGL': ['Lymphocyte', 'Large Granular Lymphocyte'],
        'BlastImm': ['Blast', 'Promyelocyte', 'Myelocyte', 'Metamyelocyte']
    }
    for target_class, source_classes in tots_dit.items():
        for mthd in ['OMR', 'DSS']:
            target_col = f'{target_class} {mthd}'
            source_cols = [f'{source_class} {mthd}' for source_class in source_classes]
            if replace_columns or target_col not in df.columns:
                # If the target column doesn't exist, create it as the sum of source_cols
                df[target_col] = df[source_cols].sum(axis=1)
            else:
                # If the target column exists, fill NaN values with the sum of source_cols
                df[target_col] = df[target_col].fillna(df[source_cols].sum(axis=1))
    return df


def deming_regression_pe(x, y, lambda_=1):
    x, y = np.array(x), np.array(y)
    x_mean, y_mean = np.mean(x), np.mean(y)
    Sxy = np.sum((x - x_mean) * (y - y_mean))
    Sxx = np.sum((x - x_mean) ** 2)
    Syy = np.sum((y - y_mean) ** 2)
    delta = Syy - lambda_ * Sxx
    slope = (delta + np.sqrt(delta ** 2 + 4 * lambda_ * Sxy ** 2)) / (2 * Sxy)
    intercept = y_mean - slope * x_mean
    return slope, intercept


def bootstrap_deming_ci(x_series, y_series, lambda_=1, n_bootstrap=1000, ci=95):
    x = np.array(x_series)
    y = np.array(y_series)
    slopes, intercepts = [], []
    n = len(x)
    for _ in range(n_bootstrap):
        indices = np.random.choice(range(n), size=n, replace=True)
        x_sample = x[indices]
        y_sample = y[indices]
        slope, intercept = deming_regression_pe(x_sample, y_sample, lambda_)
        slopes.append(slope)
        intercepts.append(intercept)
    # Calculate confidence intervals
    slope_lower = np.percentile(slopes, (100 - ci) / 2)
    slope_upper = np.percentile(slopes, (100 + ci) / 2)
    intercept_lower = np.percentile(intercepts, (100 - ci) / 2)
    intercept_upper = np.percentile(intercepts, (100 + ci) / 2)
    return (slope_lower, slope_upper), (intercept_lower, intercept_upper)


def deming_regression(x, y, lambda_=1, n_bootstrap=1000, ci=95):
    slope, intercept = deming_regression_pe(x, y, lambda_)
    slope_ci, intercept_ci = bootstrap_deming_ci(x, y, lambda_, n_bootstrap, ci)
    return {
        "regression method": f"Deming - lambda={lambda_}",
        "CI method": "Bootstrap",
        "iterations": n_bootstrap,
        "slope": slope,
        "intercept": intercept,
        "slope_ci_bottom": slope_ci[0],
        "slope_ci_top": slope_ci[1],
        "intercept_ci_bottom": intercept_ci[0],
        "intercept_ci_top": intercept_ci[1]
    }


def u_test_ci(data, confidence):
    sorted_data = np.sort(data)
    n = len(sorted_data)
    alpha = 1 - (confidence / 100)
    z = stats.norm.ppf(1 - alpha / 2)
    k_lower = int((n / 2) - z * np.sqrt(n) / 2)
    k_upper = int((n / 2) + z * np.sqrt(n) / 2)
    k_lower = max(0, k_lower)
    k_upper = min(n - 1, k_upper)
    ci_low = sorted_data[k_lower]
    ci_high = sorted_data[k_upper]
    return (ci_low, ci_high)


def bootstrap_iteration(precomputed_slopes, indices, x, y):
    # Performs a single bootstrap iteration
    # Create a set of selected pairs where both points are in the bootstrapped sample
    indices_set = set(indices)
    selected_slopes = [
        precomputed_slopes[(i, j)]
        for (i, j) in precomputed_slopes
        if i in indices_set and j in indices_set
    ]

    if selected_slopes:  # Ensure non-empty list of slopes
        boot_slope = np.median(selected_slopes)
        boot_intercept = np.median(y[indices] - boot_slope * x[indices])
        return boot_slope, boot_intercept
    return None, None


def passing_bablok_regression(x, y, confidence=95, ci_method="U test", n_bootstrap=500, vectorize=True, n_jobs=-1,
                              advanced_slopes_filtering=True, remove_slopes_smaller_than_minus_1=True, tolerance=1e-6):
    # Ensure the input arrays are numpy arrays
    x = np.array(x)
    y = np.array(y)
    n = len(x)

    # Calculate pairwise slopes
    if vectorize:  # more efficient slopes calculation
        # Get indices for upper triangular matrix excluding the diagonal
        i, j = np.triu_indices(n, k=1)
        if advanced_slopes_filtering:  # in some cases slopes close to -1 are also removed, not appropriate here
            # Filter out pairs where x differences are near-zero or y values are equal
            #  valid_indices = (np.abs(x[i] - x[j]) > tolerance) & (y[i] != y[j])
            valid_indices = (x[i] != x[j]) & (y[i] != y[j])
        else:
            # Filter out cases where the x values are equal to avoid division by zero
            valid_indices = x[i] != x[j]
        i, j = i[valid_indices], j[valid_indices]
        slopes = (y[j] - y[i]) / (x[j] - x[i])
        # index_pairs = list(zip(i, j))
        if ci_method == "U test":
            intercepts = y[i] - slopes * x[i]
        elif ci_method == "Bootstrap":
            # slopes_matrix = coo_matrix((slopes, (i, j)), shape=(n, n)).tocsr()
            # for ease of access, save slopes as dictionary with (i, j) as keys
            slopes_dict = {(i, j): slope for i, j, slope in zip(i, j, slopes)}

    else:  # Will not work with optimized bootstrapping method
        slopes = []
        intercepts = []
        index_pairs = []  # To be used for bootstrap option, precomputing pairwise slopes
        for i in range(n - 1):
            for j in range(i + 1, n):
                if x[i] != x[j]:  # Avoid division by zero
                    slope = (y[j] - y[i]) / (x[j] - x[i])
                    slopes.append(slope)
                    index_pairs.append((i, j))
                    if ci_method == "U test":
                        intercepts.append(y[i] - slope * x[i])

    # Calculate slope and intercept
    slope = np.median(slopes)
    intercept = np.median(y - slope * x)

    # Calculate confidence intervals for slope and intercept
    if ci_method == "U test":
        # We use the Mann-Whitney U test to determine the confidence interval for the slope
        slope_ci = u_test_ci(slopes, confidence)
        intercept_ci = u_test_ci(intercepts, confidence)

    elif ci_method == "Bootstrap":
        # choose indices for all bootstrap iterations in advance
        bootstrap_samples = [np.random.choice(range(n), n, replace=True) for _ in range(n_bootstrap)]
        if n_jobs > 1 and isinstance(n_jobs, int):  # Parallel bootstrapping
            bootstrap_results = Parallel(n_jobs=n_jobs)(delayed(bootstrap_iteration)(
                slopes_dict, indices, x, y
            ) for indices in bootstrap_samples)

            # Filter out empty results
            bootstrap_slopes = [res[0] for res in bootstrap_results if res[0] is not None]
            bootstrap_intercepts = [res[1] for res in bootstrap_results if res[1] is not None]

        else:
            bootstrap_slopes = []
            bootstrap_intercepts = []
            bootstrap_samples = [np.random.choice(range(n), n, replace=True) for _ in range(n_bootstrap)]
            for indices in bootstrap_samples:
                # Sample with replacement
                boot_result = bootstrap_iteration(slopes_dict, indices, x, y)
                if boot_result[0] is not None:
                    bootstrap_slopes.append(boot_result[0])
                    bootstrap_intercepts.append(boot_result[1])

        # Calculate confidence intervals from bootstrap distribution
        lower_bound = (100 - confidence) / 2
        upper_bound = 100 - lower_bound
        slope_ci = (np.percentile(bootstrap_slopes, lower_bound), np.percentile(bootstrap_slopes, upper_bound))
        intercept_ci = (np.percentile(bootstrap_intercepts, lower_bound), np.percentile(bootstrap_intercepts, upper_bound))

    else:
        raise ValueError(f"Invalid option {ci_method} for calculation of CI")

    return {
        "regression method": "Passing-Bablok",
        "CI method": ci_method,
        "iterations": n_bootstrap,
        "slope": slope,
        "intercept": intercept,
        "slope_ci_bottom": slope_ci[0],
        "slope_ci_top": slope_ci[1],
        "intercept_ci_bottom": intercept_ci[0],
        "intercept_ci_top": intercept_ci[1]
    }


def regression_comp(x, y, n_bootstrap=1000, ci=95, reg_method="deming",
                    lambda_=1,
                    ci_method="Bootstrap", n_jobs=-1,
                    res_str=True):
    """
    Calculate regression for method comparison, getting all info on the regression.
    Parameters:
        - x, y: input array_like
        - reg_method: "deming" or "passing"
        - lambda_: variations ratio, only for Deming
        - ci_method: "U test" or "Bootstrap", only for Passing
        - n_jobs: number of cores for parallel processing, only for Passing
    Returns:
        Dict
    """
    # Ensure the input arrays are numpy arrays
    x = np.array(x)
    y = np.array(y)
    n = len(x)


    if reg_method == "deming":
        ci_method = "Bootstrap"
        corr = stats.pearsonr(x, y)
        dem_dict = deming_regression(x, y, lambda_, n_bootstrap, ci)

    elif reg_method == "passing":
        corr = stats.spearmanr(x, y)
        dem_dict = passing_bablok_regression(x, y, confidence=ci, ci_method=ci_method, n_bootstrap=n_bootstrap, n_jobs=n_jobs)

    else:
        raise ValueError(f"Invalid option {ci_method} for calculation of CI")

    add_data = {"reg_method": reg_method,
                "ci_method": ci_method,
                "iterations": n_bootstrap,
                "N": n,
                "correlation_coefficient": corr[0]}
    if res_str:
        slope_str = f"{dem_dict['slope']:.2f}\n({dem_dict['slope_ci_bottom']:.2f}-{dem_dict['slope_ci_top']:.2f})"
        int_str = f"{dem_dict['intercept']:.2f}\n({dem_dict['intercept_ci_bottom']:.2f}-{dem_dict['intercept_ci_top']:.2f})"
        add_data.update({"slope_str": slope_str,
                         "intercept_str": int_str})
    dem_dict.update(add_data)
    return dem_dict