import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from joblib import Parallel, delayed


def validate_regession_input(x, y):
    """
    Validate inputs are appropriate for regression calculation.
    Args:
        x (array-like): Independent variable.
        y (array-like): Dependent variable.

    Raises:
        ValueError: If x and y are invalid for regression.
    """
    x = np.array(x)
    y = np.array(y)

    # Check if x and y are of the same length
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")
    # Check if x and y are not empty
    if len(x) == 0:
        raise ValueError("x and y cannot be empty.")
    # Check if all values in x and y are identical
    if len(np.unique(x)) == 1 and len(np.unique(y)) == 1:
        raise ValueError("x and y cannot have all identical values.")

    return x, y  # Optionally return numpy arrays for consistency


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


def precompute_pairwise_slopes(x, y):
    # input two array-like of same length (n)
    # output dict matching pairs of indices (tuple) with slope (float), (n^2 - doubles) elements in dict
    x = np.array(x)
    y = np.array(y)
    i, j = np.triu_indices(len(x), k=1)
    valid_indices = x[i] != x[j]
    i, j = i[valid_indices], j[valid_indices]
    slopes = (y[j] - y[i]) / (x[j] - x[i])
    slopes_dict = {(i, j): slope for i, j, slope in zip(i, j, slopes)}
    return slopes_dict

# def passing_bablok_regression(x, y, confidence=95, n_bootstrap=200):
#     # Ensure the input arrays are numpy arrays
#     x = np.array(x)
#     y = np.array(y)
#     n = len(x)
#
#     # Precompute all pairwise slopes and their corresponding index pairs
#     slopes_dict = precompute_pairwise_slopes(x, y)
#
#     # Calculate the median slope and intercept from all pairwise slopes
#     slopes_list = list(slopes_dict.values())
#     slope = np.median(slopes_list)
#     intercept = np.median(y - slope * x)
#
#
#     bootstrap_slopes = []
#     bootstrap_intercepts = []
#     # create list of n_boostrap items, each is a ndarray of length n (containing independent indices of datapoint)
#     bootstrap_samples = [np.random.choice(range(n), n, replace=True) for _ in range(n_bootstrap)]
#     for indices in bootstrap_samples:
#         boot_result = bootstrap_iteration(slopes_dict, indices, x, y)
#         if boot_result[0] is not None:
#             bootstrap_slopes.append(boot_result[0])
#             bootstrap_intercepts.append(boot_result[1])
#
#     # Calculate confidence intervals from bootstrap distribution
#     lower_bound = (100 - confidence) / 2
#     upper_bound = 100 - lower_bound
#     slope_ci = (np.percentile(bootstrap_slopes, lower_bound), np.percentile(bootstrap_slopes, upper_bound))
#     intercept_ci = (np.percentile(bootstrap_intercepts, lower_bound), np.percentile(bootstrap_intercepts, upper_bound))
#
#     return {
#         "regression method": "Passing-Bablok",
#         "iterations": n_bootstrap,
#         "slope": slope,
#         "intercept": intercept,
#         "slope_ci_bottom": slope_ci[0],
#         "slope_ci_top": slope_ci[1],
#         "intercept_ci_bottom": intercept_ci[0],
#         "intercept_ci_top": intercept_ci[1]
#     }


def diff_method_pass_bab(x, y, confidence=95, n_bootstrap=1000):
    x = np.array(x)
    y = np.array(y)
    points = np.stack((x, y), axis=1)
    n = len(x)

    unique_points, unique_indices = np.unique(points, axis=0, return_inverse=True)
    x_coords = unique_points[:, 0]
    y_coords = unique_points[:, 1]

    # dx = x[:, None] - x  # Difference in x-coordinates
    # dy = y[:, None] - y  # Difference in y-coordinates
    dx = x_coords[:, None] - x_coords  # Difference in x-coordinates
    dy = y_coords[:, None] - y_coords  # Difference in y-coordinates
    with np.errstate(divide='ignore', invalid='ignore'):  # Handle division by zero
        slopes = np.divide(dy, dx)  # Element-wise division
        slopes[dx == 0] = np.NaN  # Assign infinity to vertical line

    bootstrap_slopes = []
    bootstrap_intercepts = []
    bootstrap_resampled_indices = [np.random.choice(unique_indices, n, replace=True) for _ in range(n_bootstrap)]
    for resampled_indices in bootstrap_resampled_indices:
        # convert this section to a function
        i_indices, j_indices = np.meshgrid(resampled_indices, resampled_indices, indexing='ij')
        resampled_slopes = slopes[i_indices, j_indices]
        median_slope = np.nanmedian(resampled_slopes)
        median_intercept = np.nanmedian(y[resampled_indices] - median_slope * x[resampled_indices])
        bootstrap_slopes.append(median_slope)
        bootstrap_intercepts.append(median_intercept)

    # Calculate confidence intervals from bootstrap distribution
    lower_bound = (100 - confidence) / 2
    upper_bound = 100 - lower_bound
    slope_ci = (np.percentile(bootstrap_slopes, lower_bound), np.percentile(bootstrap_slopes, upper_bound))
    intercept_ci = (np.percentile(bootstrap_intercepts, lower_bound), np.percentile(bootstrap_intercepts, upper_bound))

    return slope_ci, intercept_ci


def plot_ver_reg(raw_data, x_var, y_var, cls_var=None, reg_ser=None, reg_ci=False):
    fig, ax = plt.subplots()
    if cls_var in raw_data.columns:
        groups = raw_data.groupby(cls_var)
        for name, group in groups:
            ax.plot(group[x_var], group[y_var], marker='o', linestyle='', markersize=4, label=name)
    else:
        ax.plot(raw_data[x_var], raw_data[y_var], marker='o', linestyle='', markersize=4)
    if isinstance(reg_ser, pd.Series) or isinstance(reg_ser, dict):
        # Plot the regression line
        x = raw_data[x_var]
        x_range = np.linspace(min(x) * 0.9, max(x) * 1.1, 100)
        y_line = reg_ser['slope'] * x_range + reg_ser['intercept']
        ax.plot(x_range, y_line, label="Deming Regression Line", color="red")
        #
        # Plot the confidence interval
        if reg_ci:
            y_low = reg_ser['slope_ci_bottom'] * x_range + reg_ser['intercept_ci_bottom']
            y_top = reg_ser['slope_ci_top'] * x_range + reg_ser['intercept_ci_top']
            ax.fill_between(x_range, y_low, y_top, color="pink", alpha=0.3, label="95% Confidence Interval")
    ax.set_xlabel(x_var)
    ax.set_ylabel(y_var)
    # tick_spacing = 5
    # ax.xaxis.set_major_locator(MultipleLocator(tick_spacing))
    # ax.yaxis.set_major_locator(MultipleLocator(tick_spacing))
    ax.legend()
    ax.grid(True)
    plt.show()


def plot_ver_reg_dbl(raw_data, x_var, y_var, meas1_2=[], cls_var=None, reg_ser=None, reg_ci=False):
    fig, ax = plt.subplots()
    if cls_var in raw_data.columns:
        groups = raw_data.groupby(cls_var)
        for name, group in groups:
            ax.plot(group[x_var], group[y_var], marker='o', linestyle='', markersize=4, label=name)
    else:
        ax.plot(raw_data[x_var], raw_data[y_var], marker='o', linestyle='', markersize=4)
    if isinstance(reg_ser, pd.Series) or isinstance(reg_ser, dict):
        # Plot the regression line
        x = raw_data[x_var]
        x_range = np.linspace(min(x) * 0.9, max(x) * 1.1, 100)
        y_line = reg_ser['slope'] * x_range + reg_ser['intercept']
        ax.plot(x_range, y_line, label="Deming Regression Line", color="red")
        #
        # Plot the confidence interval
        if reg_ci:
            y_low = reg_ser['slope_ci_bottom'] * x_range + reg_ser['intercept_ci_bottom']
            y_top = reg_ser['slope_ci_top'] * x_range + reg_ser['intercept_ci_top']
            ax.fill_between(x_range, y_low, y_top, color="pink", alpha=0.3, label="95% Confidence Interval")
    if len(meas1_2) == 2:
        n_sample = len(raw_data[x_var])
        meas1 = raw_data[meas1_2[0]]
        meas2 = raw_data[meas1_2[1]]
        for i in np.arange(n_sample):
            ax.plot([meas1[i], meas2[i]], [raw_data[y_var][i], raw_data[y_var][i]], 'gray', alpha=0.7)
    ax.set_xlabel(x_var)
    ax.set_ylabel(y_var)
    # tick_spacing = 5
    # ax.xaxis.set_major_locator(MultipleLocator(tick_spacing))
    # ax.yaxis.set_major_locator(MultipleLocator(tick_spacing))
    ax.legend()
    ax.grid(True)
    plt.show()



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
    try:
        # Validate inputs
        x, y = validate_regession_input(x, y)
        slope, intercept = deming_regression_pe(x, y, lambda_)
        slope_ci, intercept_ci = bootstrap_deming_ci(x, y, lambda_, n_bootstrap, ci)
    except Exception as e:
        print(f"Error in Deming regression: {e}")
        slope = intercept = 0
        slope_ci = intercept_ci = (0, 0)
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
                              advanced_slopes_filtering=True, remove_slopes_smaller_than_minus_1=True, tolerance=1e-6,
                              dtype=np.float32):
    try:
        x, y = validate_regession_input(x, y)
        # Ensure the input arrays are numpy arrays
        x = np.array(x, dtype=dtype)
        y = np.array(y, dtype=dtype)
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
            try:
                lower_bound = (100 - confidence) / 2
                upper_bound = 100 - lower_bound
                slope_ci = (np.percentile(bootstrap_slopes, lower_bound), np.percentile(bootstrap_slopes, upper_bound))
                intercept_ci = (np.percentile(bootstrap_intercepts, lower_bound), np.percentile(bootstrap_intercepts, upper_bound))
            except IndexError:
                print("No slopes can be calculated")
                slope_ci = (0, 0)
                intercept_ci = (0, 0)

        else:
            raise ValueError(f"Invalid option {ci_method} for calculation of CI")

    except Exception as e:
        print(f"Error in Deming regression: {e}")
        slope = intercept = 0
        slope_ci = intercept_ci = (0, 0)

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
        dem_dict = deming_regression(x, y, lambda_, n_bootstrap, ci)
        ci_method = "Bootstrap"
        corr = stats.pearsonr(x, y)

    elif reg_method == "passing":
        dem_dict = passing_bablok_regression(x, y,
                                             confidence=ci, ci_method=ci_method, n_bootstrap=n_bootstrap, n_jobs=n_jobs)
        corr = stats.spearmanr(x, y)

    else:
        raise ValueError(f"Invalid option {ci_method} for calculation of CI")

    add_data = {"regression method": reg_method,
                "CI method": ci_method,
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


def sen_spe(ref, test, norm_range):
    """
    Calculate sensitivity, specificity and accuracy for according to normal ranges, based on reference and test
    measurements.
    Args:
        ref (1d array-like): Ground truth values
        test (1d array-like): Measured values (same shape as ref)
        norm_range (1d array-like): Bottom and top limits of the normal range, will be included in the range.
    Returns:
        pred_vals (dict): Dictionary containing predictive values (between 0 and 1)
    """
    # Convert inputs to numpy arrays
    ref = np.asarray(ref)
    test = np.asarray(test)
    norm_range = np.asarray(norm_range)

    # Verify all inputs are appropriate.
    if ref.ndim != 1 or test.ndim != 1:
        raise ValueError(f"ref and/or test input array is not 1-dimensional.")
    if ref.shape != test.shape:
        raise ValueError(f"The ref and test input arrays should be of the exact same shape.")
    if len(norm_range) != 2 or not all(isinstance(v, (int, float, np.number)) for v in norm_range):
        raise ValueError(f"The norm_range parameter must include two numbers (bottom and top limits of normal range).")

    # Convert to binary classification (1 for abnormal, 0 for normal)
    normal_bottom = np.min(norm_range)
    normal_top = np.max(norm_range)
    gt_labels = (ref < normal_bottom) | (ref > normal_top)
    test_labels = (test < normal_bottom) | (test > normal_top)

    # Calculate confusion matrix elements
    true_positive = np.sum((gt_labels == 1) & (test_labels == 1))
    false_negative = np.sum((gt_labels == 1) & (test_labels == 0))
    true_negative = np.sum((gt_labels == 0) & (test_labels == 0))
    false_positive = np.sum((gt_labels == 0) & (test_labels == 1))

    # Calculate sensitivity, specificity and accuracy
    sensitivity = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0.0
    specificity = true_negative / (true_negative + false_positive) if (true_negative + false_positive) > 0 else 0.0
    accuracy = (true_positive + true_negative) / (true_negative + false_positive + true_positive + false_negative)

    # Return results as a dictionary
    pred_vals = {'sensitivity': sensitivity, 'specificity': specificity, 'accuracy': accuracy}
    return pred_vals


def create_conf_mtrx(x, y, labels=[0, 1, 2, 3]):
    """
    Generates a confusion matrix comparing two 1D arrays of equal size.
    Args:
        x (array-like): 1D reference / ground truth labels (horizontal axis labels)
        y (array-like): 1D test / predicted labels (vertical axis labels)
        labels (list): possible grades, in order of appearance in table's axis
    Returns:
        pd.DataFrame: confusion matrix
    """
    if len(y) != len(x):
        raise ValueError("Input arrays must be of the same length.")
    mat_size = len(labels)

    # count appearances of each combination of grades
    unique_pairs, counts = np.unique(np.vstack((x, y)).T, axis=0, return_counts=True)

    # initialize and fill matrix
    conf_mtrx = np.zeros((mat_size, mat_size), dtype=int)
    for (hor, ver), count in zip(unique_pairs, counts):
        conf_mtrx[hor, ver] = count

    # convert to DatFrame
    df_conf = pd.DataFrame(conf_mtrx, index=labels, columns=labels)
    return df_conf


# def bootstrap_CI()



# in future move to dataframe utilities or cleaning utilities
def force_num_cols(df, num_cols, droprows=True, copy=True):
    """
    Make sure all values in specified columns are numeric. Turn to nan or remove rows with any non-numeric values.
    Args:
        df (dataframe): dataframe including the needed columns
        numeric_cols (list or str): columns with values that should be numeric
        droprows (bool): if True, drop rows with any non-numeric values in specified columns
        copy (bool): if True, create copy of original df instead of adjusting original

    Returns:
        pd.DataFrame: The modified DataFrame (if copy=True) or None (if copy=False)
    """
    # Ensure columns is a list, even if a single string is passed
    if isinstance(num_cols, str):
        num_cols = [num_cols]
    # Create a new DataFrame if copy=True
    if copy:  # until creation of wrapper class or dynamic method, this will not be used
        df = df.copy()

    # apply numeric conversion
    df.loc[:, num_cols] = df.loc[:, num_cols].apply(pd.to_numeric, errors='coerce')

    # remove rows with missing/non-numeric values
    if droprows:
        df.dropna(axis='index', subset=num_cols, inplace=True)

    if copy:
        return df

# in future, consider creating wrapper class to dataframe and adding the function to it instead
# pd.DataFrame.new_funct = force_num_cols
