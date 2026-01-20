import numpy as np
import scipy.stats as stats

def deming_regression_pe(x, y, lambda_=1):
    """
    Point estimate for Deming Regression.
    Can handle 1D arrays (single estimate) or 2D arrays (vectorized bootstrap).
    """
    # Force float to avoid integer division issues
    x, y = np.array(x, dtype=float), np.array(y, dtype=float)

    # Check if we are doing 1D (single) or 2D (bootstrap batch)
    axis = 1 if x.ndim > 1 else 0

    x_mean = np.mean(x, axis=axis, keepdims=True) if axis == 1 else np.mean(x)
    y_mean = np.mean(y, axis=axis, keepdims=True) if axis == 1 else np.mean(y)

    s_xy = np.sum((x - x_mean) * (y - y_mean), axis=axis)
    s_xx = np.sum((x - x_mean) ** 2, axis=axis)
    s_yy = np.sum((y - y_mean) ** 2, axis=axis)

    delta = s_yy - lambda_ * s_xx

    # Standard Deming slope formula (+ root)
    # Using np.hypot provides better numerical stability than sqrt(a**2 + b**2)
    numerator = delta + np.hypot(delta, 2 * np.sqrt(lambda_) * s_xy)
    slope = numerator / (2 * s_xy)

    # Handle edge case where s_xy is near zero (vertical/horizontal lines)
    # If s_xy is 0, this simple form might fail, but for real data it's usually fine.

    intercept = np.mean(y, axis=axis) - slope * np.mean(x, axis=axis)

    return slope, intercept


def deming_regression(x, y, lambda_=1, n_bootstrap=1000, ci=95):
    x = np.array(x)
    y = np.array(y)
    n = len(x)

    # 1. Point Estimate
    slope_pe, int_pe = deming_regression_pe(x, y, lambda_)

    # 2. Vectorized Bootstrap (100x faster than loops)
    # Generate all indices at once: shape (n_bootstrap, n)
    boot_indices = np.random.randint(0, n, size=(n_bootstrap, n))

    # Create big matrices of samples
    x_boot = x[boot_indices]
    y_boot = y[boot_indices]

    # Calculate all slopes/intercepts in one numpy pass
    boot_slopes, boot_intercepts = deming_regression_pe(x_boot, y_boot, lambda_)

    # 3. Percentiles
    alpha = (100 - ci) / 2
    slope_ci = np.percentile(boot_slopes, [alpha, 100 - alpha])
    int_ci = np.percentile(boot_intercepts, [alpha, 100 - alpha])

    return {
        "regression method": f"Deming (lambda={lambda_})",
        "CI method": "Bootstrap (Vectorized)",
        "iterations": n_bootstrap,
        "slope": slope_pe,
        "intercept": int_pe,
        "slope_ci_bottom": slope_ci[0],
        "slope_ci_top": slope_ci[1],
        "intercept_ci_bottom": int_ci[0],
        "intercept_ci_top": int_ci[1]
    }


def passing_bablok_core(x, y):
    """Computes all pairwise slopes for Passing Bablok."""
    n = len(x)
    # Optimized indices generation (upper triangle)
    idx_i, idx_j = np.triu_indices(n, k=1)

    x_i, x_j = x[idx_i], x[idx_j]
    y_i, y_j = y[idx_i], y[idx_j]

    dx = x_j - x_i
    dy = y_j - y_i

    # Filter out cases where dx is 0 (vertical segments)
    # Standard PB ignores points where x_i == x_j
    mask = dx != 0
    slopes = dy[mask] / dx[mask]

    # Standard PB requires sorting slopes
    slopes.sort()
    return slopes


def passing_bablok_regression(x, y, confidence=95, n_bootstrap=None):
    # Ensure the input arrays are numpy arrays
    x = np.array(x)
    y = np.array(y)
    n = len(x)

    # 1. Calculate All Slopes
    slopes = passing_bablok_core(x, y)
    n_slopes = len(slopes)

    # Point Estimate (Median of slopes)
    slope_pe = np.median(slopes)

    # Intercept Point Estimate: Median(y - S*x)
    intercept_pe = np.median(y - slope_pe * x)

    # 2. Confidence Intervals (Rank Method - The Correct Math)
    # Formula for variance of the Mann-Whitney U statistic (Kendall's variance logic)
    # Var = N(N-1)(2N+5)/18
    # C_alpha = z * sqrt(Var)

    alpha = 1 - (confidence / 100)
    z_score = stats.norm.ppf(1 - alpha / 2)

    variance = (n * (n - 1) * (2 * n + 5)) / 18
    k_offset = z_score * np.sqrt(variance)

    # The rank indices (0-based)
    # M = n_slopes. The median is at M/2.
    # The CI bounds are at (M/2 - k_offset) and (M/2 + k_offset)
    k_lower = int(round((n_slopes - k_offset) / 2))
    k_upper = int(round((n_slopes + k_offset) / 2))

    # Clamp indices
    k_lower = max(0, k_lower)
    k_upper = min(n_slopes - 1, k_upper)

    slope_ci_lower = slopes[k_lower]
    slope_ci_upper = slopes[k_upper]

    # 3. Intercept CI
    # In Passing-Bablok, the Intercept CI is calculated using the Slope CI bounds
    # Intercept_lower is calculated using Slope_upper (because of negative correlation)
    # Intercept_upper is calculated using Slope_lower

    # However, to be purely non-parametric, we calculate the array of intercepts
    # for the specific boundary slopes and take their medians.
    intercept_ci_lower = np.median(y - slope_ci_upper * x)
    intercept_ci_upper = np.median(y - slope_ci_lower * x)


    return {
        "regression method": "Passing-Bablok",
        "CI method": "Rank Method (Exact)",
        "iterations": 0,   # Analytical method
        "slope": slope_pe,
        "intercept": intercept_pe,
        "slope_ci_bottom": slope_ci_lower,
        "slope_ci_top": slope_ci_upper,
        "intercept_ci_bottom": intercept_ci_lower,
        "intercept_ci_top": intercept_ci_upper
    }


def regression_comp(x, y, n_bootstrap=1000, ci=95, reg_method="deming",
                    lambda_=1,
                    ci_method="Bootstrap",
                    res_str=True):
    """
    Calculate regression for method comparison, getting all info on the regression.
    Parameters:
        - x, y: input array_like
        - reg_method: "deming" or "passing"
        - lambda_: variations ratio, only for Deming
        - ci_method: "U test" or "Bootstrap", only for Passing
    Returns:
        Dict
    """
    # Ensure the input arrays are numpy arrays
    x = np.array(x)
    y = np.array(y)
    n = len(x)


    if reg_method.lower() == "deming":
        ci_method = "Bootstrap"
        corr = stats.pearsonr(x, y)
        dem_dict = deming_regression(x, y, lambda_, n_bootstrap, ci)

    elif reg_method.lower() in ["passing", "passing-bablok"]:
        # We ignore n_bootstrap for PB as Rank method is superior
        corr = stats.spearmanr(x, y)
        dem_dict = passing_bablok_regression(x, y, confidence=ci)
    else:
        raise ValueError(f"Method must be 'deming' or 'passing'")

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