""" build_*_plot functions, each creating a single plot (which can be multi-layered """
from selectors import *
from composition import PlotLayer, CompositePlot


def build_method_comparison_plot(
    ref_method,
    test_method,
    variable=None,
    style=None
):
    layers = [
        PlotLayer(
            plot_func=scatter_basic,
            data_resolver=select_method_comparison_xy(
                ref_method, test_method, variable
            ),
            style=style or DEFAULT_SCATTER_STYLE
        ),
        PlotLayer(
            plot_func=overlay_regression_line,
            data_resolver=select_method_comparison_xy(
                ref_method, test_method, variable
            ),
            kwargs={
                "reg": select_regression(
                    ref_method, test_method, variable
                )
            }
        )
    ]

    return CompositePlot(
        layers=layers,
        title=f"{test_method} vs {ref_method}",
        xlabel=ref_method,
        ylabel=test_method
    )
