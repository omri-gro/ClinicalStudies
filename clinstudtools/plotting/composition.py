from dataclasses import dataclass, field
from typing import Callable, Optional, Dict, Sequence, Any

@dataclass
class PlotLayer:
    """ single visual element of a plot (scatter, line, bars, etc.) """
    plot_func: Callable
    data_resolver: Callable  # takes results → data
    style: Optional[Dict[str, Any]] = None
    kwargs: Dict[str, Any] = field(default_factory=dict)

    def draw(self, results, fig=None, ax=None):
        data = self.data_resolver(results)
        return self.plot_func(
            data=data,
            fig=fig,
            ax=ax,
            style=self.style,
            **self.kwargs
        )


class CompositePlot:
    def __init__(self, layers: Sequence[PlotLayer],
                 title=None, xlabel=None, ylabel=None):
        self.layers = layers
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel

    def draw(self, results, fig=None, ax=None):
        for layer in self.layers:
            fig, ax = layer.draw(results, fig=fig, ax=ax)

        if ax is not None:
            if self.title:
                ax.set_title(self.title)
            if self.xlabel:
                ax.set_xlabel(self.xlabel)
            if self.ylabel:
                ax.set_ylabel(self.ylabel)

        return fig, ax
