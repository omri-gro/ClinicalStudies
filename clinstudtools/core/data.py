import pandas as pd
from dataclasses import dataclass, field
from typing import Any, Optional, Mapping, Tuple

from ..table_integrity import *

RowFilters = Mapping[str, Any]

@dataclass
class ComparisonData:
    df: pd.DataFrame
    metadata: Optional[Any] = None
    id_cols: Tuple[str, ...] = ("SampleID", "Site")
    history: Tuple[str, ...] = field(default_factory=tuple)




