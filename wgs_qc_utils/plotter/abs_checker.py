from pandas.api import types as pd_types
from enum import Enum
import pandas as pd

class CheckerTypes(Enum):
    """
    essentially abstracts pandas/numpy type checking from types
    """
    STRING = pd_types.is_object_dtype
    INT = pd_types.is_integer_dtype
    FLOAT = pd_types.is_float_dtype
    BOOLEAN = pd_types.is_bool_dtype


def check_input_is_valid(data, type_checks):
    """
    check that input data/columns follow a specified set of types
    :param data: input data
    :param data_cols: names in data to check
    :param type_checks: list of CheckTypes to check with
    :return: assertion
    """
    if isinstance(data, pd.DataFrame):
        data = [data[col] for col in data]
     
    for type_check, col in zip(type_checks, data):
        assert type_check(col.dtype)


