import pandas as pd
import numpy as np
from scipy.stats import skewnorm, norm
import json, os
from typing import Any, NewType

InchiType = NewType('InchiType', str)

def write_jsonl(obj: Any, filename: str) -> None:
    with open(filename, 'a') as fh:
        fh.write(json.dumps(obj) + '\n')

def read_jsonl(filename: str) -> Any:
    if not os.path.exists(filename):
        return []
    with open(filename) as fh:
        data = []
        for line in fh:
            try:
                data.append(json.loads(line))
            except json.JSONDecodeError as error:
                pass  # burnt line
        return data

def get_skewnorm_params(series: pd.Series, remove_zeros=False) -> dict:
    mask = ~series.isna()
    if remove_zeros:
        mask = mask & (series != 0.)
    series: pd.Series = series.loc[mask].copy()
    # winsorise the infinites
    finite_min = series.loc[series.abs() != np.inf].min()
    finite_max = series.loc[series.abs() != np.inf].max()
    series.loc[series == -np.inf] = finite_min
    series.loc[series == +np.inf] = finite_max
    raw: np.array = series.values
    filtered = raw[np.abs(raw - np.mean(raw)) <= 3 * np.std(raw)]
    alpha, loc, scale = skewnorm.fit(filtered)
    return {'name': series.name, 'count': len(series), 'alpha': alpha, 'loc': loc, 'scale': scale}

def ultranormalize(value: float,
                   skew_loc: float=0,
                   skew_scale: float=1,
                   skew_shape:float=0,
                   nor_bound: float=2,
                   flip_sign: int=+1,
                   nan_replacement=-np.inf) -> float:
    """
    Normalise a skewed value to a Zscore, smoothed clipped at the bounds.
    Assuming it is a skew normal, unskew it by quantile transformation (by skewnormal CDF & probit mapping)
    and then smooth clip it at the bounds (`tanh`).

    :param value: the value to normalise
    :param skew_loc: skew omega
    :param skew_scale: skew xi
    :param skew_shape: skew alpha
    :param nor_bounds: normal bounds, i.e. µ = 0, σ = 1, not in skew distro
    :param flip_sign: +1 or -1
    :param nan_replacement:
    :return:
    """
    if np.isnan(value):
        return nan_replacement
    elif np.isposinf(value):
        return flip_sign * nor_bound
    elif np.isneginf(value):
        if np.isposinf(value):
            return - flip_sign * nor_bound
    else:
        # unskew by quantile transformation (by skewnormal CDF & probit mapping)
        # `skewnorm.cdf` does not have argument name for shape / alpha
        quantile = skewnorm.cdf((value - skew_loc) / skew_scale, skew_shape, loc=0, scale=1)
        qtrans = flip_sign * norm.ppf(quantile)
        # smooth clip it at the bounds
        asymptote = np.abs(1 / nor_bound)
        return np.tanh(qtrans * asymptote) / asymptote


def autopass_fun(value, bound):
    if value == 0.:
        return -bound
    else:
        return +bound