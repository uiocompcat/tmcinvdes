"""Provide utility functions used across the process steps and across
themes."""


import pandas as pd


def compare_dataframes(
    df_output: pd.DataFrame, df_expect: pd.DataFrame, columns: list[str]
) -> float:
    """Compare the reproduced and expected dataframe and calculate accuracy in
    terms of identical rows.

    Args:
        df_output (pd.DataFrame): the reproduced dataframe that would be the output written to file
        if not testing.
        df_expect (pd.DataFrame): the target output read back into a dataframe from file.
        columns (list[str]): the shared dataframe columns to use for inner and outer joins.

    Returns:
        float: accuracy as a number between 0.0 and 1.0 representing the proportion of overlapping
        identical rows between the dataframes being compared.
    """
    df_output = df_output.round(6)
    df_expect = df_expect.round(6)
    rows_intersect = pd.merge(
        df_output,
        df_expect,
        how="inner",
        on=columns,
    )
    rows_union = pd.merge(
        df_output,
        df_expect,
        how="outer",
        on=columns,
    )
    if rows_intersect.equals(rows_union):
        return 1.0
    row_accuracy = float(rows_intersect.shape[0]) / rows_union.shape[0]
    return row_accuracy
