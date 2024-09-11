"""Provide utility functions used across the process steps and across
themes."""

from pprint import pprint

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
    else:
        for output_tuple, expect_tuple in zip(
            df_output.iterrows(), df_expect.iterrows()
        ):
            output_row = output_tuple[1]
            expect_row = expect_tuple[1]
            for col in columns:
                if output_row[col] != expect_row[col]:
                    pprint(output_row)
                    pprint(expect_row)
                    break
    row_accuracy = float(rows_intersect.shape[0]) / rows_union.shape[0]
    return row_accuracy
