from pathlib import Path
from typing import NoReturn, Union

import pandas as pd
import scanpy as sc


def create_cellphone_db_meta_file(
    adata: sc.AnnData, meta_file_path: Union[Path, str], clusters_key: str = "leiden"
) -> NoReturn:
    """Convert cells clusters information to the format accepted by CellPhoneDB

    Parameters
    ----------
    adata: sc.AnnData
        AnnData object with the information about single cell experiment. Must contain information
        about cells' clusters
    meta_file_path: Union[Path, str]
        Path, where to save the file
    clusters_key: str = "leiden"
        Key in `adata.obs`, which contains information about the group of the cells (e. g. clusters)
    """
    adata.obs[clusters_key].reset_index().to_csv(
        meta_file_path, sep="\t", header=("Cell", "cell_type"), index=False
    )


def doc_string(docstring: str = "default doc"):
    """Set a documentation string for `function`. Helpful if several functions have the same doc"""

    def wrapper(function):
        function.__doc__ += docstring

        return function

    return wrapper


def cellcall_cell_id(
    adata: sc.AnnData, cell_name: str, clusters_key: str = "leiden"
) -> str:
    """Transform cell ID from adata to cellcall format

    Description of the format: https://github.com/ShellyCoder/cellcall#211-load-data
    Briefly, cells IDs should not contain any separators (comma, dash, etc.) except for underline.
    Each cell has to be named as "<cell_id>_<cell_type>". Example of the cell ID produced by this function:
    "AAACCGTGCTTCCG_1_CD14 Monocytes"

    Don't forget to set names.field = 3


    Parameters
    ----------
    adata: scanpy.AnnData
        AnnData object with cells counts and clusters
    cell_name: str
        ID of the cell in the `adata`
    clusters_key: str = "leiden"
        Key containing cell types (or clusters) in `adata.obs`

    Returns
    -------
    Cell ID in the cellcall format
    """

    cell_type = adata.obs[clusters_key][cell_name]
    clear_cell_name = cell_name.replace("-", "_")

    return f"{clear_cell_name}_{cell_type}"


def counts_to_cellcall_format(adata, clusters_key="leiden"):
    cells_to_take = []  # Raw data might not contain some cells, take indexes of cells in `adata`
    cells_ids = []

    for cell_name in adata.obs_names:
        new_cell_id = cellcall_cell_id(adata, cell_name, clusters_key)
        cells_ids.append(new_cell_id)

        cell_idx = adata.obs_names.get_loc(cell_name)
        cells_to_take.append(cell_idx)

    df = pd.DataFrame.sparse.from_spmatrix(
        adata.raw.X[cells_to_take, :], index=cells_ids, columns=adata.raw.var_names
    )

    return df.T


def matrix_to_long_data_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Convert matrix, where index and columns are the same to the long format

    Parameters
    ----------
    df: pandas.DataFrame
        matrix, where index and columns are the same (e.g. cell types)

    Returns
    -------
    A data frame with columns "index" and "variable" containing all pairs of index and columns
    from `df`. Values from `df` are saved to a column "count"

    Examples
    --------
    >>> df = pd.DataFrame({"a": [1, 2], "b": [3, 4]}, index=["a", "b"])
    >>> print(df.to_string())
       a  b
    a  1  3
    b  2  4
    >>> long_df = matrix_to_long_data_frame(df)
    >>> print(long_df.to_string())
      index variable  count
    0     a        a      1
    1     b        a      2
    2     a        b      3
    3     b        b      4
    """
    long_df: pd.DataFrame = df.melt(ignore_index=False, value_name="count").reset_index()
    long_df["count"] = long_df["count"].astype(int)

    return long_df
