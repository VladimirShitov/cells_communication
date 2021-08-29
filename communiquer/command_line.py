from pathlib import Path
import subprocess
from typing import NoReturn, Optional, Union
import warnings

from communiquer.utils import doc_string
from communiquer.docstrings import cellphone_db_parameters
# from communiquer.types import CPDB_METHOD, CPDB_OUTPUT_FORMAT, GENES_ENCODING


def python_args_to_command_line_args(command: Optional[list] = None, **kwargs) -> list:
    """Transform arguments of a python function to a command line args.
    Underscores are replaced with dashes.

    Parameters
    ----------
    command: Optional[list]
        List with command arguments to start with
    kwargs:
        Any key word arguments that you want to transform into command line arguments

    Returns
    -------
    List with command line arguments

    Examples
    --------
    >>> python_args_to_command_line_args(you_are="breathtaking")
    ['--you-are=breathtaking']

    >>> python_args_to_command_line_args(
            command=["cellphonedb", "method", "analysis"], project_name="pbmc3k")
    ['cellphonedb', 'method', 'analysis', '--project-name=pbmc3k']
    """
    args = command or []

    for key, value in kwargs.items():
        if value is not None:
            command_line_arg = key.replace("_", "-")
            args.append(f"--{command_line_arg}={value}")

    return args


@doc_string(cellphone_db_parameters)
def cellphone_db_methods_command(
    meta_file_path: Union[Path, str],
    counts_file_path: Union[Path, str],
    method: str = "statistical_analysis",  # CPDB_METHOD
    counts_data: str = "ensembl",  # GENES_ENCODING
    project_name: str = "",
    threshold: float = 0.1,
    result_precision: int = 3,
    output_path: Union[Path, str] = "",
    output_format: str = "txt",  # CPDB_OUTPUT_FORMAT
    means_result_name: str = "means",
    significant_means_result_name: str = "significant_means",
    deconvoluted_result_name: str = "deconvoluted",
    verbose: bool = True,
    database: Optional[str] = None,
    subsampling: bool = False,
    subsampling_log: bool = False,
    subsampling_num_pc: int = 100,
    subsampling_num_cells: Optional[int] = None,
    debug_seed: int = -1,
    pvalue: float = 0.05,
    pvalues_result_name: str = "pvalues",
    iterations: int = 1000,
    threads: int = -1,
) -> str:
    """Generate a bash command to run CellPhoneDB

    Parameters
    ----------
    method: str = "statistical_analysis"
        Either "analysis" or "statistical_analysis". The first option doesn't
        generate p-values, which is much faster, but less informative.
    """

    command = ["cellphonedb", "method", method]

    command.extend(
        python_args_to_command_line_args(
            counts_data=counts_data,
            project_name=project_name,
            threshold=threshold,
            result_precision=result_precision,
            output_path=output_path,
            output_format=output_format,
            means_result_name=means_result_name,
            significant_means_result_name=significant_means_result_name,
            deconvoluted_result_name=deconvoluted_result_name,
            database=database,
        )
    )

    if subsampling:
        command.append("--subsampling")
        command.extend(
            python_args_to_command_line_args(
                subsampling_log=subsampling_log,
                subsampling_num_pc=subsampling_num_pc,
                subsampling_num_cells=subsampling_num_cells,
            )
        )

    if method == "statistical_analysis":
        command.extend(
            python_args_to_command_line_args(
                debug_seed=debug_seed,
                pvalue=pvalue,
                pvalues_result_name=pvalues_result_name,
                iterations=iterations,
                threads=threads,
            )
        )

    command.append("--verbose" if verbose else "--quiet")

    command.extend([str(meta_file_path), str(counts_file_path)])

    return " ".join(command)


def cellphonedb_heatmap_plot_command(
    meta_file_path: Union[Path, str],
    pvalues_path: Union[Path, str] = "./out/pvalues.txt",
    output_path: Union[Path, str] = "./out",
    count_name: str = "heatmap_count.pdf",
    log_name: str = "heatmap_log_count.pdf",
    count_network_name: str = "count_network.txt",
    interaction_count_name: str = "interactions_count.txt",
    pvalue: float = 0.05,
    verbose: bool = True
) -> str:
    command = ["cellphonedb", "plot", "heatmap_plot"]

    command.extend(
        python_args_to_command_line_args(
            pvalues_path=pvalues_path,
            output_path=output_path,
            count_name=count_name,
            log_name=log_name,
            count_network_name=count_network_name,
            interaction_count_name=interaction_count_name,
            pvalue=pvalue
        )
    )

    command.append("--verbose" if verbose else "--quiet")
    command.append(str(meta_file_path))

    return " ".join(command)


def cellphonedb_dot_plot_command(
    means_path: Union[Path, str],
    pvalues_path: Union[Path, str],
    output_path: Union[Path, str],
    output_name: str = "plot.pdf",
    rows: Optional[Union[Path, str]] = None,
    columns: Optional[Union[Path, str]] = None,
    verbose: bool = True
) -> str:
    command = ["cellphonedb", "plot", "dot_plot"]

    command.extend(
        python_args_to_command_line_args(
            means_path=means_path,
            pvalues_path=pvalues_path,
            output_path=output_path,
            output_name=output_name,
        )
    )

    if rows is not None:
        command.extend(python_args_to_command_line_args(rows=rows))

    if columns is not None:
        command.extend(python_args_to_command_line_args(columns=columns))

    command.append("--verbose" if verbose else "--quiet")

    return " ".join(command)


def run_command_line_subprocess(command: str, verbose: bool) -> NoReturn:
    """Run a command line task

    Parameters
    ----------
    command: str
        Bash command to be called
    verbose: bool
        If True, print stdout
    """
    command_line_subprocess = subprocess.run(command, shell=True, capture_output=True)

    if verbose:
        print(command_line_subprocess.stdout.decode())

    if command_line_subprocess.stderr:
        warnings.warn(command_line_subprocess.stderr.decode())
