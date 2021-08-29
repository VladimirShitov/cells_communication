from typing import NoReturn, Optional, Union
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2.robjects.packages as rpackages
from rpy2 import robjects

from communiquer.drawing import draw_chord_diagram, draw_dotplot
from communiquer.command_line import cellphone_db_methods_command, \
    run_command_line_subprocess, cellphonedb_heatmap_plot_command, cellphonedb_dot_plot_command
from communiquer.utils import matrix_to_long_data_frame


class CellCommunicationsLauncher:
    def read_meta(self):
        self.meta_df = pd.read_csv(self.meta_file_path, sep="\t")
        self.cell_types = self.meta_df[self.cell_labels_column].unique()

    def __init__(
            self,
            meta_file_path: Union[Path, str],
            counts_file_path: Union[Path, str],
            cell_labels_column: str = "cell_type"
    ):
        self.meta_file_path: Union[Path, str] = meta_file_path
        self.counts_file_path: Union[Path, str] = counts_file_path
        self.cell_labels_column: Union[Path, str] = cell_labels_column

        self.meta_df: pd.DataFrame = None
        self.cell_types: np.array = None

        # Cellchat-like interface fields
        self.pvalues_dfs = None
        self.counts_df = None

        self.read_meta()

    def visualise_interactions(self, ax=None, plot_size=10):
        """Visualize amount of significant interactions between cell types

        Parameters
        ----------
        ax: matplotlib Axis (optional)
            Axis on which you want to draw a plot

        plot_size: int = 7
            Size of the plot sides in inches

        Returns
        -------
        Matplotlib axis with a plot
        """

        if self.counts_df is None:
            raise ValueError("Can't visualise connections. Please, fill in `counts_df` first")

        return draw_chord_diagram(matrix=self.counts_df, ax=ax, plot_size=plot_size)

    def dotplot_counts(self) -> plt.Axes:
        long_counts_df = matrix_to_long_data_frame(self.counts_df)
        return draw_dotplot(long_counts_df)


class CellPhoneDBLauncher(CellCommunicationsLauncher):
    def _df_to_cellchat_format(self, df_to_convert) -> dict:
        """Convert a data frame with CellPhoneDB output to dict of dataframes for each interaction"""
        interactions_dataframes = {}

        for _, row in df_to_convert.iterrows():
            df = pd.DataFrame(index=self.cell_types, columns=self.cell_types)

            for cell_type_i in self.cell_types:
                for cell_type_j in self.cell_types:
                    interaction_name = f"{cell_type_i}|{cell_type_j}"
                    df.loc[cell_type_i, cell_type_j] = row[interaction_name]

            interactions_dataframes[row.interacting_pair] = df

        return interactions_dataframes

    def __init__(
            self,
            meta_file_path: Union[Path, str],
            counts_file_path: Union[Path, str],
            cell_labels_column: str = "cell_type",
            method="statistical_analysis",  # : CPDB_METHOD
            counts_data: str = "ensembl",  # : GENES_ENCODING
            project_name: str = "",
            threshold: float = 0.1,
            result_precision: int = 3,
            output_path: Union[Path, str] = "out/",
            output_format: str = "txt",  # : CPDB_OUTPUT_FORMAT
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
            threads: int = -1
    ):
        super().__init__(
            meta_file_path=meta_file_path,
            counts_file_path=counts_file_path,
            cell_labels_column=cell_labels_column
        )

        self.method = method
        self.counts_data = counts_data
        self.project_name = project_name
        self.threshold = threshold
        self.result_precision = result_precision
        self.output_path = output_path
        self.output_format = output_format
        self.means_result_name = means_result_name
        self.significant_means_result_name = significant_means_result_name
        self.deconvoluted_result_name = deconvoluted_result_name
        self.verbose = verbose
        self.database = database
        self.subsampling = subsampling
        self.subsampling_log = subsampling_log
        self.subsampling_num_pc = subsampling_num_pc
        self.subsampling_num_cells = subsampling_num_cells
        self.debug_seed = debug_seed
        self.pvalue = pvalue
        self.pvalues_result_name = pvalues_result_name
        self.iterations = iterations
        self.threads = threads

        self.output_dir = Path(self.output_path) / self.project_name

        self.command = None
        self.means_df = None
        self.significant_means_df = None
        self.pvalues_df = None
        self.deconvoluted_df = None

        # Fields for cellchat output format
        self.means_dfs = None
        self.pvalues_dfs = None

    def run(self, verbose=True) -> NoReturn:
        self.command = cellphone_db_methods_command(
            meta_file_path=self.meta_file_path,
            counts_file_path=self.counts_file_path,
            method=self.method,
            counts_data=self.counts_data,
            project_name=self.project_name,
            threshold=self.threshold,
            result_precision=self.result_precision,
            output_path=self.output_path,
            output_format=self.output_format,
            means_result_name=self.means_result_name,
            significant_means_result_name=self.significant_means_result_name,
            deconvoluted_result_name=self.deconvoluted_result_name,
            verbose=self.verbose,
            database=self.database,
            subsampling=self.subsampling,
            subsampling_log=self.subsampling_log,
            subsampling_num_pc=self.subsampling_num_pc,
            subsampling_num_cells=self.subsampling_num_cells,
            debug_seed=self.debug_seed,
            pvalue=self.pvalue,
            pvalues_result_name=self.pvalues_result_name,
            iterations=self.iterations,
            threads=self.threads
        )

        if verbose:
            print("Running command", self.command)

        run_command_line_subprocess(command=self.command, verbose=verbose)

    def output_file_path(self, output_file_name) -> Path:
        return self.output_dir / f"{output_file_name}.{self.output_format}"

    def read_output(self, convert_to_cellchat_format=True) -> NoReturn:
        """Create dataframes from cellphoneDB output files"""
        self.means_df = pd.read_csv(self.output_file_path(self.means_result_name), sep="\t")

        self.significant_means_df = pd.read_csv(
            self.output_file_path(self.significant_means_result_name), sep="\t")

        self.pvalues_df = pd.read_csv(self.output_file_path(self.pvalues_result_name), sep="\t")

        self.deconvoluted_df = pd.read_csv(
            self.output_file_path(self.deconvoluted_result_name), sep="\t")

        if convert_to_cellchat_format:
            self.convert_output_to_cellchat_format()

    def convert_output_to_cellchat_format(self):
        self.means_dfs = self._df_to_cellchat_format(self.means_df)
        self.pvalues_dfs = self._df_to_cellchat_format(self.pvalues_df)

    def inbuilt_heatmap_plot(
            self,
            output_path: Union[Path, str] = "./out",
            counts_file_path: str = "heatmap_count.pdf",
            log_name: str = "heatmap_log_count.pdf",
            count_network_name: str = "count_network.txt",
            interaction_counts_file_path: str = "interactions_count.txt",
            pvalue: float = 0.05,
            verbose: bool = True
    ):
        """Plot a heatmap with counts by interaction pair

        Plot is drawn by running cellphoneDB program.
        This plot type requires pheatmap R package installed and working.

        Parameters
        ----------
        output_path: Union[Path, str] = "./out"
            Output folder for the plot
        counts_file_path: str = "heatmap_count.pdf"
            Filename of the output plot. You can choose format from pdf, png and jpeg
        log_name: str = "heatmap_log_count.pdf"
            Filename of the output plot using log-count of interactions.
            You can choose format from pdf, png and jpeg
        count_network_name: str = "count_network.txt"
            Filename of the output network file
        interaction_counts_file_path: str = "interactions_count.txt"
            Filename of the output interactions-count file
        pvalue: float = 0.05
            pvalue threshold to consider when plotting
        verbose: bool = True
            Print or hide logs
        """

        command = cellphonedb_heatmap_plot_command(
            meta_file_path=self.meta_file_path,
            pvalues_path=self.output_file_path(self.pvalues_result_name),
            output_path=output_path,
            counts_file_path=counts_file_path,
            log_name=log_name,
            count_network_name=count_network_name,
            interaction_counts_file_path=interaction_counts_file_path,
            pvalue=pvalue,
            verbose=verbose
        )

        if verbose:
            print("Running command", command)

        run_command_line_subprocess(command=command, verbose=verbose)

    def inbuilt_dot_plot(
            self,
            output_path: Union[Path, str],
            output_name: str = "plot.pdf",
            rows: Optional[Union[Path, str]] = None,
            columns: Optional[Union[Path, str]] = None,
            verbose: bool = True
    ):
        """Make a dot plot with counts by interaction pair

        Plot is drawn by running cellphoneDB program.
        This plot type requires ggplot2 R package installed and working.

        Parameters
        ----------
        output_path: Union[Path, str] = "./out"
            Output folder for the plot
        rows: Optional[Union[Path, str]] = None
            File with a list of rows to plot, one per line.
            By default, all rows are used
        columns: Optional[Union[Path, str]] = None
            File with a list of columns to plot, one per line.
            By default, all columns are used
        verbose: bool = True
            Print or hide logs
        """

        command = cellphonedb_dot_plot_command(
            means_path=self.output_file_path(self.means_result_name),
            pvalues_path=self.output_file_path(self.pvalues_result_name),
            output_path=output_path,
            output_name=output_name,
            rows=rows,
            columns=columns,
            verbose=verbose
        )

        if verbose:
            print("Running command", command)

        run_command_line_subprocess(command=command, verbose=verbose)

    def count_significant_interactions(self, alpha=None) -> NoReturn:
        """Count significant interactions between cell types

        Parameters
        ----------
        alpha: Optional[float] = None
            A significal level between 0 and 1. Interactions with p-value less than `alpha`
            are considered significant. If `alpha` is not set, `self.pvalue` is used,
            which had been set during the initialization

        Sets
        ----
        counts_df

        """
        p_value_threshold = alpha or self.pvalue

        self.counts_df = pd.DataFrame(
            np.zeros(shape=(len(self.cell_types), len(self.cell_types))),
            columns=self.cell_types,
            index=self.cell_types,
            dtype=np.int16
        )

        for interaction_df in self.pvalues_dfs.values():
            self.counts_df += (interaction_df < p_value_threshold)


class CellChatLauncher(CellCommunicationsLauncher):
    PIPELINE_CODE = cellchat_pipeline_code = """
    library(CellChat)

    cellchat_pipeline = function(
      counts.file.path,
      meta.file.path,
      cell.labels.column = "cell_type",
      use.human.db = T,
      db.subset = "",
      project.to.ppi = F,
      min.cells = 10
    ) {
      data.input <- read.csv(counts.file.path, sep="\t", row.names = 1, check.names = F)
      meta <- read.csv(meta.file.path, sep="\t", row.names = 1)

      cellchat <- createCellChat(object = as.matrix(data.input),
                                 meta = meta, 
                                 group.by = cell.labels.column)

      cellchat <- addMeta(cellchat, meta = meta)
      cellchat <- setIdent(cellchat, ident.use = cell.labels.column)

      groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

      if (use.human.db) {
        CellChatDB <- CellChatDB.human  
      } else {
        CellChatDB <- CellChatDB.mouse  
      }


      # We can choose "Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"
      if(db.subset != "") {
        CellChatDB <- subsetDB(CellChatDB, search = db.subset)   
      }

      cellchat@DB <- CellChatDB

      cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)

      if(project.to.ppi) {
        # project gene expression data onto PPI network (optional)
        cellchat <- projectData(cellchat, PPI.human)  
      }


      # Compute the communication probability and infer cellular communication network
      cellchat <- computeCommunProb(cellchat)
      # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
      cellchat <- filterCommunication(cellchat, min.cells = min.cells)

      df.net <- subsetCommunication(cellchat) 

      # Infer the cell-cell communication at a signaling pathway level
      cellchat <- computeCommunProbPathway(cellchat)

      # Calculate the aggregated cell-cell communication network
      cellchat <- aggregateNet(cellchat)

      return(cellchat)
    }
    """

    def __init__(
        self,
        counts_file_path: Union[Path, str],
        meta_file_path: Union[Path, str],
        cell_labels_column: str = "cell_type",
        use_human_db: bool = True,
        db_subset: str = "",  # CELLCHAT_DB_TYPE
        project_to_ppi: bool = False,
        min_cells: int = 10,
        install_dependencies: bool = False
    ):
        """
        Initializes launcher for CellChat

        Parameters
        ----------
        counts_file_path: Union[Path, str]
            Path to a file with normalized counts. Columns contain cells, and rows contain genes
        meta_file_path: Union[Path, str]
            Path to a meta file. It should have at least 2 columns. The first one must
            contain cells' IDs. Another column should contain cells' labels
        cell_labels_column: str = "cell_type"
            Name of the column in the `meta_file_path`, which contains cells' labels
        use_human_db: bool = True
            If True, use human interactions database. Otherwise, use a database for mouse
        db_subset: str = ""  # CELLCHAT_DB_TYPE
            If not empty, use subset of the interactions database. Can be one of
            "Secreted Signaling", "ECM-Receptor" or "Cell-Cell Contact"
        project_to_ppi: bool = False
            Project gene expression data onto PPI network
        min_cells: int = 10
            Filter out the cell-cell communication if there are
            less than this number of cells in certain cell groups
        install_dependencies: bool = False
            If True, install required R dependencies for the CellChat
        """
        super().__init__(
            meta_file_path=meta_file_path,
            counts_file_path=counts_file_path,
            cell_labels_column=cell_labels_column
        )

        if install_dependencies:
            self._install_dependencies()

        self.pipeline = robjects.r(self.PIPELINE_CODE)

        self.counts_file_path = counts_file_path
        self.meta_file_path = meta_file_path
        self.cell_labels_column = cell_labels_column
        self.use_human_db = use_human_db
        self.db_subset = db_subset
        self.project_to_ppi = project_to_ppi
        self.min_cells = min_cells

        self.cellchat = None

        self.probabilities_dfs = None
        self.weights_df = None

    @staticmethod
    def _install_dependencies() -> NoReturn:
        utils = rpackages.importr("utils")
        utils.chooseCRANmirror(ind=1)  # select the first mirror in the list

        packages = ("dplyr", "ggplot2", "patchwork", "BiocManager", "Cairo", "NMF", "devtools")

        for package in packages:
            if not rpackages.isinstalled(package):
                utils.install_packages(package)

        robjects.r('BiocManager::install("Biobase")')

        # Uncomment if there is an error with tar
        # robjects.r('Sys.setenv(TAR = "/bin/tar")')

        github_packages = ("jokergoo/circlize", "jokergoo/ComplexHeatmap", "sqjin/CellChat")

        for package in github_packages:
            robjects.r(f'devtools::install_github("{package}")')

    def _read_slot(self, slot_name: str, has_third_dimension: bool):
        net = self.cellchat.slots["net"]
        slot_names = list(net.names)

        slot_index = slot_names.index(slot_name)
        slot_matrix = net[slot_index]

        if has_third_dimension:
            index, columns, interactions = map(tuple, slot_matrix.names)
            new_shape = (len(index), len(columns), len(interactions))

            array = np.array(list(slot_matrix)).reshape(new_shape)

            data_frames = {}

            for i, interaction_name in enumerate(interactions):
                df = pd.DataFrame(array[:, :, i], index=index, columns=columns)
                data_frames[interaction_name] = df

            return data_frames

        else:
            index, columns = map(tuple, slot_matrix.names)
            new_shape = (len(index), len(columns))

            array = np.array(list(slot_matrix)).reshape(new_shape)

            df = pd.DataFrame(array, index=index, columns=columns)

            return df

    def run(self):
        self.cellchat = self.pipeline(
            counts_file_path=str(self.counts_file_path),
            meta_file_path=str(self.meta_file_path),
            cell_labels_column=self.cell_labels_column,
            use_human_db=self.use_human_db,
            db_subset=self.db_subset,
            project_to_ppi=self.project_to_ppi,
            min_cells=self.min_cells
        )

    def read_output(self):
        self.probabilities_dfs = self._read_slot(slot_name="prob", has_third_dimension=True)
        self.pvalues_dfs = self._read_slot(slot_name="pval", has_third_dimension=True)
        self.weights_df = self._read_slot(slot_name="weight", has_third_dimension=False)
        self.counts_df = self._read_slot(slot_name="count", has_third_dimension=False)
