cellphone_db_parameters = """meta_file_path: Union[Path, str]
    Path to a meta file. It has 2 columns: "Cell" and "cell_type" with cells' IDs
    and cells' clusters (or some other types) respectively
counts_file_path: Union[Path, str]
    Path to a file with normalized counts. Columns contain cells, and rows contain genes
counts-data: one of "ensembl", "gene_name", "hgnc_symbol"
    Type of genes ID that is used in input parameters
project-name: str
    Name of the project. It creates a subfolder in output folder
threshold: float = 0.1
    % of cells expressing a gene
result-precision: int = 3
    Number of decimal digits in results
output-path: Union[Path, str] = "out/"
    Directory where the results will be allocated (the directory must exist)
output-format: one of "txt", "csv", "tsv", "tab"
    Format of the output files. "txt" by default
means-result-name: str = "means"
    Means result namefile
significant-means-result-name: str = "significant_means"
    Significant result namefile
deconvoluted-result-name: str = "deconvoluted"
    Deconvoluted result namefile
verbose: bool = True
    Print cellphonedb logs
database: Optional[str]
subsampling: bool = False
    Enable subsampling
subsampling-log: bool = False
    Enable subsampling log1p for non transformed data inputs. !mandatory!
subsampling-num-pc: int = 100
    Subsampling NumPC argument
subsampling-num-cells: Optional[int]
    Number of cells to subsample (defaults to a 1/3 of cells)
debug-seed: int = -1
    Debug random seed 0 for disable it. >=0 to set it. Ignored 
    if `method` is not "statistical_analysis"
pvalue: float = 0.05                  
    Pvalue threshold. Ignored if `method` is not "statistical_analysis"
pvalues-result-name: str = "pvalues"
    Pvalues result namefile. Ignored if `method` is not "statistical_analysis"
iterations: int = 1000
    Number of pvalues analysis iterations. Ignored if `method` is not
    "statistical_analysis"
threads: int = 4
    Max of threads to process the data. Ignored if `method` is not
    "statistical_analysis"
"""

