 #!/usr/bin/env python3
import argparse
from scipy.spatial.distance import _METRICS_NAMES
from scipy.cluster.hierarchy import _LINKAGE_METHODS


class Config:

    version = "1.1.0"
    default_rank_name = "default"
    output_plots = ["overview", "samples", "heatmap", "correlation"]

    def __new__(self, argv=None):

        parser = argparse.ArgumentParser(description='grimer')

        required_group = parser.add_argument_group('required arguments')
        required_group.add_argument('-i', '--input-file', required=True, type=str, help="Main input table with counts (Observation table, Count table, Contingency Tables, ...) or .biom file. By default rows contain observations and columns contain samples (use --transpose if your file is reversed). First column and first row are used as headers.")

        main_group = parser.add_argument_group('main arguments')
        main_group.add_argument('-m', '--metadata-file', type=str, help="Input metadata file in simple tabular format with samples in rows and metadata fields in columns. QIIME 2 metadata format is also accepted, with an extra row to define categorical and numerical fields. If not provided and --input-file is a .biom files, will attempt to get metadata from it. ")
        main_group.add_argument('-t', '--taxonomy', type=str, default=None, help="Define taxonomy to convert entry and annotate samples. Will automatically download and parse or files can be provided with --tax-files.", choices=["ncbi", "gtdb", "silva", "greengenes", "ott"])
        main_group.add_argument('-b', '--tax-files', nargs="*", type=str, default=[], help="Optional specific taxonomy files to use.")
        main_group.add_argument('-r', '--ranks', nargs="*", default=[Config.default_rank_name], type=str, help="Taxonomic ranks to generate visualizations. Use '" + Config.default_rank_name + "' to use entries from the table directly. Default: " + Config.default_rank_name)
        main_group.add_argument('-c', '--config', type=str, help="Configuration file with definitions of references, controls and external tools.")

        output_group = parser.add_argument_group('output arguments')
        output_group.add_argument('-g', '--mgnify', default=False, action='store_true', help="Plot MGnify chart")
        output_group.add_argument('-d', '--decontam', default=False, action='store_true', help="Run and plot DECONTAM")
        output_group.add_argument('-l', '--title', type=str, default="", help="Title to display on the header of the report.")
        output_group.add_argument('-p', '--output-plots', nargs="*", type=str, default=Config.output_plots, help="Plots to generate. Default: " + ",".join(Config.output_plots), choices=Config.output_plots)
        output_group.add_argument('-o', '--output-html', type=str, default="output.html", help="File to output report. Default: output.html")
        output_group.add_argument('--full-offline', default=False, action='store_true', help="Embed javascript library in the output file. File will be around 1.5MB bigger but also work without internet connection. That way your report will live forever.")

        data_group = parser.add_argument_group('general data options')
        data_group.add_argument('-f', '--level-separator', default=None, type=str, help="If provided, consider --input-table to be a hierarchical multi-level table where the observations headers are separated by the indicated separator characther (usually ';' or '|')")
        data_group.add_argument('-y', '--values', default=None, type=str, help="Force 'count' or 'normalized' data parsing. Empty to auto-detect.")
        data_group.add_argument('-w', '--cumm-levels', default=False, action='store_true', help="Activate if input table has already cummulative values among levels.")
        data_group.add_argument('-s', '--transpose', default=False, action='store_true', help="Transpose --input-table (if samples are listed on columns and observations on rows)")
        data_group.add_argument('-u', '--unassigned-header', nargs="*", type=str, default=None, help="Define one or more header names containing unsassinged/unclassified counts.")
        data_group.add_argument('--obs-replace', nargs="*", type=str, default=[], help="Replace values on table observations labels/headers (support regex). Example: '_' ' ' will replace underscore with spaces, '^.+__' '' will remove the matching regex.")
        data_group.add_argument('--sample-replace', nargs="*", type=str, default=[], help="Replace values on table sample labels/headers (support regex). Example: '_' ' ' will replace underscore with spaces, '^.+__' '' will remove the matching regex.")
        data_group.add_argument('-z', '--replace-zeros', type=str, default="1000", help="INT (add 'smallest count'/INT to every raw count), FLOAT (add FLOAT to every raw count). Default: 1000")
        data_group.add_argument('--min-frequency', type=float, help="Define minimum number/percentage of samples containing an observation to keep the observation [values between 0-1 for percentage, >1 specific number].")
        data_group.add_argument('--max-frequency', type=float, help="Define maximum number/percentage of samples containing an observation to keep the observation [values between 0-1 for percentage, >1 specific number].")
        data_group.add_argument('--min-count', type=float, help="Define minimum number/percentage of counts to keep an observation [values between 0-1 for percentage, >1 specific number].")
        data_group.add_argument('--max-count', type=float, help="Define maximum number/percentage of counts to keep an observation [values between 0-1 for percentage, >1 specific number].")

        sample_group = parser.add_argument_group('Samples options')
        sample_group.add_argument('-j', '--top-obs-bars', type=int, default=20, help="Top abundant observations to show in the bars.")

        heatmap_group = parser.add_argument_group('Heatmap and clustering options')
        heatmap_group.add_argument('-a', '--transformation', type=str, default="log", help="none (counts), norm (percentage), log (log10), clr (centre log ratio). Default: log")
        heatmap_group.add_argument('-e', '--metadata-cols', type=int, default=3, help="How many metadata cols to show on the heatmap. Higher values makes plot slower to navigate.")
        heatmap_group.add_argument('--optimal-ordering', default=False, action='store_true', help="Activate optimal_ordering on linkage, takes longer for large number of samples.")
        heatmap_group.add_argument('--show-zeros', default=False, action='store_true', help="Do not skip zeros on heatmap. File will be bigger and iteraction with heatmap slower.")
        heatmap_group.add_argument('--linkage-methods', type=str, nargs="*", default=["complete"], choices=list(_LINKAGE_METHODS))
        heatmap_group.add_argument('--linkage-metrics', type=str, nargs="*", default=["euclidean"], choices=_METRICS_NAMES)
        heatmap_group.add_argument('--skip-dendrogram', default=False, action='store_true', help="Disable dendogram. Will create smaller files.")

        correlation_group = parser.add_argument_group('Correlation options')
        correlation_group.add_argument('-x', '--top-obs-corr', type=int, default=50, help="Top abundant observations to build the correlationn matrix, based on the avg. percentage counts/sample. 0 for all")

        parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + Config.version)
        parser.add_argument('-D', '--debug', default=False, action='store_true', help=argparse.SUPPRESS)

        return parser.parse_args(argv)
