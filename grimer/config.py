 #!/usr/bin/env python3
import argparse
import sys
from scipy.spatial.distance import _METRICS_NAMES
from scipy.cluster.hierarchy import _LINKAGE_METHODS


class Config:

    version = "1.0.0-alpha1"
    default_rank_name = "default"

    def __new__(self, argv=None):

        parser = argparse.ArgumentParser(description='grimer')

        parser.add_argument('-i', '--input-file', required=True, type=str, help="Main input table with counts (Observation table, Count table, Contingency Tables, ...) or .biom file. By default rows contain observations and columns contain samples (use --tranpose if your file is reversed). First column and first row are used as headers.")
        parser.add_argument('-c', '--config', type=str, help="Configuration file")
        parser.add_argument('-m', '--metadata-file', type=str, help="Input metadata file in simple tabular format. Sample identifiers will be matched with ones provided by --input-table. QIIME 2 metadata format is also accepted, with categorical and numerical fields.")
        parser.add_argument('-t', '--tax', type=str, default=None, help="Define taxonomy to use. By default, do not use any taxonomy.", choices=["ncbi", "gtdb", "silva", "greengenes", "ott"])
        parser.add_argument('-b', '--tax-files', nargs="*", type=str, default=[], help="Taxonomy files. If not provided, will automatically be downloaded.")
        parser.add_argument('-z', '--replace-zeros', type=str, default="1000", help="INT (add 'smallest count'/INT to every raw count), FLOAT (add FLOAT to every raw count). Default: 1000")
        parser.add_argument('-r', '--ranks', nargs="*", default=[Config.default_rank_name], type=str, help="Taxonomic ranks to generate visualizations. Use '" + Config.default_rank_name + "' to use entries from the table directly. Default: " + Config.default_rank_name)
        parser.add_argument('-l', '--title', type=str, default="", help="Title to display on the header of the report.")
        parser.add_argument('-o', '--output-html', type=str, default="output.html", help="File to output report. Default: output.html")
        parser.add_argument('--full-offline', default=False, action='store_true', help="Embed javascript library in the output file. File will be around 1.5MB bigger but also work without internet connection. That way your report will live forever.")

        table_group = parser.add_argument_group('Table options')
        table_group.add_argument('-f', '--level-separator', default=None, type=str, help="If provided, consider --input-table to be a hiearchical multi-level table where the observations headers are separated by the indicated separator characther (usually ';' or '|')")
        table_group.add_argument('-y', '--values', default=None, type=str, help="Force 'count' or 'normalized' data parsing. Empty to auto-detect.")
        table_group.add_argument('-s', '--transpose', default=False, action='store_true', help="Transpose --input-table (if samples are listed on columns and observations on rows)")
        table_group.add_argument('-u', '--unassigned-header', nargs="*", type=str, default=None, help="Define one or more header names containing unsassinged/unclassified counts.")
        table_group.add_argument('--obs-replace', nargs="*", type=str, default=[], help="Replace values on table observations labels/headers (support regex). Example: '_' ' ' will replace underscore with spaces, '^.+__' '' will remove the matching regex.")
        table_group.add_argument('--sample-replace', nargs="*", type=str, default=[], help="Replace values on table sample labels/headers (support regex). Example: '_' ' ' will replace underscore with spaces, '^.+__' '' will remove the matching regex.")

        filter_group = parser.add_argument_group('Observation filter options')
        filter_group.add_argument('--min-frequency', type=float, help="Define minimum number/percentage of samples containing an observation to keep the observation [values between 0-1 for percentage, >1 specific number].")
        filter_group.add_argument('--max-frequency', type=float, help="Define maximum number/percentage of samples containing an observation to keep the observation [values between 0-1 for percentage, >1 specific number].")
        filter_group.add_argument('--min-count', type=float, help="Define minimum number/percentage of counts to keep an observation [values between 0-1 for percentage, >1 specific number].")
        filter_group.add_argument('--max-count', type=float, help="Define maximum number/percentage of counts to keep an observation [values between 0-1 for percentage, >1 specific number].")

        overview_group = parser.add_argument_group('Overview options')
        overview_group.add_argument('-g', '--mgnify', default=False, action='store_true', help="Use MGNify data")
        overview_group.add_argument('-d', '--decontam', default=False, action='store_true', help="Run DECONTAM")

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

        bars_group = parser.add_argument_group('Bars options')
        bars_group.add_argument('-j', '--top-obs-bars', type=int, default=20, help="Top abundant observations to show in the bars.")

        parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + Config.version)
        parser.add_argument('-D', '--debug', default=False, action='store_true', help=argparse.SUPPRESS)

        return parser.parse_args(argv)
