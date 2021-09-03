#!/usr/bin/env python3
_debug = False

#General
import argparse
import yaml

#Internal
from grimer.table import Table
from grimer.metadata import Metadata
from grimer.mgnify import MGnify
from grimer.callbacks import *
from grimer.cds import *
from grimer.layout import *
from grimer.plots import *
from grimer.utils import *

# MultiTax
from multitax import *

#Bokeh
from bokeh.io import save
from bokeh.plotting import output_file

# Scipy
from scipy.spatial.distance import _METRICS_NAMES
from scipy.cluster.hierarchy import _LINKAGE_METHODS


def main():

    default_rank_name = "default"

    version = "1.0.0-alpha1"
    parser = argparse.ArgumentParser(description='grimer')
    parser.add_argument('-i', '--input-file', required=True, type=str, help="Main input table with counts (Observation table, Count table, Contingency Tables, ...) or .biom file. By default rows contain observations and columns contain samples (use --tranpose if your file is reversed). First column and first row are used as headers.")
    parser.add_argument('-c', '--config', required=True, type=str, help="Configuration file")
    parser.add_argument('-m', '--metadata', type=str, help="Input metadata file in simple tabular format. Sample identifiers will be matched with ones provided by --input-table. QIIME 2 metadata format is also accepted, with categorical and numerical fields.")
    parser.add_argument('-t', '--tax', type=str, default=None, help="Define taxonomy to use. By default, do not use any taxonomy.", choices=["ncbi", "gtdb", "silva", "greengenes", "ott"])
    parser.add_argument('-b', '--tax-files', nargs="*", type=str, default=None, help="Taxonomy files. If not provided, will automatically be downloaded.")
    parser.add_argument('-z', '--replace-zeros', type=str, default="1000", help="INT (add 'smallest count'/INT to every raw count), FLOAT (add FLOAT to every raw count). Default: 1000")
    parser.add_argument('-r', '--ranks', nargs="*", default=[default_rank_name], type=str, help="Taxonomic ranks to generate visualizations. Use '" + default_rank_name + "' to use entries from the table directly. Default: " + default_rank_name)
    parser.add_argument('-l', '--title', type=str, default="", help="Title to display on the header of the report.")
    parser.add_argument('-o', '--output-html', type=str, default="output.html", help="File to output report. Default: output.html")
    parser.add_argument('--full-offline', default=False, action='store_true', help="Embed javascript library in the output file. File will be around 1.5MB bigger but also work without internet connection. That way your report will live forever.")

    table_group = parser.add_argument_group('Table options')
    table_group.add_argument('-f', '--level-separator', default=None, type=str, help="If provided, consider --input-table to be a hiearchical multi-level table where the observations headers are separated by the indicated separator characther (usually ';' or '|')")
    table_group.add_argument('-s', '--transpose', default=False, action='store_true', help="Transpose --input-table (if samples are listed on columns and observations on rows)")
    table_group.add_argument('-u', '--unassigned-header', nargs="*", type=str, default=None, help="Define one or more header names containing unsassinged/unclassified counts.")
    table_group.add_argument('--obs-replace', nargs="*", type=str, default=[], help="Replace values on observations headers usin (support regex). Example: '_' ' ' will replace underscore with spaces, '^.+__' '' will remove the matching regex.")

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
    heatmap_group.add_argument('-e', '--metadata-cols', type=int, default=5, help="How many metadata cols to show on the heatmap. Higher values makes plot slower to navigate.")
    heatmap_group.add_argument('--optimal-ordering', default=False, action='store_true', help="Activate optimal_ordering on linkage, takes longer for large number of samples.")
    heatmap_group.add_argument('--show-zeros', default=False, action='store_true', help="Do not skip zeros on heatmap. File will be bigger and iteraction with heatmap slower.")
    heatmap_group.add_argument('--linkage-methods', type=str, nargs="*", default=["complete"], choices=list(_LINKAGE_METHODS))
    heatmap_group.add_argument('--linkage-metrics', type=str, nargs="*", default=["euclidean", "braycurtis"], choices=_METRICS_NAMES)
    heatmap_group.add_argument('--skip-dendrogram', default=False, action='store_true', help="Disable dendogram. Will create smaller files.")

    correlation_group = parser.add_argument_group('Correlation options')
    correlation_group.add_argument('-x', '--top-obs-corr', type=int, default=20, help="Top abundant observations to build the correlationn matrix, based on the avg. percentage counts/sample. 0 for all")

    bars_group = parser.add_argument_group('Bars options')
    bars_group.add_argument('-j', '--top-obs-bars', type=int, default=20, help="Top abundant observations to show in the bars.")

    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    parser.add_argument('-D', '--debug', default=False, action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()

    print_logo_cli(version)

    global _debug
    _debug = args.debug

    # Config file
    with open(args.config, 'r') as file:
        cfg = yaml.safe_load(file)

    # Taxonomy
    tax = None
    if args.tax:
        if args.tax_files:
            print_log("- Parsing taxonomy")
        else:
            print_log("- Downloading and parsing taxonomy")
        print_log(args.tax)
        if args.tax == "ncbi":
            tax = NcbiTx(files=args.tax_files, extended_names=True)
        elif args.tax == "gtdb":
            tax = GtdbTx(files=args.tax_files)
        elif args.tax == "silva":
            tax = SilvaTx(files=args.tax_files)
        elif args.tax == "greengenes":
            tax = GreengenesTx(files=args.tax_files)
        elif args.tax == "ott":
            tax = OttTx(files=args.tax_files, extended_names=True)
    else:
        print_log(" - No taxonomy set")
    print_log("")

    # Table of counts
    print_log("- Parsing table")
    if not args.ranks:
        args.ranks = [default_rank_name]

    if args.input_file.endswith(".biom"):
        args.level_separator = ";"
        args.transpose = True

    table_df, total, unassigned = parse_input_table(args.input_file, args.unassigned_header, args.transpose, args.min_frequency, args.max_frequency, args.min_count, args.max_count)
    if args.level_separator:
        ranked_tables, lineage = parse_multi_table(table_df, args.ranks, tax, args.level_separator, args.obs_replace)
    else:
        ranked_tables, lineage = parse_single_table(table_df, args.ranks, tax, default_rank_name)

    if not ranked_tables:
        print_log("Could not parse input table")
        return 1

    table = Table(table_df.index, total, unassigned)
    table.lineage = lineage
    print_log("Samples: " + str(len(table.samples)))
    print_log("Observations: ")
    for r, t in ranked_tables.items():
        print_log(" " + r + ":")
        if t.empty:
            print_log("Skipping without valid entries")
        else:
            # Trim table for empty zeros rows/cols
            table.add_rank(r, trim_table(t))
            print_log("  " + str(len(table.observations(r))) + " observations")

    print_log("")
    print_log("Total assigned (sum): " + str(table.total.sum()))
    print_log("Total unassigned (sum): " + str(table.unassigned.sum()))
    print_log("")

    # Zero replacement
    try:
        replace_zero_value = table_df[table_df.gt(0)].min().min() / int(args.replace_zeros)
    except:
        replace_zero_value = float(args.replace_zeros)

    # Do not allow value 1 using log
    if replace_zero_value == 1 and args.transformation == "log":
        replace_zero_value = 0.999999

    # Parse Metadata
    max_metadata_cols = args.metadata_cols
    if args.metadata:
        print_log("- Parsing metadata")
        metadata = Metadata(args.metadata, samples=table.samples.to_list())
        if metadata.data.empty:
            metadata = None
            print_log("No valid metadata")
        else:
            print_log("Samples: " + str(metadata.data.shape[0]))
            print_log("Numeric Fields: " + str(metadata.get_data("numeric").shape[1]))
            print_log("Categorical Fields: " + str(metadata.get_data("categorical").shape[1]))
            if len(metadata.get_col_headers()) < args.metadata_cols:
                max_metadata_cols = len(metadata.get_col_headers())
        print_log("")
    else:
        metadata = None

    # Sources of contamination/references/controls
    print_log("- Parsing sources (contamination/references/controls)")
    if args.tax == "ncbi":
        contaminants, references = parse_sources(cfg, tax, table.ranks())
    else:
        contaminants, references = [{}, {}]

    controls, control_samples = parse_controls(cfg, tax, table)
    print_log("")

    # Run and load decontam results
    if args.decontam:
        print_log("- Running DECONTAM")
        decontam = run_decontam(cfg, table, metadata, control_samples)
        print_log("")
    else:
        decontam = None

    # Mgnify
    if args.mgnify and "mgnify" in cfg["external"]:
        print_log("- Parsing MGNify")
        mgnify = MGnify(cfg["external"]["mgnify"], ranks=table.ranks() if args.ranks != [default_rank_name] else [])
        if tax:
            mgnify.update_taxids(update_tax_nodes([tuple(x) for x in mgnify.data[["rank", "taxa"]].to_numpy()], tax))
        print_log("")
    else:
        mgnify = None

    # Hiearchical clustering
    print_log("- Running hiearchical clustering")
    hcluster, dendro = run_hclustering(table, args.linkage_methods, args.linkage_metrics, args.transformation, replace_zero_value, args.skip_dendrogram, args.optimal_ordering)
    print_log("")

    # save max/min values to control ranges
    max_total_count = table.total.max()
    min_obs_perc = min([table.get_counts_perc(rank)[table.get_counts_perc(rank) > 0].min().min() for rank in table.ranks()])

    print_log("- Generating GRIMER report")
    ############ cds (ColumnDataSource) and dict containers: data structures loaded and parsed by bokehjs
    ############ "cds" for matrix like dataframes with fixed column sizes
    ############ "dict" for variable column sizes
    ############ _p_ : plot -> direct source of figures
    ############ _d_ : data -> auxiliar containers to be used/shared among plots
    ############               usually by copying and/or transforming values into a _p_ container

    # _p_
    # df: index (unique observations), col|...,  tax|..., aux|ref
    # this cds an exeption and contains data to plot (col|) and auxiliary data (tax|)
    cds_p_obstable = generate_cds_obstable(table, tax, contaminants, references, controls, control_samples, decontam)
    # df: index (unique sample-ids), aux|..., bar|..., tax|...
    cds_p_samplebars = generate_cds_bars(table)
    # matrix: index (unique sample-ids), concentrations, controls, counts
    cds_p_decontam = generate_cds_plot_decontam(decontam) if decontam else None
    # {x: [min,max], y_cont: [None,None], y_noncont: [None,None]}
    cds_p_decontam_models = generate_cds_plot_decontam_models(decontam) if decontam else None
    # stacked: index (taxa, level, lineage), count, perc
    cds_p_mgnify = generate_cds_mgnify(mgnify, table, tax) if mgnify else None
    # stacked: index (repeated sample-ids), obs, rank, ov, tv
    cds_p_heatmap = generate_cds_heatmap(table, args.transformation, replace_zero_value, args.show_zeros)
    # matrix: index (unique sample-ids), md0, md1, ..., md(max_metadata_cols) -> (metadata field, metadata values)
    cds_p_metadata = generate_cds_plot_metadata(metadata, max_metadata_cols) if metadata else None
    # stacked: index (repeated observations), rank, annot
    cds_p_annotations = generate_cds_annotations(table, contaminants, references, controls, decontam)
    # empty matrix {"x": [], "y": [], "c": []}
    cds_p_dendro_x, cds_p_dendro_y = generate_cds_plot_dendro() if not args.skip_dendrogram else [None, None]
    # stacked: index (repeated observations), other observation, rank, rho, pval, pval_corr
    cds_p_correlation = generate_cds_correlation(table, args.top_obs_corr)
    # matrix: index (unique sample-ids), 0, 1, ..., top_obs_bars, unassigned, others, factors
    cds_p_obsbars = generate_cds_obsbars(table, args.top_obs_bars)

    # _d_
    # matrix: index (unique sample-ids), columns (unique observations) -> raw counts
    cds_d_sampleobs = generate_cds_sampleobs(table)
    # df: index (unique sample-ids), aux|..., cnt|...,
    cds_d_samples = generate_cds_samples(table, references, contaminants, controls, decontam)
    # matrix: index (unique sample-ids) x columns (metadata fields) -> metadata values
    cds_d_metadata = generate_cds_metadata(metadata) if metadata else None
    # {taxid: (contam_y1, contam_y2, non_contam_y, pval)}
    cds_d_decontam = generate_cds_decontam(decontam, table.ranks()) if decontam else None
    # key = rank + "|" + method + "|" + metric
    # y: {"default": sorted sample-ids, key: sorted sample-ids, ...}
    # x: {"default|rank": sorted sample-ids, key: sorted sample-ids, ...}
    dict_d_hcluster_x, dict_d_hcluster_y = generate_dict_hcluster(table, hcluster)
    # {key+"|x": x-values, key+"|y": y-values , key+"|c": colors}
    dict_d_dedro_x, dict_d_dedro_y = generate_dict_dendro(table, dendro) if not args.skip_dendrogram else [None, None]
    # {taxid: name}
    dict_d_taxname = generate_dict_taxname(tax, [txid for rank in table.ranks() for txid in table.observations(rank)])
    # {rank: [taxid1,taxid2, ..., taxid(top_obs_bars)]}
    dict_d_topobs = generate_dict_topobs(table, args.top_obs_bars)
    # {taxid: {source: {desc: [refs]}}
    dict_d_refs = generate_dict_refs(table, contaminants, references)

    ############ PLOT ELEMENTS (Figures, Widgets, ...)
    ############ "fig": main figure
    ############ "wid": widgets

    ele = {}

    # obstable
    ele["obstable"] = {}
    ele["obstable"]["fig"], ele["obstable"]["widgets_filter"] = plot_obstable(cds_p_obstable, table.ranks(), contaminants.keys(), controls.keys())
    ele["obstable"]["wid"] = plot_obstable_widgets(dict_d_taxname, max(cds_p_obstable.data["col|total_counts"]))

    # infopanel
    ele["infopanel"] = {}
    ele["infopanel"]["textarea"] = plot_infopanel()

    # mgnify
    ele["mgnify"] = {}
    if cds_p_mgnify:
        ele["mgnify"]["fig"], ele["mgnify"]["filter"] = plot_mgnify(cds_p_mgnify)
    else:
        ele["mgnify"]["fig"], ele["mgnify"]["filter"] = None, None
    ele["mgnify"]["wid"] = plot_mgnify_widgets()

    # decontam
    ele["decontam"] = {}
    ele["decontam"]["wid"] = {}
    if decontam:
        ele["decontam"]["fig"] = plot_decontam(cds_p_decontam, cds_p_decontam_models, min_obs_perc)
    else:
        ele["decontam"]["fig"] = None
    ele["decontam"]["wid"] = plot_decontam_widgets()

    # samplebars
    ele["samplebars"] = {}
    ele["samplebars"]["fig"], ele["samplebars"]["legend_obs"], ele["samplebars"]["legend_bars"] = plot_samplebars(cds_p_samplebars, max_total_count, table.ranks())
    ele["samplebars"]["wid"] = plot_samplebars_widgets(table.ranks(), metadata, list(contaminants.keys()), list(references.keys()), list(controls.keys()), decontam)

    # heatmap
    tools_heatmap = "hover,save,box_zoom,reset,crosshair,box_select"
    ele["heatmap"] = {}
    ele["heatmap"]["fig"] = plot_heatmap(table, cds_p_heatmap, tools_heatmap, args.transformation, dict_d_taxname)
    ele["heatmap"]["wid"] = plot_heatmap_widgets(table.ranks(), args.linkage_methods, args.linkage_metrics, list(contaminants.keys()), list(references.keys()), list(controls.keys()), metadata, decontam)

    # metadata (heatmap)
    ele["metadata"] = {}
    ele["metadata"]["wid"] = {}
    if metadata:
        ele["metadata"]["fig"], ele["metadata"]["wid"] = plot_metadata(ele["heatmap"]["fig"], tools_heatmap, metadata, cds_d_metadata, cds_p_metadata)
    else:
        ele["metadata"]["fig"] = Spacer()
        ele["metadata"]["wid"]["metadata_multiselect"] = Spacer()

    # annotations
    ele["annotations"] = {}
    ele["annotations"]["fig"] = plot_annotations(ele["heatmap"]["fig"], tools_heatmap, cds_p_annotations, dict_d_taxname)

    # dendrograms
    ele["dendrox"] = {}
    ele["dendroy"] = {}
    if not args.skip_dendrogram:
        ele["dendrox"]["fig"], ele["dendroy"]["fig"] = plot_dendrogram(ele["heatmap"]["fig"], tools_heatmap, cds_p_dendro_x, cds_p_dendro_y)
    else:
        ele["dendrox"]["fig"] = Spacer()
        ele["dendroy"]["fig"] = Spacer()

    # correlation
    ele["correlation"] = {}
    ele["correlation"]["fig"], ele["correlation"]["rho_filter"], ele["correlation"]["pval_filter"] = plot_correlation(cds_p_correlation, table.ranks(), dict_d_taxname)
    ele["correlation"]["wid"] = plot_correlation_widgets(table.ranks(), args.top_obs_corr)

    # obsbars
    ele["obsbars"] = {}
    ele["obsbars"]["fig"], ele["obsbars"]["legend"] = plot_obsbars(cds_p_obsbars, dict_d_topobs, table.ranks(), args.top_obs_bars, dict_d_taxname)
    ele["obsbars"]["wid"] = plot_obsbars_widgets(table.ranks(), metadata, dict_d_topobs, dict_d_taxname, args.top_obs_bars)

    ############ JAVASCRIPT LINKING

    link_obstable_filter(ele, cds_p_obstable, table.ranks())

    link_obstable_samplebars(ele,
                             cds_p_obstable,
                             cds_p_samplebars,
                             cds_d_samples,
                             cds_d_sampleobs,
                             cds_d_metadata,
                             cds_p_decontam,
                             cds_p_decontam_models,
                             cds_d_decontam,
                             table.ranks(),
                             min_obs_perc,
                             max_total_count,
                             cds_p_mgnify,
                             dict_d_refs)

    link_heatmap_widgets(ele,
                         cds_d_samples,
                         cds_d_metadata,
                         dict_d_hcluster_x,
                         dict_d_hcluster_y,
                         cds_p_dendro_x,
                         cds_p_dendro_y,
                         dict_d_dedro_x,
                         dict_d_dedro_y,
                         cds_p_annotations,
                         cds_p_obstable,
                         cds_p_heatmap)

    link_metadata_widgets(ele, cds_p_metadata, cds_d_metadata, max_metadata_cols)

    link_correlation_widgets(ele, cds_p_correlation)

    link_obsbars_widgets(ele,
                         cds_p_obsbars,
                         dict_d_topobs,
                         cds_d_sampleobs,
                         cds_d_samples,
                         args.top_obs_bars,
                         dict_d_taxname,
                         cds_d_metadata)

    ############ LAYOUT

    # Define path of running script to get static files
    script_dir, _ = os.path.split(__file__)
    logo_path = os.path.join(script_dir, "img", "logo.png")

    final_layout = make_layout(ele, version, logo_path, args.title)

    template = include_scripts({os.path.join(script_dir, "js", "func.js"): "script",
                                os.path.join(script_dir, "js", "popup.js"): "script",
                                os.path.join(script_dir, "css", "popup.css"): "style"})

    if args.full_offline:
        mode = "inline"  # configure to provide entire Bokeh JS and CSS inline
    elif _debug:
        mode = "absolute-dev"  # non-minimized - configure to load from the installed Bokeh library static directory
    else:
        mode = "cdn"  # configure to load Bokeh JS and CSS from https://cdn.bokeh.org

    # setup output file and JS mode
    output_file(args.output_html, title="GRIMER" if not args.title else "GRIMER - " + args.title, mode=mode)
    save(final_layout, template=template)
    print_log("File: " + args.output_html)
    file_size_bytes = os.path.getsize(args.output_html)
    print_log("Size: " + str(file_size_bytes) + " bytes (" + '{0:.2f} MB'.format(file_size_bytes / float(1024 ** 2)) + ")")

if __name__ == "__main__":
    main()
