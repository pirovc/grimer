#!/usr/bin/env python3
_debug = False

#General
import sys

#Internal
from grimer.callbacks import *
from grimer.cds import *
from grimer.config import Config
from grimer.layout import *
from grimer.plots import *
from grimer.func import *

#Bokeh
from bokeh.io import save
from bokeh.plotting import output_file


def main(argv=sys.argv[1:]):
    """
    GRIMER code overview
    1) Load data/analysis: parse configuration, load files and run analysis into data objects
        e.g. args.input_file to Table() and decontam
    2) Generata data sources: Convert objects and analysis int cds/dict
        e.g. table to cds_m_obstable
    3) Plot elements: plot figures and widgets based on cds/dict (and some objects)
        e.g cds_m_obstable to ele["obstable"]["fig"]
    4) Link javascript: link data sources and javascript custom callbacks
    5) Draw layout: Put elements into layout scheme and generate report
    """

    # Parse CLI arguments
    args = Config(argv)
    print_logo_cli(Config.version)
    # Setup global _debug variable to be used by other files with #from grimer.grimer import _debug
    global _debug
    _debug = args.debug

    # 1) Load data/analysis
    cfg = None
    tax = None
    table = None
    metadata = None
    references = None
    controls = None
    control_samples = None
    hcluster = None
    dendro = None
    corr = None

    print_log("- Parsing configuration file")
    cfg = parse_config_file(args.config)

    print_log("- Parsing taxonomy")
    tax = parse_taxonomy(args.tax, args.tax_files)

    print_log("- Parsing input table")
    table = parse_table(args, tax)

    print_log("- Parsing metadata")
    metadata = parse_metadata(args, table)

    print_log("- Parsing references")
    references = parse_references(cfg, tax, args.tax, table.ranks())

    print_log("- Parsing controls")
    controls, control_samples = parse_controls(cfg, table)

    print_log("- Parsing MGnify database")
    mgnify = parse_mgnify(args.mgnify, cfg, tax, table.ranks())

    print_log("- Running DECONTAM")
    decontam = run_decontam(args.decontam, cfg, table, metadata, control_samples)

    print_log("- Running hiearchical clustering")
    hcluster, dendro = run_hclustering(table, args.linkage_methods, args.linkage_metrics, args.transformation, args.skip_dendrogram, args.optimal_ordering)

    print_log("- Running correlation")
    corr = run_correlation(table, args.top_obs_corr)

    # 2) Generata data sources:
    # cds (ColumnDataSource) and dict containers: data structures loaded and parsed by bokehjs
    # "cds" for matrix like dataframes with fixed column sizes
    # "dict" for variable column sizes
    # _p_ : plot -> direct source of figures either pre-loaded or empty
    # _d_ : data -> auxiliar containers to be used/shared among plots
    #               usually by copying and/or transforming values into a _p_ container
    # _m_ : mixed -> contain both plot and data properties

    print_log("- Generating data sources")
    # _m_
    # df: index (unique observations), col|...,  tax|..., aux|ref
    cds_m_obstable = cds_obstable(table, tax, references, controls, control_samples, decontam)
    # _p_
    # df: index (unique sample-ids), aux|..., bar|..., tax|...
    cds_p_samplebars = cds_samplebars(table)
    # stacked: index (repeated observations), rank, ref, direct, parent
    cds_p_references = cds_plot_references(table, tax, references)
    # matrix: index (unique sample-ids), concentrations, controls, counts
    cds_p_decontam = cds_plot_decontam(decontam) if decontam else None
    # {x: [min,max], y_cont: [None,None], y_noncont: [None,None]}
    cds_p_decontam_models = cds_plot_decontam_models(decontam) if decontam else None
    # stacked: index (taxa, level, lineage), count, perc
    cds_p_mgnify = cds_mgnify(mgnify, table, tax) if mgnify else None
    # stacked: index (repeated sample-ids), obs, rank, ov, tv
    cds_p_heatmap = cds_heatmap(table, args.transformation, args.show_zeros)
    # matrix: index (unique sample-ids), md0, md1, ..., md(args.metadata_cols) -> (metadata field, metadata values)
    cds_p_metadata = cds_plot_metadata(metadata, args.metadata_cols) if metadata else None
    # stacked: index (repeated observations), rank, annot
    cds_p_annotations = cds_annotations(table, references, controls, decontam, control_samples)
    # empty matrix {"x": [], "y": [], "c": []}
    cds_p_dendro_x, cds_p_dendro_y = cds_plot_dendro() if not args.skip_dendrogram else [None, None]
    # stacked: index (repeated observations), other observation, rank, rho
    cds_p_correlation = cds_correlation(table, corr)
    # matrix: index (unique sample-ids), 0, 1, ..., top_obs_bars, unassigned, others, factors
    cds_p_obsbars = cds_obsbars(table, args.top_obs_bars)
    # df: index (unique sample-ids), col|...
    cds_p_sampletable = cds_sampletable(table)
    # _d_
    # df: index (unique sample-ids), aux|..., cnt|...,
    cds_d_samples = cds_samples(table, references, controls, decontam)
    # matrix: index (unique sample-ids) x columns (metadata fields) -> metadata values
    cds_d_metadata = cds_metadata(metadata) if metadata else None
    # {taxid: (contam_y1, contam_y2, non_contam_y, pval)}
    cds_d_decontam = cds_decontam(decontam, table.ranks()) if decontam else None
    # key = rank + "|" + method + "|" + metric
    # y: {"default": sorted sample-ids, key: sorted sample-ids, ...}
    # x: {"default|rank": sorted sample-ids, key: sorted sample-ids, ...}
    dict_d_hcluster_x, dict_d_hcluster_y = dict_hcluster(table, hcluster)
    # {key+"|x": x-values, key+"|y": y-values , key+"|c": colors}
    dict_d_dedro_x, dict_d_dedro_y = dict_dendro(table, dendro) if not args.skip_dendrogram else [None, None]
    # {taxid: name}
    dict_d_taxname = dict_taxname(tax, [txid for rank in table.ranks() for txid in table.observations(rank)])
    # {rank: [taxid1,taxid2, ..., taxid(top_obs_bars)]}
    dict_d_topobs = dict_topobs(table, args.top_obs_bars)
    # {taxid: {source: {desc: [refs]}}
    dict_d_refs = dict_refs(table, references)
    # dict: {rank: {obs: {sample: count}}}
    dict_d_sampleobs = dict_sampleobs(table)

    # 3) Plot elements
    print_log("- Plotting elements")

    # Defined fixed layout and plot sizes
    sizes = {}
    sizes["overview_top_panel_height"] = 300
    sizes["overview_top_panel_width_left"] = 250
    sizes["overview_top_panel_width_right"] = 450

    # Elements to plot
    # ele[name]["fig"] -> main figure/element
    # ele[name]["filter"] -> filter to the figure
    # ele[name]["wid"][widget1] -> widgets to the figure
    ele = {}

    # obstable
    ele["obstable"] = {}
    ele["obstable"]["fig"], ele["obstable"]["filter"] = plot_obstable(sizes, cds_m_obstable, table.ranks(), references, controls)
    ele["obstable"]["wid"] = plot_obstable_widgets(sizes, dict_d_taxname, max(cds_m_obstable.data["col|total_counts"]))

    # infopanel
    ele["infopanel"] = {}
    ele["infopanel"]["textarea"] = plot_infopanel()

    # references
    ele["references"] = {}
    if references:
        ele["references"]["fig"], ele["references"]["filter"] = plot_references(sizes, table, cds_p_references, dict_d_taxname)
    else:
        ele["references"]["fig"], ele["references"]["filter"] = None, None
    ele["references"]["wid"] = plot_references_widgets(sizes, references)

    # mgnify
    ele["mgnify"] = {}
    if cds_p_mgnify:
        ele["mgnify"]["fig"], ele["mgnify"]["filter"] = plot_mgnify(sizes, cds_p_mgnify)
    else:
        ele["mgnify"]["fig"], ele["mgnify"]["filter"] = None, None
    ele["mgnify"]["wid"] = plot_mgnify_widgets()

    # decontam
    ele["decontam"] = {}
    ele["decontam"]["wid"] = {}
    if decontam:
        ele["decontam"]["fig"] = plot_decontam(sizes, cds_p_decontam, cds_p_decontam_models, table.get_min_valid_count_perc())
    else:
        ele["decontam"]["fig"] = None
    ele["decontam"]["wid"] = plot_decontam_widgets(sizes)

    # samplebars
    ele["samplebars"] = {}
    ele["samplebars"]["fig"], ele["samplebars"]["legend_obs"], ele["samplebars"]["legend_bars"] = plot_samplebars(cds_p_samplebars, table)
    ele["samplebars"]["wid"] = plot_samplebars_widgets(table.ranks(), metadata, references, controls, decontam, table.normalized)

    # sampletable
    ele["sampletable"] = {}
    ele["sampletable"]["fig"] = plot_sampletable(cds_p_sampletable, sizes, table.ranks())
    ele["sampletable"]["wid"] = plot_sampletable_widgets(sizes, max(cds_p_sampletable.data["col|total"]), metadata)

    # heatmap
    tools_heatmap = "hover,save,box_zoom,reset,crosshair,box_select"
    ele["heatmap"] = {}
    ele["heatmap"]["fig"] = plot_heatmap(table, cds_p_heatmap, tools_heatmap, args.transformation, dict_d_taxname)
    ele["heatmap"]["wid"] = plot_heatmap_widgets(table.ranks(), args.linkage_methods, args.linkage_metrics, references, controls, metadata, decontam)

    # metadata (heatmap)
    ele["metadata"] = {}
    ele["metadata"]["wid"] = {}
    if metadata:
        ele["metadata"]["fig"], ele["metadata"]["wid"] = plot_metadata(ele["heatmap"]["fig"], tools_heatmap, metadata, cds_d_metadata, cds_p_metadata)
    else:
        ele["metadata"]["fig"] = Spacer()
        ele["metadata"]["wid"]["metadata_multiselect"] = Spacer()
        ele["metadata"]["wid"]["legend_colorbars"] = Spacer()
        ele["metadata"]["wid"]["toggle_legend"] = Spacer()

    # annotations
    ele["annotations"] = {}
    if cds_p_annotations.data["index"].size:
        ele["annotations"]["fig"] = plot_annotations(ele["heatmap"]["fig"], tools_heatmap, cds_p_annotations, dict_d_taxname)
    else:
        ele["annotations"]["fig"] = Spacer()

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
    ele["correlation"]["fig"], ele["correlation"]["filter"] = plot_correlation(cds_p_correlation, table.ranks(), dict_d_taxname)
    ele["correlation"]["wid"] = plot_correlation_widgets(table.ranks(), args.top_obs_corr)

    # obsbars
    ele["obsbars"] = {}
    ele["obsbars"]["wid"] = plot_obsbars_widgets(table.ranks(), metadata, dict_d_topobs, dict_d_taxname, args.top_obs_bars)
    ele["obsbars"]["fig"], ele["obsbars"]["legend"] = plot_obsbars(cds_p_obsbars, dict_d_topobs, table.ranks(), args.top_obs_bars, dict_d_taxname, ele["obsbars"]["wid"]["rank_select"])

    #4) Link javascript:
    print_log("- Linking javascript")

    link_obstable_filter(ele, cds_m_obstable, table.ranks())
    link_obstable_samplebars(ele,
                             cds_m_obstable,
                             cds_p_samplebars,
                             cds_d_samples,
                             dict_d_sampleobs,
                             cds_d_metadata,
                             cds_p_decontam,
                             cds_p_decontam_models,
                             cds_d_decontam,
                             cds_p_references,
                             table.ranks(),
                             table.get_min_valid_count_perc(),
                             table.get_total().max(),
                             cds_p_mgnify,
                             dict_d_refs,
                             dict_d_taxname)
    link_heatmap_widgets(ele,
                         cds_d_samples,
                         cds_d_metadata,
                         cds_p_metadata,
                         dict_d_hcluster_x,
                         dict_d_hcluster_y,
                         cds_p_dendro_x,
                         cds_p_dendro_y,
                         dict_d_dedro_x,
                         dict_d_dedro_y,
                         cds_p_annotations,
                         cds_m_obstable,
                         cds_p_heatmap,
                         table.ranks(),
                         dict_d_taxname)
    link_metadata_widgets(ele, cds_p_metadata, cds_d_metadata, args.metadata_cols)
    link_correlation_widgets(ele, cds_p_correlation)
    link_obsbars_widgets(ele,
                         cds_p_obsbars,
                         dict_d_topobs,
                         dict_d_sampleobs,
                         cds_d_samples,
                         args.top_obs_bars,
                         dict_d_taxname,
                         cds_d_metadata,
                         cds_p_sampletable)
    link_sampletable_select(ele, cds_p_sampletable, cds_d_metadata)

    # 5) Draw layout
    print_log("- Drawing layout")
    # Define path of running script to get static files
    script_dir, _ = os.path.split(__file__)
    logo_path = os.path.join(script_dir, "img", "logo.png")

    final_layout = make_layout(ele, sizes, Config.version, logo_path, args.title)

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
    print_log("- Saving report")
    output_file(args.output_html, title="GRIMER" if not args.title else "GRIMER - " + args.title, mode=mode)
    save(final_layout, template=template)
    print_log("File: " + args.output_html)
    file_size_bytes = os.path.getsize(args.output_html)
    print_log("Size: " + str(file_size_bytes) + " bytes (" + '{0:.2f} MB'.format(file_size_bytes / float(1024 ** 2)) + ")")

if __name__ == "__main__":
    main()
