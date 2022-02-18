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
from grimer.config import Config
from grimer.layout import *
from grimer.plots import *
from grimer.utils import *

# MultiTax
from multitax import *

#Bokeh
from bokeh.io import save
from bokeh.plotting import output_file


def main(argv=sys.argv[1:]):

    args = Config(argv)
    print_logo_cli(Config.version)

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
        print_log("- No taxonomy set")
    print_log("")

    # Table of counts
    print_log("- Parsing table")
    if not args.ranks:
        args.ranks = [Config.default_rank_name]

    if args.input_file.endswith(".biom"):
        args.level_separator = ";"
        args.transpose = True

    table_df, total, unassigned = parse_input_table(args.input_file, args.unassigned_header, args.transpose, args.sample_replace)
    if args.level_separator:
        ranked_tables, lineage = parse_multi_table(table_df, args.ranks, tax, args.level_separator, args.obs_replace)
    else:
        ranked_tables, lineage = parse_single_table(table_df, args.ranks, tax, Config.default_rank_name)

    if not ranked_tables:
        print_log("Could not parse input table")
        return 1

    table = Table(table_df.index, total, unassigned)
    table.lineage = lineage

    print_log("")
    print_log("Total valid samples: " + str(len(table.samples)))
    # Check for long sample headers, break some plots
    long_sample_headers = [h for h in table_df.index if len(h) > 70]
    if long_sample_headers:
        print_log("Long sample labels/headers detected, plots may break: ")
        print_log("\n".join(long_sample_headers))
    print_log("")

    for r, t in ranked_tables.items():
        print_log("--- " + r + " ---")
        filtered_trimmed_t = trim_table(filter_input_table(t, total, args.min_frequency, args.max_frequency, args.min_count, args.max_count))
        if t.empty:
            print_log("No valid entries, skipping")
        else:
            # Trim table for empty zeros rows/cols
            table.add_rank(r, filtered_trimmed_t)
            print_log("Total valid observations: " + str(len(table.observations(r))))

    print_log("")
    print_log("Total assigned (sum): " + str(table.total.sum() - table.unassigned.sum()))
    print_log("Total unassigned (sum): " + str(table.unassigned.sum()))
    print_log("")

    # Zero replacement
    try:
        replace_zero_value = table_df[table_df.gt(0)].min().min() / int(args.replace_zeros)
    except:
        replace_zero_value = float(args.replace_zeros)
    if replace_zero_value == 1 and args.transformation == "log":
        replace_zero_value = 0.999999  # Do not allow value 1 using log

    # Metadata
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

    # References (only possible with ncbi identifiers)
    references = {}
    if "references" in cfg and args.tax == "ncbi":
        print_log("- Parsing references")
        references = parse_references(cfg, tax, table.ranks())
        print_log("")

    controls, control_samples = [{}, {}]
    if "controls" in cfg:
        print_log("- Parsing controls")
        # Controls
        controls, control_samples = parse_controls(cfg, table)
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
        mgnify = MGnify(cfg["external"]["mgnify"], ranks=table.ranks() if args.ranks != [Config.default_rank_name] else [])
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
    cds_p_obstable = generate_cds_obstable(table, tax, references, controls, control_samples, decontam)
    # df: index (unique sample-ids), aux|..., bar|..., tax|...
    cds_p_samplebars = generate_cds_samplebars(table)
    # stacked: index (repeated observations), rank, ref, direct, parent
    cds_p_references = generate_cds_plot_references(table, tax, references)
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
    cds_p_annotations = generate_cds_annotations(table, references, controls, decontam)
    # empty matrix {"x": [], "y": [], "c": []}
    cds_p_dendro_x, cds_p_dendro_y = generate_cds_plot_dendro() if not args.skip_dendrogram else [None, None]
    # stacked: index (repeated observations), other observation, rank, rho
    cds_p_correlation = generate_cds_correlation(table, args.top_obs_corr, replace_zero_value)
    # matrix: index (unique sample-ids), 0, 1, ..., top_obs_bars, unassigned, others, factors
    cds_p_obsbars = generate_cds_obsbars(table, args.top_obs_bars)
    # df: index (unique sample-ids), col|...
    cds_p_sampletable = generate_cds_sampletable(table)

    # _d_
    # dict: {rank: {obs: {sample: count}}}
    dict_d_sampleobs = generate_dict_sampleobs(table)
    # df: index (unique sample-ids), aux|..., cnt|...,
    cds_d_samples = generate_cds_samples(table, references, controls, decontam)
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
    dict_d_refs = generate_dict_refs(table, references)

    ############ PLOT ELEMENTS (Figures, Widgets, ...)
    ############ "fig": main figure
    ############ "wid": widgets

    # Layout and plot sizes
    sizes = {}
    sizes["overview_top_panel_height"] = 300
    sizes["overview_top_panel_width_left"] = 250
    sizes["overview_top_panel_width_right"] = 450

    ele = {}

    # obstable
    ele["obstable"] = {}
    ele["obstable"]["fig"], ele["obstable"]["widgets_filter"] = plot_obstable(sizes, cds_p_obstable, table.ranks(), references.keys(), controls.keys())
    ele["obstable"]["wid"] = plot_obstable_widgets(sizes, dict_d_taxname, max(cds_p_obstable.data["col|total_counts"]))

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
        ele["decontam"]["fig"] = plot_decontam(sizes, cds_p_decontam, cds_p_decontam_models, min_obs_perc)
    else:
        ele["decontam"]["fig"] = None
    ele["decontam"]["wid"] = plot_decontam_widgets()

    # samplebars
    ele["samplebars"] = {}
    ele["samplebars"]["fig"], ele["samplebars"]["legend_obs"], ele["samplebars"]["legend_bars"] = plot_samplebars(cds_p_samplebars, max_total_count, table.ranks())
    ele["samplebars"]["wid"] = plot_samplebars_widgets(table.ranks(), metadata, list(references.keys()), list(controls.keys()), decontam)

    # sampletable
    ele["sampletable"] = {}
    ele["sampletable"]["fig"] = plot_sampletable(cds_p_sampletable, sizes, table.ranks())
    ele["sampletable"]["wid"] = plot_sampletable_widgets(sizes, max(cds_p_sampletable.data["col|total"]), metadata)

    # heatmap
    tools_heatmap = "hover,save,box_zoom,reset,crosshair,box_select"
    ele["heatmap"] = {}
    ele["heatmap"]["fig"] = plot_heatmap(table, cds_p_heatmap, tools_heatmap, args.transformation, dict_d_taxname)
    ele["heatmap"]["wid"] = plot_heatmap_widgets(table.ranks(), args.linkage_methods, args.linkage_metrics, list(references.keys()), list(controls.keys()), metadata, decontam)

    # metadata (heatmap)
    ele["metadata"] = {}
    ele["metadata"]["wid"] = {}
    if metadata:
        ele["metadata"]["fig"], ele["metadata"]["wid"] = plot_metadata(ele["heatmap"]["fig"], tools_heatmap, metadata, cds_d_metadata, cds_p_metadata)
    else:
        ele["metadata"]["fig"] = Spacer()
        ele["metadata"]["wid"]["metadata_multiselect"] = Spacer()
        ele["metadata"]["wid"]["legend_colorbars"] = Spacer()

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
    ele["correlation"]["fig"], ele["correlation"]["rho_filter"] = plot_correlation(cds_p_correlation, table.ranks(), dict_d_taxname)
    ele["correlation"]["wid"] = plot_correlation_widgets(table.ranks(), args.top_obs_corr)

    # obsbars
    ele["obsbars"] = {}
    ele["obsbars"]["wid"] = plot_obsbars_widgets(table.ranks(), metadata, dict_d_topobs, dict_d_taxname, args.top_obs_bars)
    ele["obsbars"]["fig"], ele["obsbars"]["legend"] = plot_obsbars(cds_p_obsbars, dict_d_topobs, table.ranks(), args.top_obs_bars, dict_d_taxname, ele["obsbars"]["wid"]["rank_select"])

    ############ JAVASCRIPT LINKING

    link_obstable_filter(ele, cds_p_obstable, table.ranks())

    link_obstable_samplebars(ele,
                             cds_p_obstable,
                             cds_p_samplebars,
                             cds_d_samples,
                             dict_d_sampleobs,
                             cds_d_metadata,
                             cds_p_decontam,
                             cds_p_decontam_models,
                             cds_d_decontam,
                             cds_p_references,
                             table.ranks(),
                             min_obs_perc,
                             max_total_count,
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
                         cds_p_obstable,
                         cds_p_heatmap,
                         table.ranks(),
                         dict_d_taxname)

    link_metadata_widgets(ele, cds_p_metadata, cds_d_metadata, max_metadata_cols)

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

    ############ LAYOUT

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
    output_file(args.output_html, title="GRIMER" if not args.title else "GRIMER - " + args.title, mode=mode)
    save(final_layout, template=template)
    print_log("File: " + args.output_html)
    file_size_bytes = os.path.getsize(args.output_html)
    print_log("Size: " + str(file_size_bytes) + " bytes (" + '{0:.2f} MB'.format(file_size_bytes / float(1024 ** 2)) + ")")

if __name__ == "__main__":
    main()
