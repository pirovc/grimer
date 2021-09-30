import markdown
# Bokeh
from bokeh.models import AdaptiveTicker, Button, CategoricalColorMapper, CDSView, CheckboxGroup, ColorBar, ColumnDataSource, CustomJS, CustomJSHover, Div, FactorRange, FuncTickFormatter, HoverTool, Legend, LinearAxis, LinearColorMapper, MultiChoice, MultiSelect, NumberFormatter, Panel, Paragraph, Range1d, RangeSlider, Select, Spinner, Tabs, TextAreaInput, TextInput
from bokeh.models.filters import IndexFilter, GroupFilter
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.palettes import Blues, Category10, Category20, Colorblind, Dark2, linear_palette, Magma256, Reds, Turbo256
from bokeh.plotting import figure
from bokeh.transform import cumsum, factor_cmap


def plot_samplebars(cds_p_samplebars, max_total_count, ranks):
    # Bar plots has 3 main stacks: selection, others, unassigned
    # stacks can be annotated with references and controls
    samplebars_fig = figure(x_range=FactorRange(factors=cds_p_samplebars.data["aux|factors"]),
                            y_range=Range1d(start=0, end=max_total_count),
                            plot_height=400,
                            sizing_mode="stretch_width",
                            tools="box_zoom,reset,save")

    samplebars_fig.add_tools(HoverTool(
        tooltips=[
            ("Sample", "@index"),
            ("Value", "@$name")
        ],
        mode="mouse",
        point_policy="follow_mouse",
    ))

    fixed_bar_options = ["selected", "others", "unassigned"]
    vbar_ren = samplebars_fig.vbar_stack(["bar|" + f for f in fixed_bar_options],
                                         x="aux|factors",
                                         width=1,
                                         source=cds_p_samplebars,
                                         line_color=None,  # to avoid printing small border for zeros
                                         color=make_color_palette(len(fixed_bar_options)))

    # Second y-axis to plot observations
    samplebars_fig.extra_y_ranges = {"obs": Range1d(start=0, end=100)}
    samplebars_fig.add_layout(LinearAxis(y_range_name="obs"), 'right')

    # Plot obs ranks
    obs_palette = make_color_palette(len(ranks), palette=Dark2)
    legend_obs_items = []
    for i, rank in enumerate(ranks):
        ren = samplebars_fig.scatter(x="aux|factors", y="tax|" + rank,
                                     y_range_name="obs",
                                     name="tax|" + rank, #to work with hover properly
                                     source=cds_p_samplebars,
                                     marker="circle", size=7, line_color="navy", alpha=0.6,
                                     fill_color=obs_palette[i])
        legend_obs_items.append((rank, [ren]))

    # Legend counts (vbars)
    legend_bars_items = [(f, [vbar_ren[i]]) for i, f in enumerate([ranks[0] + "|assigned"] + fixed_bar_options[1:])]
    legend_bars = Legend(items=legend_bars_items)
    legend_bars.margin = 0
    legend_bars.border_line_width = 0
    legend_bars.spacing = 0
    legend_bars.padding = 0
    legend_bars.orientation = "horizontal"
    legend_bars.click_policy = "hide"
    legend_bars.location = "bottom_left"
    samplebars_fig.add_layout(legend_bars, 'above')

    # Legend observations (markers)
    legend_obs = Legend(items=legend_obs_items)
    legend_obs.margin = 0
    legend_obs.border_line_width = 0
    legend_obs.spacing = 0
    legend_obs.padding = 0
    legend_obs.orientation = "horizontal"
    legend_obs.location = "bottom_right"
    legend_obs.click_policy = "hide"
    samplebars_fig.add_layout(legend_obs, "above")

    samplebars_fig.xaxis.major_label_text_font_size = '0pt'
    samplebars_fig.xgrid.grid_line_color = None
    samplebars_fig.xaxis.major_tick_line_color = None
    samplebars_fig.xaxis.minor_tick_line_color = None
    samplebars_fig.yaxis.minor_tick_line_color = None

    samplebars_fig.xaxis.group_label_orientation = "horizontal"
    samplebars_fig.xaxis.subgroup_label_orientation = "vertical"

    samplebars_fig.xaxis.axis_label = "samples"
    samplebars_fig.yaxis[0].axis_label = "# counts"
    samplebars_fig.yaxis[1].axis_label = "% observations"

    return samplebars_fig, legend_obs, legend_bars


def plot_obsbars(cds_p_obsbars, dict_d_topobs, ranks, top_obs_bars, dict_d_taxname):
    obsbars_fig = figure(x_range=FactorRange(factors=cds_p_obsbars.data["factors"]),
                         y_range=Range1d(start=0, end=100),
                         height=450,
                         sizing_mode="stretch_width",
                         tools="box_zoom,reset,hover,save",
                         tooltips=[("Sample", "@index"),
                                   ("Label", "$name"),
                                   ("Value", "@$name")])

    # TODO Need to know which rank to get the correct set of top taxa
    # taxid_name_custom = CustomJSHover(
    #     args=dict(dict_d_taxname=ColumnDataSource(dict(dict_d_taxname=[dict_d_taxname])),
    #               dict_d_topobs=ColumnDataSource(dict(dict_d_topobs=[dict_d_topobs]))),
    #     code='''
    #     //var taxid = dict_d_topobs.data.dict_d_topobs[0][RANK][value];
    #     //console.log(dict_d_topobs);
    #     //return dict_d_taxname.data.dict_d_taxname[0][taxid]; // value holds the @taxid
    #     ''')
    # # Add custom tooltip for heatmap (taxid->name)
    # obsbars.add_tools(HoverTool(
    #     tooltips=[
    #         ('Sample', '@index'),
    #         ('Taxa', '$name{custom}'),
    #         ("Value", "@$name")
    #     ],
    #     formatters={"$name": taxid_name_custom}
    # ))

    bars = [str(i) for i in range(top_obs_bars)] + ["others", "unassigned"]
    # Plot stacked bars with counts
    vbar_ren = obsbars_fig.vbar_stack(bars, x="factors",
                                      source=cds_p_obsbars,
                                      width=1,
                                      #line_color=None,  # to avoid printing small border for zeros
                                      #color=make_color_palette(top_obs_bars, linear=True) + ("#868b8e", "#eeede7"))
                                      color=make_color_palette(top_obs_bars) + ("#868b8e", "#eeede7"))

    obsbars_fig.xaxis.major_label_orientation = "vertical"
    obsbars_fig.xaxis.group_label_orientation = "horizontal"
    obsbars_fig.xaxis.subgroup_label_orientation = "vertical"
    obsbars_fig.xaxis.major_label_text_font_size = "8px"
    obsbars_fig.xgrid.grid_line_color = None
    obsbars_fig.xaxis.axis_label = "samples"
    obsbars_fig.yaxis.axis_label = "% counts"

    # Fixed legend
    fixed_legend_bars_items = []
    fixed_legend_bars_items.append(("others", [vbar_ren[-2]]))
    fixed_legend_bars_items.append(("unassigned", [vbar_ren[-1]]))
    fixed_legend_obsbars = Legend(items=fixed_legend_bars_items)
    fixed_legend_obsbars.border_line_width = 0
    fixed_legend_obsbars.margin = 0
    fixed_legend_obsbars.spacing = 0
    fixed_legend_obsbars.padding = 0
    fixed_legend_obsbars.orientation = "horizontal"
    fixed_legend_obsbars.click_policy = "hide"
    fixed_legend_obsbars.location = "bottom_right"
    obsbars_fig.add_layout(fixed_legend_obsbars, 'above')

    # Legend counts (vbars)
    legend_bars_items = []
    for i in range(top_obs_bars):
        if i < len(dict_d_topobs[ranks[0]]):
            label = str(i) + ") " + dict_d_taxname[dict_d_topobs[ranks[0]][i]]
        else:
            label = None
        legend_bars_items.append((label, [vbar_ren[i]]))

    legend_obsbars = Legend(items=legend_bars_items)
    # legend_bars.label_text_font_size="9px"
    # legend_bars.glyph_height=9
    # legend_bars.glyph_width=9
    # legend_bars.label_height=9
    legend_obsbars.border_line_width = 0
    legend_obsbars.margin = 0
    legend_obsbars.spacing = 0
    legend_obsbars.padding = 0
    legend_obsbars.orientation = "vertical"
    legend_obsbars.click_policy = "hide"
    legend_obsbars.location = "top_left"
    obsbars_fig.add_layout(legend_obsbars, 'right')

    return obsbars_fig, legend_obsbars


def plot_obsbars_widgets(ranks, metadata, dict_d_topobs, dict_d_taxname, top_obs_bars):
    rank_select = Select(title="Taxonomic rank:", value=ranks[0], options=ranks)

    sort_options = {}
    sort_options["Default"] = [("input_order", "input order")]
    sort_options["Default"].append(("col|others", "others"))
    sort_options["Default"].append(("col|unassigned", "unassigned"))

    #sort_options["Observation"] = [("col|" + str(i), dict_d_taxname[t]) for i, t in enumerate(dict_d_topobs[ranks[0]])]
    sort_options["Observation"] = [("col|" + str(i), str(i)) for i in range(top_obs_bars)]

    sort_options["Numeric Metadata"] = []
    if metadata:
        numeric_md_data = metadata.get_data(metadata_type="numeric").columns.to_list()
        if numeric_md_data:
            sort_options["Numeric Metadata"] = [("metadata_num|" + md, md) for md in numeric_md_data]
    sort_select = Select(title="Sort samples by", value="input_order", options=sort_options, sizing_mode="stretch_width")

    groupby_options = {}
    groupby_options["Default"] = ["none"]
    groupby_options["Categorical Metadata"] = []
    if metadata:
        categorical_md_data = metadata.get_data(metadata_type="categorical").columns.to_list()
        if categorical_md_data:
            groupby_options["Categorical Metadata"] = [("metadata_cat|" + md, md) for md in categorical_md_data]
    groupby1_select = Select(title="1) Group samples by", value="none", options=groupby_options, sizing_mode="stretch_width")
    groupby2_select = Select(title="2) Group samples by", value="none", options=groupby_options, sizing_mode="stretch_width")

    toggle_label = CheckboxGroup(labels=["Hide sample labels"], active=[])

    help_text = """
Observation bars showing proportions of top """ + str(top_obs_bars) + """ most abundant observations.

"others" list all other observation not on top list but also observations classified in higher taxonomic ranks.

Samples can be grouped and sorted. When sorting by numeric metadata, labels will change to show the value of each sample.
"""

    return {"rank_select": rank_select,
            "groupby1_select": groupby1_select,
            "groupby2_select": groupby2_select,
            "sort_select": sort_select,
            "toggle_label": toggle_label,
            "help_button": help_button(title="Observation bars", text=help_text)}


def plot_samplebars_widgets(ranks, metadata, reference_names, control_names, decontam):
    annotbar_rank_select = Select(title="Annotate bars at rank:", value=ranks[0], options=[r for r in ranks])

    annotbar_options = {}
    annotbar_options["Default"] = ["assigned"]
    annotbar_options["References"] = [r for r in reference_names]
    annotbar_options["Controls"] = [c for c in control_names]
    if decontam:
        annotbar_options["Decontam"] = ["decontam"]
    annotbar_select = Select(title="Annotate bars by:", value="assigned", options=annotbar_options)

    y1_select = Select(title="Counts", value="#", options=["#", "%"], width=80)
    y2_select = Select(title="Observations", value="%", options=["#", "%", "log10(#)", "log10(%)"], width=80)

    sort_options = {}
    sort_options["Default"] = [("input_order", "input order"), ("#", "# counts"), ("selected_annotation", "selected annotation")]
    sort_options["Selected Rank"] = [("tax|" + r, r) for r in ranks]
    sort_options["Numeric Metadata"] = []
    if metadata:
        numeric_md_data = metadata.get_data(metadata_type="numeric").columns.to_list()
        if numeric_md_data:
            sort_options["Numeric Metadata"] = [("metadata_num|" + md, md) for md in numeric_md_data]
    sort_select = Select(title="Sort samples by", value="input_order", options=sort_options, sizing_mode="stretch_width")

    groupby_options = {}
    groupby_options["Default"] = ["none"]
    groupby_options["Categorical Metadata"] = []
    if metadata:
        categorical_md_data = metadata.get_data(metadata_type="categorical").columns.to_list()
        if categorical_md_data:
            groupby_options["Categorical Metadata"] = [("metadata_cat|" + md, md) for md in categorical_md_data]
    groupby1_select = Select(title="1) Group samples by", value="none", options=groupby_options, sizing_mode="stretch_width")
    groupby2_select = Select(title="2) Group samples by", value="none", options=groupby_options, sizing_mode="stretch_width")

    help_text = """
Bars showing total counts (left y-axis) for each sample (x-axis).

Selecting elements on the table above will plot specific observation counts for each sample (right y-axis).

Samples (x-axis) can be grouped and sorted.

Bars can be annotated by taxonomic rank. The annotation will show the proportion of counts matching the selected annotation.

Raw counts and observations (#) can be normalized (%) and/or log transformed (log10) for better visualization.
    """

    return {"y1_select": y1_select,
            "y2_select": y2_select,
            "annotbar_rank_select": annotbar_rank_select,
            "annotbar_select": annotbar_select,
            "sort_select": sort_select,
            "groupby1_select": groupby1_select,
            "groupby2_select": groupby2_select,
            "help_button": help_button(title="Sample bars", text=help_text)}


def plot_obstable(sizes, cds_p_obstable, ranks, reference_names, control_names):
    # General filter for widgets
    widgets_filter = IndexFilter()

    # Generate tables (on tabs) for each rank
    obstable_tabs = []
    # Create table with view for each rank
    for rank in ranks:
        rank_filter = GroupFilter(column_name='col|rank', group=rank)
        cds_view = CDSView(source=cds_p_obstable, filters=[rank_filter, widgets_filter])

        table_cols = []
        table_cols.append(TableColumn(field="col|name", title="Name"))
        table_cols.append(TableColumn(field="col|frequency_perc", title="Frequency", default_sort="descending", formatter=NumberFormatter(format="0.00%")))
        table_cols.append(TableColumn(field="col|counts_perc_avg", title="Avg. counts/sample", default_sort="descending", formatter=NumberFormatter(format="0.00%")))
        table_cols.append(TableColumn(field="col|total_counts", title="Total counts", default_sort="descending"))

        for ctrl_name in control_names:
            table_cols.append(TableColumn(field="col|" + ctrl_name, title="(F) " + ctrl_name, default_sort="descending", formatter=NumberFormatter(format="0.00%")))

        for ref_name in reference_names:
            table_cols.append(TableColumn(field="col|" + ref_name, title=ref_name, default_sort="descending"))

        if "col|decontam" in cds_p_obstable.data:
            table_cols.append(TableColumn(field="col|decontam", title="DECONTAM", default_sort="descending"))

        datatable = DataTable(height=sizes["overview_top_panel_height"],
                              sizing_mode="stretch_width",
                              index_position=None,
                              autosize_mode="fit_viewport",
                              #reorderable=True,
                              #selectable="checkbox",
                              frozen_columns=1,
                              columns=table_cols,
                              source=cds_p_obstable,
                              view=cds_view)

        obstable_tabs.append(Panel(child=datatable, title=rank))

    obstable = Tabs(tabs=obstable_tabs)

    return obstable, widgets_filter


def plot_obstable_widgets(sizes, dict_d_taxname, max_count_rank):
    # Filtering options
    spinner_width = sizes["overview_top_panel_width_left"] - 20
    frequency_spinner = Spinner(title="Frequency", low=0, high=100, value=0, step=1, width=spinner_width, height=50)
    counts_perc_avg_spinner = Spinner(title="Avg. counts/sample", low=0, high=100, value=0, step=0.1, width=spinner_width, height=50)
    total_counts_spinner = Spinner(title="Total counts", low=1, high=max_count_rank, step=1, value=1, width=spinner_width, height=50)
    # Create unique list of names with taxids for filtering. map to str and set to get unique
    unique_dict_d_taxname_tuples = set(zip(dict_d_taxname.keys(), map(str, dict_d_taxname.values())))
    name_multichoice = MultiChoice(title="Obs. name or id",
                                   options=list(unique_dict_d_taxname_tuples),
                                   sizing_mode="fixed",
                                   width=sizes["overview_top_panel_width_left"] - 20, height=60)

    help_text = """
Summary of observations among all samples. If taxonomy is provided, panels will show different taxonomic ranks.

Clicking on the entries will load further information of the observation in the other plots/panels.

The table contain the following fixed columns:

- **Name**: Taxonomic or given observation name
- **Frequency**: How often the observation is occurring among all samples
- **Avg. counts/sample**: Averge percentage of the observation among all samples
- **Total counts**: Sum of counts of this observation in all samples
- **(F) Controls**: (F)requency for controls: how often the observation is occurring in the given control samples
- **References**:  How many times the observation was reported in the references
- **DECONTAM**: Final contamination output from DECONTAM method

Widgets can filter entries of the table. "Obs. name or id" filters the lineage of the entries, if taxonomy is provided. With that is possible to, for example, filter a certain genus and the table will show only children species.
    """

    return {"frequency_spinner": frequency_spinner,
            "counts_perc_avg_spinner": counts_perc_avg_spinner,
            "total_counts_spinner": total_counts_spinner,
            "name_multichoice": name_multichoice,
            "help_button": help_button(title="Observation table", text=help_text, align="start")}


def plot_infopanel():
    return TextAreaInput(value="Click on the table items to load more information",
                         sizing_mode="stretch_both",
                         disabled=False)


def plot_decontam(sizes, cds_p_decontam, cds_p_decontam_lines, min_obs_perc):
    decontam_fig = figure(x_axis_type="log",
                          y_axis_type="log",
                          height=sizes["overview_top_panel_height"] - 50,
                          width=sizes["overview_top_panel_width_right"],
                          sizing_mode="stretch_width",
                          tools="save")

    palette = make_color_palette(2)
    factors = list(set(cds_p_decontam.data["controls"]))

    # Add legend on top
    decontam_fig.add_layout(Legend(), 'above')
    decontam_fig.circle(x="concentration", y="counts",
                        source=cds_p_decontam,
                        color=factor_cmap('controls', palette=palette, factors=factors),
                        legend_group="controls",
                        size=3)

    # If there are controls, show legend
    if len(factors) > 1:
        decontam_fig.legend.margin = 0
        decontam_fig.legend.border_line_width = 0
        decontam_fig.legend.spacing = 0
        decontam_fig.legend.padding = 0
        decontam_fig.legend.orientation = "horizontal"
        decontam_fig.legend.location = "bottom_right"
    else:
        decontam_fig.legend.visible = False

    decontam_fig.line(x="x",
                      y="y_cont",
                      source=cds_p_decontam_lines,
                      color="red")
    decontam_fig.line(x="x",
                      y="y_noncont",
                      source=cds_p_decontam_lines,
                      color="black",
                      line_dash="dashed")

    decontam_fig.xaxis.axis_label = 'DNA Concentration/Total counts'
    decontam_fig.yaxis.axis_label = 'obs. counts'
    decontam_fig.y_range.start = min_obs_perc
    decontam_fig.y_range.end = 1

    return decontam_fig


def plot_decontam_widgets():
    pvalue_text = Paragraph(text="P-value")
    pvalue_input = TextInput(value="", width=180, align='end')

    help_text = """
Plot to verify the DECONTAM [1] output. Proportion of counts (y-axis) against DNA Concentration (if provided) or Total number of counts (x-axis), both in log10 scale. If provided, controls samples are displayed in a different color.

A indication of contamination is when counts are inversely proportional to DNA concentration. The red and black dotted lines are the expected models for contamination and non-contamination, respectively. A good indication for contamination is when the dots (excluding control samples) "fit" the red line model.

More details can be found in the [DECONTAM Introduction guide](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)

[1] Davis, N. M. et al. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome (2018) 10.1186/s40168-018-0605-2.
"""
    return {"pvalue_text": pvalue_text,
            "pvalue_input": pvalue_input, 
            "help_button": help_button(title="DECONTAM", text=help_text, align="start")}


def plot_references(sizes, table, cds_p_references, dict_d_taxname):
    references_fig = figure(x_range=table.ranks(),
                            height=sizes["overview_top_panel_height"] - 50,
                            width=sizes["overview_top_panel_width_right"],
                            tools="save,reset")

    # Need to pass dict_d_taxname inside a one column data
    taxid_name_custom = CustomJSHover(
        args=dict(dict_d_taxname=ColumnDataSource(dict(dict_d_taxname=[dict_d_taxname]))),
        code="return dict_d_taxname.data.dict_d_taxname[0][value]; // value holds the @taxid"
    )
    # Add custom tooltip for heatmap (taxid->name)
    references_fig.add_tools(HoverTool(
        tooltips=[
            ('Observation', '@obs{custom}'),
            ('# reported (directly)', '@direct'),
            ('as child', '@child'),
            ('as parent', '@parent'),
        ],
        mode="mouse",
        point_policy="follow_mouse",
        formatters={"@obs": taxid_name_custom}
    ))

    references_filter = IndexFilter(indices=[])
    cds_view_references = CDSView(source=cds_p_references, filters=[references_filter])

    references_fig.add_layout(Legend(), 'above')

    fixed_bar_options = ["direct", "child", "parent"]
    palette = ["red", "orange", "black"]
    references_fig.vbar_stack(fixed_bar_options,
                              x="rank",
                              width=1,
                              source=cds_p_references,
                              view=cds_view_references,
                              color=palette,
                              line_color=None,  # to avoid printing small border for zeros
                              fill_alpha=[1, 0.3, 0.3],
                              legend_label=fixed_bar_options)

    references_fig.y_range.start = 0

    references_fig.xaxis.major_label_orientation = "vertical"
    references_fig.xgrid.grid_line_color = None
    references_fig.xaxis.minor_tick_line_color = None
    references_fig.yaxis.minor_tick_line_color = None
    references_fig.xaxis.major_tick_line_color = None
    references_fig.yaxis.major_tick_line_color = None
    references_fig.yaxis.axis_label = "# reported"

    references_fig.legend.margin = 0
    references_fig.legend.border_line_width = 0
    references_fig.legend.spacing = 0
    references_fig.legend.padding = 0
    references_fig.legend.orientation = "horizontal"
    references_fig.legend.location = "bottom_right"

    return references_fig, references_filter


def plot_references_widgets(sizes, references):
    references_select = Select(value=list(references.keys())[0], width=sizes["overview_top_panel_width_right"] - 70, options=list(references.keys()))
    help_text = """
Plot of number of occurences of provided references for each observation and its lineage.

**direct** counts represent direct matches with reference identifiers

**child** counts accounts for the number of times a related parent on the lineage of the selected observation node was reported among references

**parent** counts accounts for the number of times related children (not necessarily reported) on the lineage of the selected observation node was reported among references
"""

    return {"references_select": references_select,
            "help_button": help_button(title="References", text=help_text, align="start")}


def plot_mgnify(sizes, cds_p_mgnify):
    mgnify_fig = figure(height=sizes["overview_top_panel_height"] - 50,
                        width=sizes["overview_top_panel_width_right"],
                        tools="save,wheel_zoom,reset")

    mgnify_filter = IndexFilter(indices=[])
    cds_view_mgnify = CDSView(source=cds_p_mgnify, filters=[mgnify_filter])

    # make mapping of lineages for each level
    level_lineage = {}
    for level, lineage in zip(cds_p_mgnify.data["level"], cds_p_mgnify.data["lineage"]):
        if level not in level_lineage:
            level_lineage[level] = set()
        level_lineage[level].add(lineage)

    # build unique pallete for each level of biome
    factors = []
    palette = []
    for level, lineages in level_lineage.items():
        factors.extend(lineages)
        palette.extend(make_color_palette(len(lineages)))

    # Add custom tooltip to show percentage (based on angle)
    mgnify_fig.add_tools(HoverTool(
        tooltips=[("Biome", "@lineage"),
                  ("Studies", "@count"),
                  ("Percentage", "@angle{custom}%")],
        mode="mouse",
        point_policy="follow_mouse",
        formatters={"@angle": CustomJSHover(code="return ((value/6.2831853071795)*100).toFixed(2);")}
    ))

    #mgnify_fig.text(0, 1, text=["No data"], text_baseline="middle", text_align="center")
    mgnify_fig.wedge(x=0, y=1, radius=0.5,
                     start_angle=cumsum('angle', include_zero=True),
                     end_angle=cumsum('angle'),
                     #line_color="white",
                     line_width=0,
                     fill_color=factor_cmap('lineage', palette=palette, factors=factors),
                     source=cds_p_mgnify,
                     view=cds_view_mgnify)

    mgnify_fig.axis.axis_label = None
    mgnify_fig.axis.visible = False
    mgnify_fig.grid.grid_line_color = None

    return mgnify_fig, mgnify_filter


def plot_mgnify_widgets():
    biome_spinner = Spinner(title="Biome level", low=1, high=5, value=1, step=1, width=100, height=50, orientation="horizontal")
    help_text = """
Pie chart with the number of occurrences of the selected taxa in the table by environment (biome) in other microbiome studies analyzed and publicly available at the [MGNify](https://www.ebi.ac.uk/metagenomics) [1] resource.

The biomes are hierarchically classified in 5 different levels, from general (1) to specific (5). For example: 1) Host-associated > 2) Human > 3) Digestive system > 4) Large intestine > 5) Fecal

[1] Mitchell, A. L. et al. MGnify: the microbiome analysis resource in 2020. Nucleic Acids Res (2020) 10.1093/nar/gkz1035.
    """

    return {"biome_spinner": biome_spinner,
            "help_button": help_button(title="MGNify", text=help_text)}


def plot_heatmap(table, cds_p_heatmap, tools_heatmap, transformation, dict_d_taxname):
    heatmap = figure(x_range=table.observations(table.ranks()[0]).to_list(),
                     y_range=table.samples.to_list(),
                     height=500,
                     sizing_mode="stretch_width",
                     x_axis_location="above",
                     tools=tools_heatmap,
                     toolbar_location='right',
                     tooltips="")

    # Need to pass dict_d_taxname inside a one column data
    taxid_name_custom = CustomJSHover(
        args=dict(dict_d_taxname=ColumnDataSource(dict(dict_d_taxname=[dict_d_taxname]))),
        code="return dict_d_taxname.data.dict_d_taxname[0][value]; // value holds the @taxid"
    )
    # Add custom tooltip for heatmap (taxid->name)
    heatmap.add_tools(HoverTool(
        tooltips=[
            ('Sample', '@index'),
            ('Observation', '@obs{custom}'),
            ('Original value', '@ov'),
            ('Transformed value (' + transformation + ')', '@tv')
        ],
        formatters={"@obs": taxid_name_custom}
    ))

    color_palette = Magma256[::-1]
    color_mapper = LinearColorMapper(palette=color_palette)
    color_mapper.low = min(cds_p_heatmap.data["tv"])
    color_mapper.high = max(cds_p_heatmap.data["tv"])

    # Convert taxid ticks to taxa names on client-side
    heatmap.xaxis.formatter = FuncTickFormatter(args=dict(dict_d_taxname=dict_d_taxname), code='''
        return dict_d_taxname[tick];
    ''')

    heatmap.rect(x="obs", y="index", width=1, height=1,
                 source=cds_p_heatmap,
                 fill_color={'field': 'tv', 'transform': color_mapper},
                 line_color=None)

    color_bar = ColorBar(color_mapper=color_mapper, label_standoff=12, border_line_color=None, location="center", orientation="horizontal")
    heatmap.add_layout(color_bar, 'above')

    # Convert taxid ticks to taxa names on client-side
    heatmap.xaxis.formatter = FuncTickFormatter(args=dict(dict_d_taxname=dict_d_taxname), code='''
        return dict_d_taxname[tick];
    ''')

    heatmap.xaxis.major_label_orientation = "vertical"
    heatmap.xgrid.grid_line_color = None
    heatmap.ygrid.grid_line_color = None
    heatmap.xaxis.minor_tick_line_color = None
    heatmap.yaxis.minor_tick_line_color = None
    heatmap.xaxis.major_tick_line_color = None
    heatmap.yaxis.major_tick_line_color = None
    heatmap.xaxis.major_label_text_font_size = "0px"
    heatmap.yaxis.major_label_text_font_size = "0px"
    heatmap.xaxis.axis_label = "observations"
    heatmap.yaxis.axis_label = "samples"

    return heatmap


def plot_heatmap_widgets(ranks, linkage_methods, linkage_metrics, reference_names, controls_names, metadata, decontam):

    rank_select = Select(title="Taxonomic rank:", value=ranks[0], options=ranks)

    x_sort_options = {}
    x_sort_options["Clustering Metric"] = [("metric|" + lm, lm) for lm in linkage_metrics]
    x_sort_options["Default order"] = [("none", "none"), ("counts", "counts"), ("observations", "observations")]
    x_sort_options["Sort by References"] = [("annot|" + r, r) for r in reference_names]
    if controls_names:
        x_sort_options["Sort by Controls"] = [("annot|" + c, c) for c in controls_names]
    if decontam:
        x_sort_options["Sort by DECONTAM"] = [("annot|decontam", "decontam")]

    y_sort_options = {}
    y_sort_options["Clustering Metric"] = [("metric|" + lm, lm) for lm in linkage_metrics]
    y_sort_options["Default order"] = [("none", "none"), ("counts", "counts"), ("samples", "samples")]

    if metadata:
        numeric_md_data = metadata.get_data(metadata_type="numeric").columns.to_list()
        if numeric_md_data:
            y_sort_options["Sort by Numeric Metadata"] = [("metadata_num|" + md, md) for md in numeric_md_data]
        categorical_md_data = metadata.get_data(metadata_type="categorical").columns.to_list()
        if categorical_md_data:
            y_sort_options["Sort by Categorical Metadata"] = [("metadata_cat|" + md, md) for md in categorical_md_data]


    x_sort_select = Select(title="Observation cluster/sort:", value="none", options=x_sort_options)
    x_method_select = Select(title="Observation clustering method:", value=linkage_methods[0], options=linkage_methods, disabled=True)
    #x_group_select = Select(title="Observation group by:", value="none", options={"Default": ["none"], "Taxonomic ranks": ranks})
    
    y_sort_select = Select(title="Sample cluster/sort:", value="none", options=y_sort_options)
    y_method_select = Select(title="Sample clustering method:", value=linkage_methods[0], options=linkage_methods, disabled=True)

    toggle_label_text = Div(text="show/hide labels")
    toggle_label_heatmap = CheckboxGroup(labels=["Observations", "Samples"], active=[])

    help_text = """
The heatmap shows [transformed] values from the input table (color bar on top). If taxonomy is provided, one heatmap for each taxonomic rank is generated.

Values on axis can be independently clustered or sorted. Hierarchical clustering is done with [scipy linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html). If clustering is selected, dendrograms will be plotted in the panels around the heatmap.

The right-most panel will show metadata values related to each sample (y-axis), if provided. Colors are automatically generated for categorical (distinct colors) and numeric (sequential colors) fields. Multiple metadata fields can be select with Ctrl + click in the metadata selection widget.

The bottom-most panel shows annotations for each observation/taxa values (x-axis).

The metadata and annotation plots are automatically sorted to reflect the clustering/sort of the heatmap.
"""

    return {"rank_select": rank_select,
            "x_method_select": x_method_select,
            "x_sort_select": x_sort_select,
            #"x_group_select": x_group_select,
            "y_method_select": y_method_select,
            "y_sort_select": y_sort_select,
            "toggle_label_heatmap": toggle_label_heatmap,
            "toggle_label_text": toggle_label_text,
            "help_button": help_button(title="Heatmap/Clustering", text=help_text)}


def plot_dendrogram(heatmap, tools_heatmap, cds_p_dendro_x, cds_p_dendro_y):
    
    dendrox_fig = figure(x_range=heatmap.x_range,
                         tools="save",
                         height=80,
                         sizing_mode="stretch_width",
                         tooltips=[("y", "$y{(0.00)}"), ("c", "$swatch:c")])
    dendroy_fig = figure(y_range=heatmap.y_range,
                         tools="save",
                         width=80,
                         height=heatmap.height,
                         tooltips=[("x", "$x{(0.00)}"), ("c", "$swatch:c")])

    dendroy_fig.multi_line(xs="x", ys="y", color="c",
                           source=cds_p_dendro_y)

    dendrox_fig.multi_line(xs="x", ys="y", color="c",
                           source=cds_p_dendro_x)

    dendroy_fig.xaxis.major_tick_line_color = None
    dendroy_fig.xaxis.minor_tick_line_color = None
    dendroy_fig.yaxis.major_tick_line_color = None
    dendroy_fig.yaxis.minor_tick_line_color = None
    dendroy_fig.xaxis.major_label_text_font_size = '0pt'  # preferred method for removing tick labels
    dendroy_fig.yaxis.major_label_text_font_size = '0pt'  # preferred method for removing tick labels
    dendroy_fig.xaxis.axis_line_color = None
    dendroy_fig.yaxis.axis_line_color = None
    dendroy_fig.xgrid.grid_line_color = None
    dendroy_fig.ygrid.grid_line_color = None
    dendrox_fig.xaxis.major_tick_line_color = None
    dendrox_fig.xaxis.minor_tick_line_color = None
    dendrox_fig.yaxis.major_tick_line_color = None
    dendrox_fig.yaxis.minor_tick_line_color = None
    dendrox_fig.xaxis.major_label_text_font_size = "0pt"
    dendrox_fig.yaxis.major_label_text_font_size = "0pt"
    dendrox_fig.xaxis.axis_line_color = None
    dendrox_fig.yaxis.axis_line_color = None
    dendrox_fig.xgrid.grid_line_color = None
    dendrox_fig.ygrid.grid_line_color = None

    return dendrox_fig, dendroy_fig


def plot_metadata(heatmap, tools_heatmap, metadata, cds_d_metadata, cds_p_metadata):
    # Number of cols added to the plot cds (user input)
    ncols = len(cds_p_metadata.data) - 1

    metadata_fig = figure(x_range=["md" + str(c) for c in range(ncols)],
                          y_range=heatmap.y_range,
                          tools=tools_heatmap,
                          x_axis_location="above",
                          width=80, height=heatmap.height,
                          tooltips="")
    metadata_fields = metadata.get_col_headers().to_list()

    # A bit of a hack of the "proper" user of a colormap, since I need mixed types (categorical and numeric)
    # and different scales for multiple columns
    # Build colormap for all unique entries in the metadata for each metadata field
    # (metadata header, str(metadata value)) -> color
    # Add them to a CategoricalColorMapper which will be applied for the whole plot
    # Use different palettes for numeric types, but convert to string to be treated as a category
    # Need to be careful with int and float, since the value of str(0.0) 
    # will not match the 0 in the javascript data conversion, therefore use the numeric to calculate palette
    # but make get_formatted_unique_values on the dictionary
    factors = []
    palette = []
    for i, md_header in enumerate(metadata_fields):
        unique_values = metadata.get_unique_values(md_header)
        if unique_values:
            n = len(unique_values)
            if metadata.get_type(md_header) == "numeric":
                unique_palette = make_color_palette(n, linear=True)
            else:
                unique_palette = make_color_palette(n)
            assert len(unique_palette) == n, 'Wrong number of colors on palette'
            palette.extend(unique_palette)
            factors.extend([(md_header, md_value) for md_value in metadata.get_formatted_unique_values(md_header)])
    metadata_colormap = CategoricalColorMapper(palette=palette, factors=factors)

    # Custom tooltip to show metadata field and value
    md_custom = CustomJSHover(code='return value[0] ? "(" + value[0] + ") " + value[1] : "";')
    tooltips = [('Sample', '@index')]
    formatters = {}
    for i in range(ncols):
        tooltips.append(("md" + str(i), "@md" + str(i) + "{custom}"))
        formatters["@md" + str(i)] = md_custom
    metadata_fig.add_tools(HoverTool(tooltips=tooltips, formatters=formatters))

    for i in range(ncols):
        metadata_fig.rect(x=i + 0.5, y="index",
                          width=1, height=1,
                          source=cds_p_metadata,
                          fill_color={'field': "md" + str(i), 'transform': metadata_colormap},
                          line_color=None)

    metadata_fig.xaxis.axis_label = "metadata"
    metadata_fig.xaxis.major_label_orientation = "vertical"
    metadata_fig.xaxis.major_label_text_font_size = "11px"
    metadata_fig.xaxis.minor_tick_line_color = None

    metadata_fig.yaxis.major_tick_line_color = None
    metadata_fig.yaxis.minor_tick_line_color = None
    metadata_fig.yaxis.major_label_text_font_size = '0pt'
    metadata_fig.yaxis.axis_line_color = None
    metadata_fig.ygrid.grid_line_color = None

    metadata_multiselect = MultiSelect(title="Metadata to show ('md#', max. " + str(ncols) + "):", value=[metadata_fields[c] for c in range(ncols)], options=metadata_fields)

    return metadata_fig, {"metadata_multiselect": metadata_multiselect}


def plot_annotations(heatmap, tools_heatmap, cds_p_annotations, dict_d_taxname):
    # Sorted list of annotations available in the CDS
    rows = list(sorted(set(cds_p_annotations.data["annot"])))[::-1]

    annot_fig = figure(x_range=heatmap.x_range,
                       y_range=rows,
                       tools=tools_heatmap,
                       height=120, sizing_mode="stretch_width",
                       tooltips="")

    # Need to pass dict_d_taxname inside a one column data
    taxid_name_custom = CustomJSHover(
        args=dict(dict_d_taxname=ColumnDataSource(dict(dict_d_taxname=[dict_d_taxname]))),
        code="return dict_d_taxname.data.dict_d_taxname[0][value]; // value holds the @taxid"
    )
    # Add custom tooltip for heatmap (taxid->name)
    annot_fig.add_tools(HoverTool(
        tooltips=[('Annotation', '@annot'), ('Observation', '@index{custom}')],
        formatters={"@index": taxid_name_custom}
    ))

    annot_fig.rect(x="index", y="annot",
                   width=1, height=1,
                   source=cds_p_annotations,
                   fill_color="black",
                   line_color=None)

    annot_fig.yaxis.axis_label = "annotations"
    annot_fig.xgrid.grid_line_color = None
    annot_fig.ygrid.grid_line_color = None
    annot_fig.xaxis.minor_tick_line_color = None
    annot_fig.yaxis.minor_tick_line_color = None
    annot_fig.xaxis.major_tick_line_color = None
    annot_fig.xaxis.major_label_text_font_size = "0px"

    return annot_fig


def plot_correlation(cds_p_correlation, ranks, dict_d_taxname):
    taxids = set()
    for i, rank in enumerate(cds_p_correlation.data["rank"]):
        if rank == ranks[0]:
            taxids.add(cds_p_correlation.data["index"][i])

    taxids = sorted(taxids)
    corr_fig = figure(x_range=taxids,
                      y_range=list(reversed(taxids)),
                      tools="hover,save,reset,crosshair,tap,box_zoom",
                      tooltips="",
                      sizing_mode="scale_height")

    # Need to pass dict_d_taxname inside a one column data
    taxid_name_custom = CustomJSHover(
        args=dict(dict_d_taxname=ColumnDataSource(dict(dict_d_taxname=[dict_d_taxname]))),
        code="return dict_d_taxname.data.dict_d_taxname[0][value]; // value holds the @taxid"
    )
    # Add custom tooltip for heatmap (taxid->name)
    corr_fig.add_tools(HoverTool(
        tooltips=[
            ('x', '@taxid{custom}'),
            ('y', '@index{custom}'),
            ('Correlation coefficient (rho)', '@rho'),
            ('Raw p-value', '@pval'),
            ('Corrected p-value', '@pval_corr')
        ],
        formatters={"@taxid": taxid_name_custom, "@index": taxid_name_custom}
    ))

    color_palette = Blues[9][:-1] + ("#ffffff",) + Reds[9][:-1][::-1]
    color_mapper = LinearColorMapper(palette=color_palette, low=-1, high=1)

    rho_filter = IndexFilter()
    pval_filter = IndexFilter()
    cds_view_correlation = CDSView(source=cds_p_correlation, filters=[rho_filter, pval_filter])
    corr_fig.rect(x="taxid", y="index",
                  width=1, height=1,
                  source=cds_p_correlation,
                  view=cds_view_correlation,
                  fill_color={'field': 'rho', 'transform': color_mapper},
                  line_color=None)

    color_bar = ColorBar(color_mapper=color_mapper,
                         label_standoff=12,
                         ticker=AdaptiveTicker(),
                         border_line_color=None,
                         location="center",
                         orientation="horizontal")
    corr_fig.add_layout(color_bar, 'above')

    # Convert taxid ticks to taxa names on client-side
    corr_fig.xaxis.formatter = FuncTickFormatter(args=dict(dict_d_taxname=dict_d_taxname), code='''
        return dict_d_taxname[tick];
    ''')
    corr_fig.yaxis.formatter = FuncTickFormatter(args=dict(dict_d_taxname=dict_d_taxname), code='''
        return dict_d_taxname[tick];
    ''')
    corr_fig.xgrid.grid_line_color = None
    corr_fig.ygrid.grid_line_color = None
    corr_fig.xaxis.minor_tick_line_color = None
    corr_fig.yaxis.minor_tick_line_color = None
    corr_fig.xaxis.major_label_orientation = "vertical"

    return corr_fig, rho_filter, pval_filter


def plot_correlation_widgets(ranks, top_obs_corr):
    rank_select = Select(title="Taxonomic rank:", value=ranks[0], options=ranks)
    neg_slider = RangeSlider(start=-1, end=0, value=(-1, 0), step=.01, title="Negative correlation")
    pos_slider = RangeSlider(start=0, end=1, value=(0, 1), step=.01, title="Positive correlation")
    pval_spinner = Spinner(title="Corrected P-value", low=0, high=1, step=0.01, value=1, width=100, height=50)

    help_text = """
Spearman correlation coefficient with associated and corrected (Benjamini/Hochberg) p-values between the top """ + str(top_obs_corr) + """ most abundant observations. [spearmanr](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html) function from scipy is used.

- Negative correlation values [-1 .. 0] are displayed in blue.
- Positive correlation values [0 .. 1] are displayed in red.

Only half matrix is displayed, since the values are symmetric. Widgets can control level of correlation and corrected p-value to show.
"""
    return {"rank_select": rank_select,
            "neg_slider": neg_slider,
            "pos_slider": pos_slider,
            "pval_spinner": pval_spinner,
            "help_button": help_button(title="Correlation", text=help_text)}


def help_button(title: str="", text: str="", align: str="end"):
    hb = Button(width=32, height=32, label="?", align=align, button_type="warning")

    # Convert markdown to html
    html_text = markdown.markdown(text)
    # Open links on new page
    html_text = html_text.replace("<a ", "<a target=\"_blank\" ")
    # Put code in one line
    html_text = html_text.replace("\n", "")

    hb.js_on_click(CustomJS(code="pop.open('" + title + "', '" + html_text + "');"))
    return hb


def make_color_palette(n_colors, linear: bool=False, palette: dict=None):
    if isinstance(palette, dict) and n_colors <= max(palette.keys()):
        # Special case for 1 and 2 (not in palettes)
        palette = palette[3 if n_colors < 3 else n_colors]

    if linear or n_colors > 20:
        if not palette:
            palette = Turbo256
        if n_colors <= 256:
            return linear_palette(palette, n_colors)
        else:
            # Repeat colors
            return [palette[int(i * 256.0 / n_colors)] for i in range(n_colors)]
    else:
        # Select color palette based on number of requested colors
        # Return the closest palette with most distinc set of colors
        if not palette:
            if n_colors <= 8:
                palette = Colorblind[8]
            elif n_colors <= 10:
                palette = Category10[10]
            elif n_colors <= 20:
                palette = Category20[20]
            else:
                palette = Turbo256

        return palette[:n_colors]
