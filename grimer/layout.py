from bokeh.layouts import column, row, gridplot
from bokeh.models import Spacer, Tabs, Panel, Div
from grimer.func import print_log
import base64


def make_layout(ele, sizes, version, logo_path, title, output_plots):

    main_panels = {}
    if "overview" in output_plots:
        filterwidgets = column(ele["obstable"]["wid"]["frequency_spinner"],
                               ele["obstable"]["wid"]["counts_perc_avg_spinner"],
                               ele["obstable"]["wid"]["total_counts_spinner"],
                               ele["obstable"]["wid"]["name_multichoice"],
                               ele["obstable"]["wid"]["help_button"])
        filterwidgetstabs = Tabs(tabs=[Panel(child=filterwidgets, title="Filter")],
                                 sizing_mode="fixed",
                                 height=sizes["overview_top_panel_height"] + 20,
                                 width=sizes["overview_top_panel_width_left"])
        info_tabs = [Panel(child=ele["infopanel"]["textarea"], title="Info")]
        if ele["references"]["fig"]:
            info_tabs.append(Panel(child=column(ele["references"]["fig"],
                                                row(ele["references"]["wid"]["references_select"],
                                                    ele["references"]["wid"]["help_button"])
                                                ), title="References"))
        if ele["mgnify"]["fig"]:
            info_tabs.append(Panel(child=column(ele["mgnify"]["fig"],
                                                row(ele["mgnify"]["wid"]["biome_spinner"],
                                                    ele["mgnify"]["wid"]["help_button"])
                                                ), title="MGNify"))
        if ele["decontam"]["fig"]:
            info_tabs.append(Panel(child=column(ele["decontam"]["fig"],
                                                row(ele["decontam"]["wid"]["pscore_text"],
                                                    ele["decontam"]["wid"]["pscore_input"],
                                                    ele["decontam"]["wid"]["help_button"])
                                                ), title="DECONTAM"))
        infotabs = Tabs(tabs=info_tabs,
                        sizing_mode="fixed",
                        height=sizes["overview_top_panel_height"] + 20,
                        width=sizes["overview_top_panel_width_right"])
        row_obstable = row(filterwidgetstabs,
                           ele["obstable"]["fig"],
                           infotabs,
                           sizing_mode="stretch_width")
        row_barpot = column(row(ele["samplebars"]["fig"]),
                            row(ele["samplebars"]["wid"]["y1_select"],
                                ele["samplebars"]["wid"]["annotbar_rank_select"],
                                ele["samplebars"]["wid"]["annotbar_select"],
                                ele["samplebars"]["wid"]["groupby1_select"],
                                ele["samplebars"]["wid"]["groupby2_select"],
                                ele["samplebars"]["wid"]["sort_select"],
                                ele["samplebars"]["wid"]["y2_select"],
                                ele["samplebars"]["wid"]["help_button"]),
                            ele["samplebars"]["wid"]["toggle_label"])
        main_panels["overview"] = Panel(child=column(row_obstable, row_barpot, sizing_mode="stretch_width"), title="Overview")

    if "samples" in output_plots:
        selectwidgets = column(ele["sampletable"]["wid"]["total_counts_spinner"],
                               ele["sampletable"]["wid"]["assigned_spinner"],
                               ele["sampletable"]["wid"]["metadata_multichoice"],
                               ele["sampletable"]["wid"]["help_button"])
        selectwidgetstabs = Tabs(tabs=[Panel(child=selectwidgets, title="Select")],
                                 sizing_mode="fixed",
                                 height=sizes["overview_top_panel_height"] + 20,
                                 width=sizes["overview_top_panel_width_left"])
        row_sampletable = row(selectwidgetstabs,
                              ele["sampletable"]["fig"],
                              sizing_mode="stretch_width")
        row_obsbars = column(row(ele["obsbars"]["fig"]),
                             row(ele["obsbars"]["wid"]["rank_select"],
                                 ele["obsbars"]["wid"]["groupby1_select"],
                                 ele["obsbars"]["wid"]["groupby2_select"],
                                 ele["obsbars"]["wid"]["sort_select"],
                                 ele["obsbars"]["wid"]["help_button"]),
                             ele["obsbars"]["wid"]["toggle_label"])
        main_panels["samples"] = Panel(child=column(row_sampletable, row_obsbars, sizing_mode="stretch_width"), title="Samples")

    if "heatmap" in output_plots:
        row_heatmap = gridplot([[ele["heatmap"]["fig"], ele["dendroy"]["fig"], ele["metadata"]["fig"]],
                               [ele["dendrox"]["fig"]],
                               [ele["annotations"]["fig"], None, ele["heatmap"]["wid"]["help_button"]]],
                               toolbar_location='right',
                               merge_tools=True)
        row_heatmap_widgets = row(column(ele["heatmap"]["wid"]["rank_select"],
                                         ele["heatmap"]["wid"]["toggle_label"],
                                         width=300),
                                  row(column(ele["heatmap"]["wid"]["x_groupby_select"],
                                             ele["heatmap"]["wid"]["x_sort_select"]),
                                      column(ele["heatmap"]["wid"]["y_groupby_select"],
                                             ele["heatmap"]["wid"]["y_sort_select"]),
                                      sizing_mode="stretch_width"),
                                  column(ele["metadata"]["wid"]["metadata_multiselect"],
                                         ele["metadata"]["wid"]["toggle_legend"],
                                         sizing_mode="stretch_height",
                                         width=300))
        main_panels["heatmap"] = Panel(child=column(row_heatmap, row_heatmap_widgets, sizing_mode="stretch_width"), title="Heatmap")

    if "correlation" in output_plots:
        row_correlation = row(column(ele["correlation"]["wid"]["rank_select"],
                                     ele["correlation"]["wid"]["neg_slider"],
                                     ele["correlation"]["wid"]["pos_slider"],
                                     ele["correlation"]["wid"]["toggle_label"],
                                     ele["correlation"]["wid"]["help_button"]),
                              ele["correlation"]["fig"])
        main_panels["correlation"] = Panel(child=column(row_correlation, sizing_mode="stretch_width"), title="Correlation")

    if not main_panels:
        print_log("No valid plots to output")
        return None
    else:
        # Add plots in user chosen order
        tabs = [main_panels[p] for p in output_plots]

    main_tab = Tabs(tabs=tabs)
    logo_base64 = base64.b64encode(open(logo_path, 'rb').read())  # encode to base64
    logo_base64 = logo_base64.decode()    # convert to string
    logo_div = Div(text='<img src="data:image/png;base64,' + logo_base64 + '">' + '<a target="_blank" style="color: black" href="https://github.com/pirovc/grimer">v' + version + '</a>', width=300, height=40, sizing_mode="fixed")
    if title:
        title_div = Div(text='<h2>' + title + '</h2>', height=40, sizing_mode="stretch_width")
    else:
        title_div = Spacer()
    final = column([row(logo_div, title_div), main_tab], sizing_mode="stretch_width")

    return final
