from bokeh.layouts import column, row, gridplot
from bokeh.models import Spacer, Tabs, Panel, Div
import base64


def make_layout(ele, sizes, version, logo_path, title):

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
                                            row(ele["decontam"]["wid"]["pvalue_text"],
                                                ele["decontam"]["wid"]["pvalue_input"],
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
                            ele["samplebars"]["wid"]["help_button"]))

    row_heatmap = gridplot([[ele["heatmap"]["fig"], ele["dendroy"]["fig"], ele["metadata"]["fig"]],
                           [ele["dendrox"]["fig"]],
                           [ele["annotations"]["fig"], ele["heatmap"]["wid"]["help_button"]]],
                           toolbar_location='right',
                           merge_tools=True)

    row_heatmap_widgets = row(column(ele["heatmap"]["wid"]["rank_select"],
                                     row(ele["heatmap"]["wid"]["toggle_label_text"],
                                         ele["heatmap"]["wid"]["toggle_label_heatmap"]),
                                     sizing_mode="stretch_height",
                                     width=300),
                              column(row(ele["heatmap"]["wid"]["x_sort_select"],
                                         ele["heatmap"]["wid"]["y_sort_select"],
                                         sizing_mode="stretch_width"),
                                     row(ele["heatmap"]["wid"]["x_method_select"],
                                         ele["heatmap"]["wid"]["y_method_select"],
                                         sizing_mode="stretch_width"),
                                     #row(ele["heatmap"]["wid"]["x_group_select"],
                                     #    Spacer(),
                                     #    sizing_mode="stretch_width"),
                                     sizing_mode="stretch_width"),
                              column(ele["metadata"]["wid"]["metadata_multiselect"],
                                     sizing_mode="stretch_height",
                                     width=300))

    row_correlation = row(column(ele["correlation"]["wid"]["rank_select"],
                                 ele["correlation"]["wid"]["neg_slider"],
                                 ele["correlation"]["wid"]["pos_slider"],
                                 ele["correlation"]["wid"]["help_button"]),
                          ele["correlation"]["fig"])

    row_obsbars = column(row(ele["obsbars"]["fig"]),
                         ele["obsbars"]["wid"]["toggle_label"],
                         row(ele["obsbars"]["wid"]["rank_select"],
                             ele["obsbars"]["wid"]["groupby1_select"],
                             ele["obsbars"]["wid"]["groupby2_select"],
                             ele["obsbars"]["wid"]["sort_select"],
                             ele["obsbars"]["wid"]["help_button"]))

    main_panels = []
    main_panels.append(Panel(child=column(row_obstable, row_barpot, sizing_mode="stretch_width"), title="Overview"))
    main_panels.append(Panel(child=column(row_heatmap, row_heatmap_widgets, sizing_mode="stretch_width"), title="Heatmap"))
    main_panels.append(Panel(child=column(row_correlation, sizing_mode="stretch_width"), title="Correlation"))
    main_panels.append(Panel(child=column(row_obsbars, sizing_mode="stretch_width"), title="Bars"))
    main_tab = Tabs(tabs=main_panels)

    logo_base64 = base64.b64encode(open(logo_path, 'rb').read())  # encode to base64
    logo_base64 = logo_base64.decode()    # convert to string
    logo_div = Div(text='<img src="data:image/png;base64,' + logo_base64 + '">' + '<a target="_blank" style="color: black" href="https://github.com/pirovc/grimer">v' + version + '</a>', width=300, height=40, sizing_mode="fixed")
    if title:
        title_div = Div(text='<h2>' + title + '</h2>', height=40, sizing_mode="stretch_width")
    else:
        title_div = Spacer()
    final = column([row(logo_div, title_div), main_tab], sizing_mode="stretch_width")

    return final