from bokeh.models import CustomJS


def link_obstable_samplebars(ele,
                             cds_p_obstable,
                             cds_p_samplebars,
                             cds_d_samples,
                             cds_d_sampleobs,
                             cds_d_metadata,
                             cds_p_decontam,
                             cds_p_decontam_models,
                             cds_d_decontam,
                             cds_p_references,
                             active_ranks,
                             min_obs_perc,
                             max_total_count,
                             cds_p_mgnify,
                             dict_d_refs):

    bar_select_callback = CustomJS(
        args=dict(y1_select=ele["samplebars"]["wid"]["y1_select"],
                  annotbar_select=ele["samplebars"]["wid"]["annotbar_select"],
                  annotbar_rank_select=ele["samplebars"]["wid"]["annotbar_rank_select"],
                  cds_p_samplebars=cds_p_samplebars,
                  cds_d_samples=cds_d_samples,
                  y_range=ele["samplebars"]["fig"].y_range,
                  max_total_count=max_total_count),
        code='''
        console.log("bar_select_callback");

        const pdata = cds_p_samplebars.data;
        const ddata = cds_d_samples.data;
        const total = cds_d_samples.data["cnt|total"];

        const key = "cnt|" + annotbar_rank_select.value + "|" + annotbar_select.value;
        console.log(key);
        if (y1_select.value=="%"){
            for (var i = 0; i < total.length; ++i) {
                pdata['bar|selected'][i] = (ddata[key][i]/total[i])*100;
                pdata['bar|unassigned'][i] = (ddata['cnt|unassigned'][i]/total[i])*100;
                pdata['bar|others'][i] = ((ddata['cnt|assigned'][i] - ddata[key][i])/total[i])*100;
            }
        }else{
            for (var i = 0; i < total.length; ++i) {
                pdata['bar|selected'][i] = ddata[key][i];
                pdata['bar|unassigned'][i] = ddata['cnt|unassigned'][i];
                pdata['bar|others'][i] = ddata['cnt|assigned'][i] - ddata[key][i];
            }
        }
        if (y1_select.value=="%")
            y_range.end = 100;
        else
            y_range.end = max_total_count;
        cds_p_samplebars.change.emit();
        ''')

    change_y_obs_label_callback = CustomJS(
        args=dict(yaxis=ele["samplebars"]["fig"].yaxis[1]),
        code='''
        console.log("change_y_obs_label_callback");
        yaxis.axis_label = this.value + " observations";
        ''')

    change_y_counts_label_callback = CustomJS(
        args=dict(yaxis=ele["samplebars"]["fig"].yaxis[0]),
        code='''
        console.log("change_y_counts_label_callback");
        yaxis.axis_label = this.value + " counts";
        ''')

    sort_groupby_callback = CustomJS(
        args=dict(sort_select=ele["samplebars"]["wid"]["sort_select"],
                  groupby1_select=ele["samplebars"]["wid"]["groupby1_select"],
                  groupby2_select=ele["samplebars"]["wid"]["groupby2_select"],
                  cds_p_samplebars=cds_p_samplebars,
                  cds_d_samples=cds_d_samples,
                  cds_d_metadata=cds_d_metadata,
                  samplebars=ele["samplebars"]["fig"]),
        code='''
        console.log("sort_groupby_callback");
        const samples = cds_d_samples.data["index"];

        // Define value from Sort by select
        var sort_col;
        if (sort_select.value=="input_order"){
            sort_col = cds_d_samples.data["aux|input_order"];
        }else if (sort_select.value=="#"){
            sort_col = cds_d_samples.data["cnt|total"];
        }else if (sort_select.value=="selected_annotation"){
            sort_col = cds_p_samplebars.data["bar|selected"];
        }else if (sort_select.value.startsWith("metadata_num|")){
            sort_col = cds_d_metadata.data[sort_select.value.replace('metadata_num|','')];
        }else if (sort_select.value.startsWith("tax|")){
            sort_col = cds_p_samplebars.data[sort_select.value];
        }

        // Factors can be: index (sort_col|index), [md1, index] or [md1, md2, index]
        var factors;
        var sorted_factors;
        // If group by is selected, use as first sort factor
        if(groupby1_select.value!="none"){
            var groupby_col1 = cds_d_metadata.data[groupby1_select.value.replace('metadata_cat|','')];

            // Zip sample index and metadata field to create nested factors
            factors = groupby_col1.map(function(m, i) {
              return [m, samples[i]];
            });

            // second grouping level
            if(groupby2_select.value!="none" && groupby2_select.value!=groupby1_select.value){

                var groupby_col2 = cds_d_metadata.data[groupby2_select.value.replace('metadata_cat|','')];

                factors = groupby_col2.map(function(m, i) {
                  return [m, groupby_col1[i], samples[i]];
                });

                sorted_factors = grimer_sort(factors, sort_col, "numeric", false, groupby_col1, groupby_col2);
            }else{
                sorted_factors = grimer_sort(factors, sort_col, "numeric", false, groupby_col1);
            }

        }else{
            // Single factors, just use the sample index
            factors = samples;
            sorted_factors = grimer_sort(factors, sort_col, "numeric", false);
        }

        // Change value of the factors on the samplebars cds
        cds_p_samplebars.data["aux|factors"] = factors;

        // Plot sorted factors
        samplebars.x_range.factors = sorted_factors;
        ''')

    plot_obs_callback = CustomJS(
        args=dict(y2_select=ele["samplebars"]["wid"]["y2_select"],
                  cds_p_samplebars=cds_p_samplebars,
                  cds_d_samples=cds_d_samples,
                  cds_p_obstable=cds_p_obstable,
                  cds_d_sampleobs=cds_d_sampleobs,
                  samplebars=ele["samplebars"]["fig"],
                  y_range=ele["samplebars"]["fig"].extra_y_ranges['obs'],
                  min_obs_perc=min_obs_perc,
                  max_total_count=max_total_count,
                  active_ranks=active_ranks),
        code='''
        console.log("plot_obs_callback");
        // get selected row from obstable [0 to get just the first]
        var row = cds_p_obstable.selected.indices[0];
        if (row!=undefined){
            // get totals
            const total = cds_d_samples.data["cnt|total"];
            // for each rank
            for(let r = 0; r < active_ranks.length; r++){
                // get taxid of the rank
                let rank_taxid = cds_p_obstable.data["tax|"+active_ranks[r]][row];
                // for each sample
                for (var i = 0; i < cds_d_sampleobs.length; i++) {
                    let val = 0;
                    // if taxid for the rank exists, [transform and] copy  values over to the cds_p_samplebars
                    if (rank_taxid){
                        val = cds_d_sampleobs.data[rank_taxid][i];
                        if(val>0){
                            if (y2_select.value=="#"){
                                val = cds_d_sampleobs.data[rank_taxid][i];
                            }else if (y2_select.value=="%"){
                                val = (cds_d_sampleobs.data[rank_taxid][i]/total[i])*100;
                            }else if (y2_select.value=="log10(%)"){
                                val = Math.log10((cds_d_sampleobs.data[rank_taxid][i]/total[i])*100);
                            }else if (y2_select.value=="log10(#)"){
                                val = Math.log10(cds_d_sampleobs.data[rank_taxid][i]);
                            }
                        }
                    }
                    if(val==0) val=NaN; //do not print if 0
                    cds_p_samplebars.data["tax|" + active_ranks[r]][i] = val;
                }
            }
            // Adjust ranges
            if (y2_select.value=="#"){
                y_range.start = 0;
                y_range.end = max_total_count;
            }else if (y2_select.value=="%"){
                y_range.start = 0;
                y_range.end = 100;
            }else if (y2_select.value=="log10(%)"){
                y_range.start = Math.log10(min_obs_perc*100);
                y_range.end = Math.log10(100);
            }else if (y2_select.value=="log10(#)"){
                y_range.start = 0;
                y_range.end = Math.log10(max_total_count);
            }
            cds_p_samplebars.change.emit();
        }
    ''')

    change_text_legend_obs_callback = CustomJS(
        args=dict(cds_p_obstable=cds_p_obstable,
                  legend_obs=ele["samplebars"]["legend_obs"],
                  active_ranks=active_ranks),
        code='''
        console.log("change_text_legend_obs_callback");
        // selected row
        var row = cb_obj.indices[0];
        for(let r = 0; r < active_ranks.length; r++){
            let taxid = cds_p_obstable.data["tax|"+active_ranks[r]][row];
            if (taxid){
                legend_obs.items[r].label = active_ranks[r] + "|" + cds_p_obstable.data['col|name'][cds_p_obstable.data['index'].indexOf(taxid)];
            }else{
                legend_obs.items[r].label = active_ranks[r];
            }
        }
        ''')

    change_text_legend_bars_callback = CustomJS(
        args=dict(annotbar_select=ele["samplebars"]["wid"]["annotbar_select"],
                  annotbar_rank_select=ele["samplebars"]["wid"]["annotbar_rank_select"],
                  legend_bars=ele["samplebars"]["legend_bars"]),
        code='''
        console.log("change_text_legend_bars_callback");
        legend_bars.items[0].label = annotbar_rank_select.value + "|" + annotbar_select.value;
        ''')

    load_infopanel = CustomJS(
        args=dict(infopanel=ele["infopanel"]["textarea"],
                  cds_p_obstable=cds_p_obstable,
                  dict_d_refs=dict_d_refs),
        code='''
        console.log("load_infopanel");

        // selected row
        var row = cb_obj.indices[0];
        const name = cds_p_obstable.data['col|name'][row];
        const rank = cds_p_obstable.data['col|rank'][row];
        const taxid = cds_p_obstable.data['index'][row];

        var text = "";
        text+="Name: " + name;
        text+="\\n";
        text+="Id: " + taxid;
        text+="\\n";
        text+="Rank: " + rank;
        text+="\\n";
        text+="\\n"

        for (var source in dict_d_refs[taxid]){
            text+="[ "+source+" ]";
            text+="\\n";
            for (var desc in dict_d_refs[taxid][source]){
                text+=desc+":";
                text+="\\n";
                for (var ref in dict_d_refs[taxid][source][desc]){
                    text+=" (" + dict_d_refs[taxid][source][desc][ref][0] + ") ";
                    text+=dict_d_refs[taxid][source][desc][ref][1].replace("{}", taxid);
                    text+="\\n";
                }
            }
            text+="\\n";
        }
        infopanel.value = text;
        ''')

    decontam_callback = CustomJS(
        args=dict(cds_d_samples=cds_d_samples,
                  cds_p_obstable=cds_p_obstable,
                  cds_d_sampleobs=cds_d_sampleobs,
                  cds_p_decontam=cds_p_decontam,
                  cds_p_decontam_models=cds_p_decontam_models,
                  cds_d_decontam=cds_d_decontam,
                  pvalue_input=ele["decontam"]["wid"]["pvalue_input"]),
        code='''
        console.log("decontam_callback");
        // selected row
        const row = cb_obj.indices[0];
        const taxid = cds_p_obstable.data["index"][row];
        const total = cds_d_samples.data["cnt|total"];
        for(let i = 0; i < cds_p_decontam.data["counts"].length; i++){
            cds_p_decontam.data["counts"][i] = cds_d_sampleobs.data[taxid][i]/total[i];
        }
        cds_p_decontam.change.emit();

        const lines = cds_d_decontam.data[taxid];
        if (lines!=undefined){
            cds_p_decontam_models.data["y_cont"] = [lines[0], lines[1]];
            cds_p_decontam_models.data["y_noncont"] = [lines[2], lines[2]];
            pvalue_input.value = lines[3].toString();
        }else{
            cds_p_decontam_models.data["y_cont"] = [];
            cds_p_decontam_models.data["y_noncont"] = [];
            pvalue_input.value = "";
        }
        cds_p_decontam_models.change.emit();
        ''')

    mgnify_callback = CustomJS(
        args=dict(mgnify_fig=ele["mgnify"]["fig"],
                  biome_spinner=ele["mgnify"]["wid"]["biome_spinner"],
                  mgnify_filter=ele["mgnify"]["filter"],
                  cds_p_obstable=cds_p_obstable,
                  cds_p_mgnify=cds_p_mgnify),
        code='''
        console.log("mgnify_callback");
        // selected row
        const row = cds_p_obstable.selected.indices[0];
        const indices = [];
        if (row!=undefined){
            const taxid = cds_p_obstable.data["index"][row];
            for(let i = 0; i < cds_p_mgnify.length; i++){
                if(cds_p_mgnify.data["taxa"][i]==taxid &&
                   cds_p_mgnify.data["level"][i]==biome_spinner.value.toString()){
                    indices.push(i);
                }
            }
        }
        mgnify_filter.indices = indices;
        cds_p_mgnify.change.emit();
        ''')

    references_callback = CustomJS(
        args=dict(references_fig=ele["references"]["fig"],
                  references_filter=ele["references"]["filter"],
                  references_select=ele["references"]["wid"]["references_select"],
                  cds_p_obstable=cds_p_obstable,
                  cds_p_references=cds_p_references,
                  active_ranks=active_ranks),
        code='''
        console.log("references_callback");
        // selected row
        const row = cds_p_obstable.selected.indices[0];
        const indices = [];
        if (row!=undefined){
            for(let i = 0; i < cds_p_references.length; i++){
                // for each rank
                for(let r = 0; r < active_ranks.length; r++){
                    // get taxid of the rank
                    let rank_obs = cds_p_obstable.data["tax|"+active_ranks[r]][row];
                    if(cds_p_references.data["obs"][i]==rank_obs &&
                       cds_p_references.data["rank"][i]==active_ranks[r] &&
                       cds_p_references.data["ref"][i]==references_select.value){
                        indices.push(i);
                    }
                }
            }
        }
        references_filter.indices = indices;
        cds_p_references.change.emit();
        ''')

    obstable_callbacks = [plot_obs_callback, change_text_legend_obs_callback, sort_groupby_callback, load_infopanel, references_callback]
    if cds_p_decontam:
        obstable_callbacks.append(decontam_callback)
    if cds_p_mgnify:
        obstable_callbacks.append(mgnify_callback)
    cds_p_obstable.selected.js_on_change('indices', *obstable_callbacks)

    ele["samplebars"]["wid"]["sort_select"].js_on_change('value', sort_groupby_callback)
    ele["samplebars"]["wid"]["groupby1_select"].js_on_change('value', sort_groupby_callback)
    ele["samplebars"]["wid"]["groupby2_select"].js_on_change('value', sort_groupby_callback)
    ele["samplebars"]["wid"]["annotbar_select"].js_on_change('value', bar_select_callback, change_text_legend_bars_callback, sort_groupby_callback)
    ele["samplebars"]["wid"]["annotbar_rank_select"].js_on_change('value', bar_select_callback, change_text_legend_bars_callback, sort_groupby_callback)
    ele["samplebars"]["wid"]["y1_select"].js_on_change('value', bar_select_callback, change_y_counts_label_callback, sort_groupby_callback)
    ele["samplebars"]["wid"]["y2_select"].js_on_change('value', plot_obs_callback, change_y_obs_label_callback, sort_groupby_callback)
    ele["mgnify"]["wid"]["biome_spinner"].js_on_change('value', mgnify_callback)
    ele["references"]["wid"]["references_select"].js_on_change('value', references_callback)

def link_heatmap_widgets(ele,
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
                         cds_p_heatmap):

    x_dendro_callback = CustomJS(
        args=dict(rank_select=ele["heatmap"]["wid"]["rank_select"],
                  x_sort_select=ele["heatmap"]["wid"]["x_sort_select"],
                  x_method_select=ele["heatmap"]["wid"]["x_method_select"],
                  cds_p_dendro_x=cds_p_dendro_x,
                  dict_d_dedro_x=dict_d_dedro_x),
        code='''
        console.log("x_dendro_callback");
        if (x_sort_select.value.startsWith("metric|")){
            const key = rank_select.value+"|"+x_method_select.value+"|"+x_sort_select.value.replace("metric|","");
            cds_p_dendro_x.data = {"x": dict_d_dedro_x[key+"|x"],
                                 "y": dict_d_dedro_x[key+"|y"],
                                 "c": dict_d_dedro_x[key+"|c"]};
            // Enable method select
            x_method_select.disabled=false;
        }else{
            cds_p_dendro_x.data = {"x": [], "y": [], "c": []};
            // Disable method select
            x_method_select.disabled=true;
        }
        ''')

    x_select_callback = CustomJS(
        args=dict(heatmap=ele["heatmap"]["fig"],
                  rank_select=ele["heatmap"]["wid"]["rank_select"],
                  x_method_select=ele["heatmap"]["wid"]["x_method_select"],
                  x_sort_select=ele["heatmap"]["wid"]["x_sort_select"],
                  dict_d_hcluster_x=dict_d_hcluster_x,
                  cds_p_annotations=cds_p_annotations,
                  cds_p_obstable=cds_p_obstable),
        code='''
        console.log("x_select_callback");
        const rank = rank_select.value;
        var sorted_factors = [];
        if (x_sort_select.value=="none"){
            // None
            sorted_factors = dict_d_hcluster_x["default|" + rank];
        }else if (x_sort_select.value.startsWith("metric|")){
            // Clustering
            // Get sorted elements based on rank|method|metric
            const key = rank+"|"+x_method_select.value+"|"+x_sort_select.value.replace("metric|","");
            sorted_factors = dict_d_hcluster_x[key];
        }else{
            // Sorting
            var factors =  [];
            var sort_col = [];
            if (x_sort_select.value=="counts"){
                for (let i = 0; i < cds_p_obstable.data["index"].length; i++) {
                    if(cds_p_obstable.data["col|rank"][i]==rank){
                        factors.push(cds_p_obstable.data["index"][i])
                        sort_col.push(cds_p_obstable.data["col|total_counts"][i]);
                    }
                }
                sorted_factors = grimer_sort(factors, sort_col, "numeric", false);
            }else if (x_sort_select.value=="observations"){
                for (let i = 0; i < cds_p_obstable.data["index"].length; i++) {
                    if(cds_p_obstable.data["col|rank"][i]==rank){
                        factors.push(cds_p_obstable.data["index"][i])
                        sort_col.push(cds_p_obstable.data["col|name"][i]);
                    }
                }
                sorted_factors = grimer_sort(factors, sort_col, "string", false);
            }else{
                // copy array of factors
                factors =  [...dict_d_hcluster_x["default|" + rank]];
                const annot = x_sort_select.value.replace("annot|","");
                for (let i = 0; i < cds_p_annotations.data["index"].length; i++) {
                    if (cds_p_annotations.data["rank"][i]==rank && cds_p_annotations.data["annot"][i]==annot) {
                        // if annotation is present in selected rank + annot
                        const index = factors.indexOf(cds_p_annotations.data["index"][i]);
                        if (index > -1) {
                          sorted_factors.push(factors[index]); // add to the sorted_factors
                          factors.splice(index, 1); //remove from factors
                        }
                    }
                }
                // join annotated and left-overs
                sorted_factors = sorted_factors.concat(factors);
            }
        }
        heatmap.x_range.factors = sorted_factors;
        ''')

    y_dendro_callback = CustomJS(
        args=dict(rank_select=ele["heatmap"]["wid"]["rank_select"],
                  y_sort_select=ele["heatmap"]["wid"]["y_sort_select"],
                  y_method_select=ele["heatmap"]["wid"]["y_method_select"],
                  cds_p_dendro_y=cds_p_dendro_y,
                  dict_d_dedro_y=dict_d_dedro_y),
        code='''
        console.log("y_dendro_callback");
        if (y_sort_select.value.startsWith("metric|")){
            const key = rank_select.value+"|"+y_method_select.value+"|"+y_sort_select.value.replace("metric|","");
            cds_p_dendro_y.data = {"x": dict_d_dedro_y[key+"|x"],
                                 "y": dict_d_dedro_y[key+"|y"],
                                 "c": dict_d_dedro_y[key+"|c"]};
            // Enable method select
            y_method_select.disabled=false;
        }else{
            cds_p_dendro_y.data = {"x": [], "y": [], "c": []};
            // Disable method select
            y_method_select.disabled=true;
        }
        ''')

    y_select_callback = CustomJS(
        args=dict(heatmap=ele["heatmap"]["fig"],
                  cds_d_samples=cds_d_samples,
                  cds_d_metadata=cds_d_metadata,
                  rank_select=ele["heatmap"]["wid"]["rank_select"],
                  y_method_select=ele["heatmap"]["wid"]["y_method_select"],
                  y_sort_select=ele["heatmap"]["wid"]["y_sort_select"],
                  dict_d_hcluster_y=dict_d_hcluster_y),
        code='''
        console.log("y_select_callback");
        var sorted_factors = [];
        if (y_sort_select.value=="none"){
            // None
            sorted_factors = dict_d_hcluster_y["default"];
        }else if (y_sort_select.value.startsWith("metric|")){
            // Clustering
            // Get sorted elements based on rank|method|metric
            const key = rank_select.value+"|"+y_method_select.value+"|"+y_sort_select.value.replace("metric|","");
            sorted_factors = dict_d_hcluster_y[key];
        }else{
            // Sorting
            if (y_sort_select.value=="counts"){
                sorted_factors = grimer_sort(cds_d_samples.data["index"], cds_d_samples.data["cnt|total"], "numeric", false);
            }else if (y_sort_select.value=="samples"){
                sorted_factors = grimer_sort(cds_d_samples.data["index"], cds_d_samples.data["index"], "string", false);
            }else if (y_sort_select.value.startsWith("metadata_cat|")){
                sorted_factors = grimer_sort(cds_d_samples.data["index"], cds_d_metadata.data[y_sort_select.value.replace("metadata_cat|","")], "string", false);
            }else if (y_sort_select.value.startsWith("metadata_num|")){
                sorted_factors = grimer_sort(cds_d_samples.data["index"], cds_d_metadata.data[y_sort_select.value.replace("metadata_num|","")], "numeric", false);
            }
        }
        heatmap.y_range.factors = sorted_factors;
        ''')

    toggle_label_callback = CustomJS(
        args=dict(cds_p_heatmap=cds_p_heatmap,
                  xaxis=ele["heatmap"]["fig"].xaxis[0],
                  yaxis=ele["heatmap"]["fig"].yaxis[0]),
        code='''
        if(this.active.includes(0)){
            xaxis.major_label_text_font_size = "12px";
            xaxis.major_tick_line_color="black";
        }else{
            xaxis.major_label_text_font_size = "0px";
            xaxis.major_tick_line_color=null;
        }
        if(this.active.includes(1)){
            yaxis.major_label_text_font_size = "12px";
            yaxis.major_tick_line_color="black";
        }else{
            yaxis.major_label_text_font_size = "0px";
            yaxis.major_tick_line_color=null;
        }
        ''')

    ele["heatmap"]["wid"]["toggle_label_heatmap"].js_on_click(toggle_label_callback)
    ele["heatmap"]["wid"]["rank_select"].js_on_change('value', x_select_callback, x_dendro_callback, y_select_callback, y_dendro_callback)
    ele["heatmap"]["wid"]["x_method_select"].js_on_change('value', x_select_callback, x_dendro_callback)
    ele["heatmap"]["wid"]["x_sort_select"].js_on_change('value', x_select_callback, x_dendro_callback)
    ele["heatmap"]["wid"]["y_method_select"].js_on_change('value', y_select_callback, y_dendro_callback)
    ele["heatmap"]["wid"]["y_sort_select"].js_on_change('value', y_select_callback, y_dendro_callback)


def link_metadata_widgets(ele, cds_p_metadata, cds_d_metadata, max_metadata_cols):
    metadata_multiselect_callback = CustomJS(
        args=dict(metadata_heatmap=ele["metadata"]["fig"],
                  max_metadata_cols=max_metadata_cols,
                  cds_p_metadata=cds_p_metadata,
                  cds_d_metadata=cds_d_metadata),
        code='''
        const index_len = cds_d_metadata.data["index"].length;
        console.log(cds_d_metadata.data)
        var x_factors = [];
        var empty_y_values = new Array(index_len);
        for (var i = 0; i < index_len; ++i) empty_y_values[i]=["", ""];
        for(var s=0; s < max_metadata_cols; ++s){
            if (s<this.value.length){
                var selected = this.value[s];
                var y_values = new Array(index_len);
                for (var i = 0; i < index_len; ++i){
                    y_values[i]=[selected, cds_d_metadata.data[selected][i].toString()];
                }
                console.log(y_values);
                cds_p_metadata.data["md" + s.toString()] = y_values;
                x_factors.push("md" + s.toString());
            }else{
                cds_p_metadata.data["md" + s.toString()] = empty_y_values;
            }
        }
        metadata_heatmap.x_range.factors = x_factors;
        cds_p_metadata.change.emit();
        ''')

    if cds_d_metadata:
        ele["metadata"]["wid"]["metadata_multiselect"].js_on_change('value', metadata_multiselect_callback)


def link_obstable_filter(ele, cds_p_obstable, active_ranks):
    filter_callback = CustomJS(
        args=dict(cds_p_obstable=cds_p_obstable,
                  active_ranks=active_ranks,
                  widgets_filter=ele["obstable"]["widgets_filter"],
                  frequency_spinner=ele["obstable"]["wid"]["frequency_spinner"],
                  counts_perc_avg_spinner=ele["obstable"]["wid"]["counts_perc_avg_spinner"],
                  total_counts_spinner=ele["obstable"]["wid"]["total_counts_spinner"],
                  name_multichoice=ele["obstable"]["wid"]["name_multichoice"],
                  ),
        code='''
        const indices = [];
        for (var i = 0; i < cds_p_obstable.length; i++) {
            if (cds_p_obstable.data['col|frequency_perc'][i] < (frequency_spinner.value/100)){
                continue;
            }
            if (cds_p_obstable.data['col|counts_perc_avg'][i] < (counts_perc_avg_spinner.value/100)){
                continue;
            }
            if (cds_p_obstable.data['col|total_counts'][i] < (total_counts_spinner.value)){
                continue;
            }
            if (name_multichoice.value.length > 0 ){
                var found = false;
                for(let r = 0; r < active_ranks.length; r++){
                    // Compare all names on multichoice (array) against cell
                    if (name_multichoice.value.indexOf(cds_p_obstable.data["tax|"+active_ranks[r]][i]) >= 0){
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    continue;
                }
            }
            indices.push(i);
        }
        console.log(cds_p_obstable);
        console.log(widgets_filter);
        widgets_filter.indices = indices;
        cds_p_obstable.change.emit();
        console.log(cds_p_obstable);
        console.log(widgets_filter);
    ''')
    ele["obstable"]["wid"]["frequency_spinner"].js_on_change('value', filter_callback)
    ele["obstable"]["wid"]["counts_perc_avg_spinner"].js_on_change('value', filter_callback)
    ele["obstable"]["wid"]["total_counts_spinner"].js_on_change('value', filter_callback)
    ele["obstable"]["wid"]["name_multichoice"].js_on_change('value', filter_callback)


def link_correlation_widgets(ele, cds_p_correlation):
    rank_select_callback = CustomJS(
        args=dict(correlation=ele["correlation"]["fig"],
                  cds_p_correlation=cds_p_correlation),
        code='''
        console.log("rank_select_callback");
        const factors = new Set();
        for(let i = 0; i < cds_p_correlation.data["index"].length; i++){
            if(cds_p_correlation.data["rank"][i]==this.value){
                factors.add(cds_p_correlation.data["index"][i]);
            }
        }
        correlation.x_range.factors = [...factors];
        correlation.y_range.factors = [...factors].reverse();
        ''')

    filter_callback = CustomJS(
        args=dict(rho_filter=ele["correlation"]["rho_filter"],
                  neg_slider=ele["correlation"]["wid"]["neg_slider"],
                  pos_slider=ele["correlation"]["wid"]["pos_slider"],
                  cds_p_correlation=cds_p_correlation),
        code='''
        console.log("filter_callback");
        const indices = [];
        for (var i = 0; i < cds_p_correlation.data["index"].length; i++) {
            const rho = cds_p_correlation.data["rho"][i];
            if ((rho >= neg_slider.value[0] && rho <= neg_slider.value[1]) ||
                (rho >= pos_slider.value[0] && rho <= pos_slider.value[1]))
            {
                indices.push(i)
            }
        }
        rho_filter.indices = indices;
        cds_p_correlation.change.emit();
        ''')

    ele["correlation"]["wid"]["pos_slider"].js_on_change('value', filter_callback)
    ele["correlation"]["wid"]["neg_slider"].js_on_change('value', filter_callback)
    ele["correlation"]["wid"]["rank_select"].js_on_change('value', rank_select_callback)


def link_obsbars_widgets(ele, cds_p_obsbars, dict_d_topobs, cds_d_sampleobs, cds_d_samples, top_obs_bars, dict_d_taxname, cds_d_metadata):
    rank_select_callback = CustomJS(
        args=dict(sort_select=ele["obsbars"]["wid"]["sort_select"],
                  legend=ele["obsbars"]["legend"],
                  cds_p_obsbars=cds_p_obsbars,
                  cds_d_sampleobs=cds_d_sampleobs,
                  cds_d_samples=cds_d_samples,
                  dict_d_topobs=dict_d_topobs,
                  dict_d_taxname=dict_d_taxname,
                  top_obs_bars=top_obs_bars),
        code='''
        console.log("rank_select_callback");
        const rank = this.value;
        const n_sample = cds_p_obsbars.data["index"].length;
        const total = cds_d_samples.data["cnt|total"];
        const unassigned = cds_d_samples.data["cnt|unassigned"];

        var empty_values = new Array(n_sample);
        for (var i = 0; i < n_sample; ++i) empty_values[i]=0;

        var sum_assigned = new Array(n_sample);
        for (var i = 0; i < n_sample; ++i) sum_assigned[i]=0;

        sort_select.options["Observations"] = [];

        for(let i = 0; i < top_obs_bars; i++){
            var taxid = dict_d_topobs[rank][i];
            var sum_bars = 0;
            if (taxid!=undefined){
                for(let s = 0; s<n_sample; s++){
                    cds_p_obsbars.data[i.toString()][s] = (cds_d_sampleobs.data[taxid][s]/total[s])*100;
                    // sum counts for sample
                    sum_assigned[s]+=cds_d_sampleobs.data[taxid][s];
                    // update legend label
                    legend.items[i].label = i.toString() + ") " + dict_d_taxname[taxid];
                }
                // not sync with gui
                // https://github.com/bokeh/bokeh/issues/10211
                // sort_select.options["Observations"].push(["col|" + i.toString(), dict_d_taxname[taxid]])
            }else{
                cds_p_obsbars.data[i.toString()] = [...empty_values];
                legend.items[i].label = null;
            }
        }
        for(let s = 0; s<n_sample; s++){
            cds_p_obsbars.data["others"][s] = ((total[s] - unassigned[s] - sum_assigned[s])/total[s])*100;
        }
        cds_p_obsbars.change.emit();
        ''')

    sort_groupby_callback = CustomJS(
        args=dict(obsbars=ele["obsbars"]["fig"],
                  sort_select=ele["obsbars"]["wid"]["sort_select"],
                  groupby1_select=ele["obsbars"]["wid"]["groupby1_select"],
                  groupby2_select=ele["obsbars"]["wid"]["groupby2_select"],
                  cds_p_obsbars=cds_p_obsbars,
                  cds_d_samples=cds_d_samples,
                  cds_d_metadata=cds_d_metadata),
        code='''
        console.log("sort_groupby_callback");

        // Define value from Sort by select
        var sort_col;
        var annot_samples = cds_d_samples.data["index"];
        if (sort_select.value=="input_order"){
            sort_col = cds_d_samples.data["aux|input_order"];
        }else if (sort_select.value.startsWith("metadata_num|")){
            sort_col = cds_d_metadata.data[sort_select.value.replace('metadata_num|','')];

            // Annotate label with value
            var annot_samples = annot_samples.map(function(s, i) {
              return sort_col[i] + " | " + s;
            });

        }else if (sort_select.value.startsWith("col|")){
            sort_col = cds_p_obsbars.data[sort_select.value.replace('col|','')];
        }

        // Factors can be: index (sort_col|index), [md1, index] or [md1, md2, index]
        var factors;
        var sorted_factors;
        // If group by is selected, use as first sort factor
        if(groupby1_select.value!="none"){
            var groupby_col1 = cds_d_metadata.data[groupby1_select.value.replace('metadata_cat|','')];

            // Zip sample index and metadata field to create nested factors
            factors = groupby_col1.map(function(m, i) {
              return [m, annot_samples[i]];
            });

            // second grouping level
            if(groupby2_select.value!="none" && groupby2_select.value!=groupby1_select.value){

                var groupby_col2 = cds_d_metadata.data[groupby2_select.value.replace('metadata_cat|','')];

                factors = groupby_col2.map(function(m, i) {
                  return [m, groupby_col1[i], annot_samples[i]];
                });

                sorted_factors = grimer_sort(factors, sort_col, "numeric", false, groupby_col1, groupby_col2);
            }else{
                sorted_factors = grimer_sort(factors, sort_col, "numeric", false, groupby_col1);
            }

        }else{
            // Single factors, just use the sample index
            factors = annot_samples;
            sorted_factors = grimer_sort(factors, sort_col, "numeric", false);
        }

        // Change value of the factors on the obsbars cds
        cds_p_obsbars.data["factors"] = factors;

        // Plot sorted factors
        obsbars.x_range.factors = sorted_factors;
        ''')

    toggle_label_callback = CustomJS(
        args=dict(xaxis=ele["obsbars"]["fig"].xaxis[0]),
        code='''
        if(this.active.includes(0)){
            xaxis.major_label_text_font_size = "0px";
            xaxis.major_tick_line_color=null;
        }else{
            xaxis.major_label_text_font_size = "12px";
            xaxis.major_tick_line_color="black";
        }
        ''')

    ele["obsbars"]["wid"]["toggle_label"].js_on_click(toggle_label_callback)
    ele["obsbars"]["wid"]["groupby1_select"].js_on_change('value', sort_groupby_callback)
    ele["obsbars"]["wid"]["groupby2_select"].js_on_change('value', sort_groupby_callback)
    ele["obsbars"]["wid"]["sort_select"].js_on_change('value', sort_groupby_callback)
    ele["obsbars"]["wid"]["rank_select"].js_on_change('value', rank_select_callback, sort_groupby_callback)
