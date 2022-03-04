from bokeh.models import CustomJS


def link_obstable_samplebars(ele,
                             cds_m_obstable,
                             cds_p_samplebars,
                             cds_d_samples,
                             dict_d_sampleobs,
                             cds_d_metadata,
                             cds_p_decontam,
                             cds_p_decontam_models,
                             cds_d_decontam,
                             cds_p_references,
                             active_ranks,
                             min_obs_perc,
                             max_total_count,
                             cds_p_mgnify,
                             dict_d_refs,
                             dict_d_taxname):

    bar_select_callback = CustomJS(
        args=dict(y1_select=ele["samplebars"]["wid"]["y1_select"],
                  annotbar_select=ele["samplebars"]["wid"]["annotbar_select"],
                  annotbar_rank_select=ele["samplebars"]["wid"]["annotbar_rank_select"],
                  cds_p_samplebars=cds_p_samplebars,
                  cds_d_samples=cds_d_samples,
                  y_range=ele["samplebars"]["fig"].y_range,
                  max_total_count=max_total_count),
        code='''

        const pdata = cds_p_samplebars.data;
        const ddata = cds_d_samples.data;
        const total = cds_d_samples.data["cnt|total"];

        const key = "cnt|" + annotbar_rank_select.value + "|" + annotbar_select.value;

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
        yaxis.axis_label = this.value + " observations";
        ''')

    change_y_counts_label_callback = CustomJS(
        args=dict(yaxis=ele["samplebars"]["fig"].yaxis[0]),
        code='''
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
        // Define value from Sort by select
        var sort_col;
        var annot_samples = cds_d_samples.data["index"];
        if (sort_select.value=="input_order"){
            sort_col = cds_d_samples.data["aux|input_order"];
        }else if (sort_select.value=="counts"){
            sort_col = cds_d_samples.data["cnt|total"];
        }else if (sort_select.value=="selected_annotation"){
            sort_col = cds_p_samplebars.data["bar|selected"];
        }else if (sort_select.value.startsWith("metadata_num|")){
            sort_col = cds_d_metadata.data[sort_select.value.replace('metadata_num|','')];

            // Annotate label with value
            var annot_samples = annot_samples.map(function(s, i) {
              return sort_col[i] + " | " + s;
            });

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
              return [m, annot_samples[i]];
            });

            // second grouping level
            if(groupby2_select.value!="none" && groupby2_select.value!=groupby1_select.value){

                var groupby_col2 = cds_d_metadata.data[groupby2_select.value.replace('metadata_cat|','')];

                factors = groupby_col2.map(function(m, i) {
                  return [groupby_col1[i], m, annot_samples[i]];
                });

                sorted_factors = grimer_sort(factors, sort_col, "numeric", false, groupby_col2, groupby_col1);
            }else{
                sorted_factors = grimer_sort(factors, sort_col, "numeric", false, groupby_col1);
            }

        }else{
            // Single factors, just use the sample index
            factors = annot_samples;
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
                  cds_m_obstable=cds_m_obstable,
                  dict_d_sampleobs=dict_d_sampleobs,
                  y_range=ele["samplebars"]["fig"].extra_y_ranges['obs'],
                  min_obs_perc=min_obs_perc,
                  max_total_count=max_total_count,
                  active_ranks=active_ranks),
        code='''
        // get selected row from obstable [0 to get just the first]
        var row = cds_m_obstable.selected.indices[0];
        if (row!=undefined){
            // get totals
            const total = cds_d_samples.data["cnt|total"];
            // for each rank
            for(let r = 0; r < active_ranks.length; r++){
                // get rank
                let rank = active_ranks[r];
                // get taxid of the rank
                let taxid = cds_m_obstable.data["tax|"+rank][row];
                // for each sample
                for (var i = 0; i < cds_d_samples.length; i++) {
                    let sample = cds_d_samples.data["index"][i];
                    let val = 0;
                    // if taxid exists in the lineage, [transform and] copy values over to the cds_p_samplebars
                    if (taxid){
                        val = dict_d_sampleobs[rank][taxid][sample];
                        if(val>0){
                            if (y2_select.value=="%"){
                                val = (val/total[i])*100;
                            }else if (y2_select.value=="log10(%)"){
                                val = Math.log10((val/total[i])*100);
                            }else if (y2_select.value=="log10(#)"){
                                val = Math.log10(val);
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
        args=dict(cds_m_obstable=cds_m_obstable,
                  legend_obs=ele["samplebars"]["legend_obs"],
                  samplebars=ele["samplebars"]["fig"],
                  active_ranks=active_ranks),
        code='''
        // selected row
        const row = cb_obj.indices[0];
        const selected_rank = cds_m_obstable.data['col|rank'][row];

        for(let r = 0; r < active_ranks.length; r++){
            let taxid = cds_m_obstable.data["tax|"+active_ranks[r]][row];
            if (taxid){
                legend_obs.items[r].label = active_ranks[r] + "|" + cds_m_obstable.data['col|name'][cds_m_obstable.data['index'].indexOf(taxid)];
            }else{
                legend_obs.items[r].label = active_ranks[r];
            }
            // activate only selected rank
            if(active_ranks[r]==selected_rank){
                samplebars.renderers[r+3].visible=true;
            //}else{
            //    samplebars.renderers[r+3].visible=false;
            }
        }
        ''')

    change_text_legend_bars_callback = CustomJS(
        args=dict(annotbar_select=ele["samplebars"]["wid"]["annotbar_select"],
                  annotbar_rank_select=ele["samplebars"]["wid"]["annotbar_rank_select"],
                  legend_bars=ele["samplebars"]["legend_bars"]),
        code='''
        legend_bars.items[0].label = annotbar_rank_select.value + "|" + annotbar_select.value;
        ''')

    load_infopanel = CustomJS(
        args=dict(infopanel=ele["infopanel"]["textarea"],
                  cds_m_obstable=cds_m_obstable,
                  dict_d_refs=dict_d_refs,
                  dict_d_taxname=dict_d_taxname,
                  active_ranks=active_ranks),
        code='''
        // selected row
        var row = cb_obj.indices[0];

        const name = cds_m_obstable.data['col|name'][row];
        const rank = cds_m_obstable.data['col|rank'][row];
        const taxid = cds_m_obstable.data['index'][row];

        var text = "";
        text+="[ Obs ]";
        text+="\\n";
        text+=name;
        if(taxid!=name){
            text+="\\n";
            text+="taxid: " + taxid;
        }
        text+="\\n";
        text+="[ Rank ]";
        text+="\\n";
        text+=rank;
        text+="\\n";

        var lineage = "";
        for(let r = 0; r < active_ranks.length; r++){
            var obs_lin = cds_m_obstable.data["tax|" + active_ranks[r]][row];
            if(taxid!=name){
                if(dict_d_taxname[obs_lin])
                    lineage+=dict_d_taxname[obs_lin]+" | ";
                else
                    lineage+=" | ";
            }else{
                lineage+=obs_lin+" | ";
            }
            if(active_ranks[r]==rank)
                break;
        }
        text+="[ Lineage ]";
        text+="\\n";
        text+=lineage.slice(0, -3);
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
                  cds_m_obstable=cds_m_obstable,
                  dict_d_sampleobs=dict_d_sampleobs,
                  cds_p_decontam=cds_p_decontam,
                  cds_p_decontam_models=cds_p_decontam_models,
                  cds_d_decontam=cds_d_decontam,
                  pscore_input=ele["decontam"]["wid"]["pscore_input"]),
        code='''
        // selected row
        const row = cb_obj.indices[0];
        const taxid = cds_m_obstable.data["index"][row];
        const rank = cds_m_obstable.data["col|rank"][row];
        const total = cds_d_samples.data["cnt|total"];
        for(let i = 0; i < cds_p_decontam.length; i++){
            let sample = cds_p_decontam.data["index"][i];
            if (dict_d_sampleobs[rank][taxid][sample]!=undefined){
                cds_p_decontam.data["counts"][i] = dict_d_sampleobs[rank][taxid][sample]/total[i];
            }else{
                cds_p_decontam.data["counts"][i] = 0;
            }
        }
        cds_p_decontam.change.emit();

        const lines = cds_d_decontam.data[taxid];
        if (lines!=undefined){
            cds_p_decontam_models.data["y_cont"] = [lines[0], lines[1]];
            cds_p_decontam_models.data["y_noncont"] = [lines[2], lines[2]];
            pscore_input.value = lines[3].toString();
        }else{
            cds_p_decontam_models.data["y_cont"] = [];
            cds_p_decontam_models.data["y_noncont"] = [];
            pscore_input.value = "";
        }
        cds_p_decontam_models.change.emit();
        ''')

    mgnify_callback = CustomJS(
        args=dict(mgnify_fig=ele["mgnify"]["fig"],
                  biome_spinner=ele["mgnify"]["wid"]["biome_spinner"],
                  mgnify_filter=ele["mgnify"]["filter"],
                  cds_m_obstable=cds_m_obstable,
                  cds_p_mgnify=cds_p_mgnify),
        code='''
        // selected row
        const row = cds_m_obstable.selected.indices[0];
        const indices = [];
        if (row!=undefined){
            const taxid = cds_m_obstable.data["index"][row];
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
                  cds_m_obstable=cds_m_obstable,
                  cds_p_references=cds_p_references,
                  active_ranks=active_ranks),
        code='''
        // selected row
        const row = cds_m_obstable.selected.indices[0];
        const indices = [];
        if (row!=undefined){
            for(let i = 0; i < cds_p_references.length; i++){
                // for each rank
                for(let r = 0; r < active_ranks.length; r++){
                    // get taxid of the rank
                    let rank_obs = cds_m_obstable.data["tax|"+active_ranks[r]][row];
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

    toggle_label_callback = CustomJS(
        args=dict(xaxis=ele["samplebars"]["fig"].xaxis[0]),
        code='''
        if(this.active.includes(0)){
            xaxis.major_label_text_font_size = "10px";
            xaxis.major_tick_line_color="black";
        }else{
            xaxis.major_label_text_font_size = "0px";
            xaxis.major_tick_line_color=null;
        }
        ''')

    obstable_callbacks = [plot_obs_callback, change_text_legend_obs_callback, sort_groupby_callback, load_infopanel]
    if cds_p_decontam:
        obstable_callbacks.append(decontam_callback)
    if cds_p_mgnify:
        obstable_callbacks.append(mgnify_callback)
    if ele["references"]["filter"]:
        obstable_callbacks.append(references_callback)

    cds_m_obstable.selected.js_on_change('indices', *obstable_callbacks)

    ele["samplebars"]["wid"]["sort_select"].js_on_change('value', sort_groupby_callback)
    ele["samplebars"]["wid"]["groupby1_select"].js_on_change('value', sort_groupby_callback)
    ele["samplebars"]["wid"]["groupby2_select"].js_on_change('value', sort_groupby_callback)
    ele["samplebars"]["wid"]["annotbar_select"].js_on_change('value', bar_select_callback, change_text_legend_bars_callback, sort_groupby_callback)
    ele["samplebars"]["wid"]["annotbar_rank_select"].js_on_change('value', bar_select_callback, change_text_legend_bars_callback, sort_groupby_callback)
    ele["samplebars"]["wid"]["y1_select"].js_on_change('value', bar_select_callback, change_y_counts_label_callback, sort_groupby_callback)
    ele["samplebars"]["wid"]["y2_select"].js_on_change('value', plot_obs_callback, change_y_obs_label_callback, sort_groupby_callback)
    ele["samplebars"]["wid"]["toggle_label"].js_on_click(toggle_label_callback)
    ele["mgnify"]["wid"]["biome_spinner"].js_on_change('value', mgnify_callback)
    ele["references"]["wid"]["references_select"].js_on_change('value', references_callback)


def link_heatmap_widgets(ele,
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
                         active_ranks,
                         dict_d_taxname):

    x_dendro_callback = CustomJS(
        args=dict(rank_select=ele["heatmap"]["wid"]["rank_select"],
                  x_groupby_select=ele["heatmap"]["wid"]["x_groupby_select"],
                  x_sort_select=ele["heatmap"]["wid"]["x_sort_select"],
                  dendrox=ele["dendrox"]["fig"],
                  cds_p_dendro_x=cds_p_dendro_x,
                  dict_d_dedro_x=dict_d_dedro_x),
        code='''
        if (x_groupby_select.value.startsWith("cluster|")){
            const key = rank_select.value+"|"+x_groupby_select.value.replace("cluster|","");
            cds_p_dendro_x.data = {"x": dict_d_dedro_x[key+"|x"],
                                 "y": dict_d_dedro_x[key+"|y"],
                                 "c": dict_d_dedro_x[key+"|c"]};
            x_sort_select.value="none";
            x_sort_select.disabled=true;
            dendrox.visible=true;
        }else{
            cds_p_dendro_x.data = {"x": [], "y": [], "c": []};
            x_sort_select.disabled=false;
            dendrox.visible=false;
        }
        ''')

    x_select_callback = CustomJS(
        args=dict(heatmap=ele["heatmap"]["fig"],
                  active_ranks=active_ranks,
                  rank_select=ele["heatmap"]["wid"]["rank_select"],
                  x_sort_select=ele["heatmap"]["wid"]["x_sort_select"],
                  x_groupby_select=ele["heatmap"]["wid"]["x_groupby_select"],
                  dict_d_hcluster_x=dict_d_hcluster_x,
                  cds_p_annotations=cds_p_annotations,
                  cds_m_obstable=cds_m_obstable,
                  cds_p_heatmap=cds_p_heatmap,
                  dict_d_taxname=dict_d_taxname),
        code='''
        // selected rank
        const rank = rank_select.value;

        // get index to access data from observations from cds_m_obstable
        var obs_index = [];
        for (let i = 0; i < cds_m_obstable.data["index"].length; i++) {
            if(cds_m_obstable.data["col|rank"][i]==rank){
                obs_index.push(i);
            }
        }

        var annot_obs = obs_index.map( s => cds_m_obstable.data["index"][s] );

        var sorted_factors = [];
        var dict_factors = {};
        if (x_groupby_select.value.startsWith("cluster|")){
            // Default dict_factors
            for(let i = 0; i < annot_obs.length; i++){
                dict_factors[annot_obs[i]] = annot_obs[i];
            }
            // Clustering - Get sorted elements based on rank|method|metric
            sorted_factors = dict_d_hcluster_x[rank+"|"+x_groupby_select.value.replace("cluster|","")];
        }else{
            // Define value from Sort by select
            var sort_col = [];
            var sort_col_type = "numeric";
            if (x_sort_select.value=="none"){
                sort_col = obs_index;
            }else if (x_sort_select.value=="observations"){
                sort_col = obs_index.map( s => cds_m_obstable.data["col|name"][s] );
                sort_col_type = "string";
            }else if (x_sort_select.value=="counts"){
                sort_col = obs_index.map( s => cds_m_obstable.data["col|total_counts"][s] );
            }else if (x_sort_select.value.startsWith("annot|")){
                const annot = x_sort_select.value.replace("annot|","");
                // create array with zeros, mark with one if annotation is present
                sort_col = new Array(annot_obs.length); for (let i=0; i<annot_obs.length; ++i) sort_col[i] = NaN;
                for (let i = 0; i < cds_p_annotations.data["index"].length; i++) {
                    if (cds_p_annotations.data["rank"][i]==rank && cds_p_annotations.data["annot"][i]==annot) {
                        sort_col[annot_obs.indexOf(cds_p_annotations.data["index"][i])] = cds_p_annotations.data["tv"][i];
                    }
                }
            }

            if(x_groupby_select.value=="none"){
                // Default dict_factors
                for(let i = 0; i < annot_obs.length; i++){
                    dict_factors[annot_obs[i]] = annot_obs[i];
                }
                sorted_factors = grimer_sort(annot_obs, sort_col, sort_col_type, false);
            }else if (x_groupby_select.value.startsWith("tax|")){
                const group_rank = x_groupby_select.value.replace("tax|","");
                // if grouping with a higher rank
                if(active_ranks.indexOf(rank) > active_ranks.indexOf(group_rank)){
                    // group entries without selected rank with space " "
                    var groupby_col = obs_index.map(function(s) { return cds_m_obstable.data["tax|" + group_rank][s] == "" ? " " : dict_d_taxname[cds_m_obstable.data["tax|" + group_rank][s]]; });
                    var factors = [];
                    for(let i = 0; i < annot_obs.length; i++){
                        dict_factors[annot_obs[i]] = [groupby_col[i], annot_obs[i]];
                        factors.push([groupby_col[i], annot_obs[i]]);
                    }
                    sorted_factors = grimer_sort(factors, sort_col, sort_col_type, false, groupby_col);
                }else{
                    // Default dict_factors
                    for(let i = 0; i < annot_obs.length; i++){
                        dict_factors[annot_obs[i]] = annot_obs[i];
                    }
                    // normal sort
                    sorted_factors = grimer_sort(annot_obs, sort_col, sort_col_type, false);
                }
            }
        }

        // update factors on heatmap col otherwise remove
        for (let i = 0; i < cds_p_heatmap.data["index"].length; i++) {
            if(cds_p_heatmap.data["rank"][i]==rank){
                cds_p_heatmap.data["factors_obs"][i] = dict_factors[cds_p_heatmap.data["obs"][i]];
            }else{
                cds_p_heatmap.data["factors_obs"][i] = "";
            }
        }

        for (let i = 0; i < cds_p_annotations.data["index"].length; i++) {
            if(cds_p_annotations.data["rank"][i]==rank){
                cds_p_annotations.data["factors"][i] = dict_factors[cds_p_annotations.data["index"][i]];
            }else{
                cds_p_annotations.data["factors"][i] = "";
            }
        }

        heatmap.x_range.factors = sorted_factors;
        ''')

    y_dendro_callback = CustomJS(
        args=dict(rank_select=ele["heatmap"]["wid"]["rank_select"],
                  y_groupby_select=ele["heatmap"]["wid"]["y_groupby_select"],
                  y_sort_select=ele["heatmap"]["wid"]["y_sort_select"],
                  dendroy=ele["dendroy"]["fig"],
                  cds_p_dendro_y=cds_p_dendro_y,
                  dict_d_dedro_y=dict_d_dedro_y),
        code='''
        if (y_groupby_select.value.startsWith("cluster|")){
            const key = rank_select.value+"|"+y_groupby_select.value.replace("cluster|","");
            cds_p_dendro_y.data = {"x": dict_d_dedro_y[key+"|x"],
                                   "y": dict_d_dedro_y[key+"|y"],
                                   "c": dict_d_dedro_y[key+"|c"]};
            y_sort_select.value="none";
            y_sort_select.disabled=true;
            dendroy.visible=true;
        }else{
            cds_p_dendro_y.data = {"x": [], "y": [], "c": []};
            y_sort_select.disabled=false;
            dendroy.visible=false;
        }
        ''')

    y_select_callback = CustomJS(
        args=dict(heatmap=ele["heatmap"]["fig"],
                  cds_d_samples=cds_d_samples,
                  cds_d_metadata=cds_d_metadata,
                  cds_p_metadata=cds_p_metadata,
                  cds_p_heatmap=cds_p_heatmap,
                  rank_select=ele["heatmap"]["wid"]["rank_select"],
                  y_sort_select=ele["heatmap"]["wid"]["y_sort_select"],
                  y_groupby_select=ele["heatmap"]["wid"]["y_groupby_select"],
                  dict_d_hcluster_y=dict_d_hcluster_y),
        code='''
        // selected rank
        const rank = rank_select.value;
        var annot_samples = cds_d_samples.data["index"];

        var sorted_factors = [];
        var dict_factors = {};
        if (y_groupby_select.value.startsWith("cluster|")){
            // Clustering - Get sorted elements based on rank|method|metric
            sorted_factors = dict_d_hcluster_y[rank+"|"+y_groupby_select.value.replace("cluster|","")];
            // Default dict_factors
            for(let i = 0; i < annot_samples.length; i++){
                dict_factors[annot_samples[i]] = annot_samples[i];
            }
        }else{
            // Define value from Sort by select
            var sort_col = [];
            var sort_col_type = "string";
            if (y_sort_select.value=="none"){
                sort_col = dict_d_hcluster_y["default"];
            }else if (y_sort_select.value=="samples"){
                sort_col = annot_samples;
            }else if (y_sort_select.value=="counts"){
                sort_col = cds_d_samples.data["cnt|total"];
                sort_col_type = "numeric";
            }else if (y_sort_select.value.startsWith("metadata_num|")){
                sort_col = cds_d_metadata.data[y_sort_select.value.replace("metadata_num|","")];
                sort_col_type = "numeric";
            }else if (y_sort_select.value.startsWith("metadata_cat|")){
                sort_col = cds_d_metadata.data[y_sort_select.value.replace("metadata_cat|","")];
            }

            if(y_groupby_select.value=="none"){
                sorted_factors = grimer_sort(annot_samples, sort_col, sort_col_type, false);
                // Default dict_factors
                for(let i = 0; i < annot_samples.length; i++){
                    dict_factors[annot_samples[i]] = annot_samples[i];
                }
            }else if (y_groupby_select.value.startsWith("group_metadata|")){
                const group_metadata = y_groupby_select.value.replace("group_metadata|","");

                // group entries and replace empty with space " "
                var groupby_col = cds_d_metadata.data[group_metadata].map(function(m) { return m == "" ? " " : m; });

                var factors = [];
                for(let i = 0; i < annot_samples.length; i++){
                    dict_factors[annot_samples[i]] = [groupby_col[i], annot_samples[i]];
                    factors.push([groupby_col[i], annot_samples[i]]);
                }
                sorted_factors = grimer_sort(factors, sort_col, sort_col_type, false, groupby_col);
            }
        }

        // update factors on heatmap col otherwise remove
        for (let i = 0; i < cds_p_heatmap.data["index"].length; i++) {
            if(cds_p_heatmap.data["rank"][i]==rank){
                cds_p_heatmap.data["factors_sample"][i] = dict_factors[cds_p_heatmap.data["index"][i]];
            }else{
                cds_p_heatmap.data["factors_sample"][i] = "";
            }
        }

        if (cds_p_metadata){
            for (let i = 0; i < cds_p_metadata.data["index"].length; i++) {
                cds_p_metadata.data["factors"][i] = dict_factors[cds_p_metadata.data["index"][i]];
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
            xaxis.major_label_text_font_size = "10px";
            xaxis.major_tick_line_color="black";
        }else{
            xaxis.major_label_text_font_size = "0px";
            xaxis.major_tick_line_color=null;
        }
        if(this.active.includes(1)){
            yaxis.major_label_text_font_size = "10px";
            yaxis.major_tick_line_color="black";
        }else{
            yaxis.major_label_text_font_size = "0px";
            yaxis.major_tick_line_color=null;
        }
        ''')

    ele["heatmap"]["wid"]["toggle_label"].js_on_click(toggle_label_callback)
    ele["heatmap"]["wid"]["rank_select"].js_on_change('value', y_select_callback, x_select_callback, x_dendro_callback, y_dendro_callback)
    ele["heatmap"]["wid"]["x_sort_select"].js_on_change('value', x_select_callback, x_dendro_callback)
    ele["heatmap"]["wid"]["x_groupby_select"].js_on_change('value', x_select_callback, x_dendro_callback)
    ele["heatmap"]["wid"]["y_sort_select"].js_on_change('value', y_select_callback, y_dendro_callback)
    ele["heatmap"]["wid"]["y_groupby_select"].js_on_change('value', y_select_callback, y_dendro_callback)


def link_metadata_widgets(ele, cds_p_metadata, cds_d_metadata, max_metadata_cols):
    metadata_multiselect_callback = CustomJS(
        args=dict(metadata_heatmap=ele["metadata"]["fig"],
                  metadata_heatmap_xaxis=ele["metadata"]["fig"].xaxis[0] if cds_p_metadata else None,
                  metadata_multiselect=ele["metadata"]["wid"]["metadata_multiselect"],
                  legend_colorbars=ele["metadata"]["wid"]["legend_colorbars"],
                  toggle_legend=ele["metadata"]["wid"]["toggle_legend"],
                  max_metadata_cols=max_metadata_cols,
                  cds_p_metadata=cds_p_metadata,
                  cds_d_metadata=cds_d_metadata),
        code='''
        const index_len = cds_d_metadata.data["index"].length;

        var x_factors = [];
        var empty_y_values = new Array(index_len);
        for (var i = 0; i < index_len; ++i) empty_y_values[i]=["", ""];
        // hide all legends
        for (let md_header in legend_colorbars) legend_colorbars[md_header].visible = false;

        // set legend orientation
        if(metadata_heatmap_xaxis){
            if(metadata_multiselect.value.length==1)
                metadata_heatmap_xaxis.major_label_orientation = "horizontal";
            else
                metadata_heatmap_xaxis.major_label_orientation = 0.7;
        }

        for(var s=0; s < max_metadata_cols; ++s){
            if (s<metadata_multiselect.value.length){
                var selected = metadata_multiselect.value[s];
                var y_values = new Array(index_len);
                for (var i = 0; i < index_len; ++i){
                    var val = cds_d_metadata.data[selected][i];
                    // fix conversion from float to integer (30.0).toString() = "30"
                    if (Number.isInteger(val)){
                        val = val + ".0";
                    } else {
                        val = val.toString();
                    }
                    y_values[i]=[selected, val];
                }
                cds_p_metadata.data[(s+1).toString()] = y_values;
                x_factors.push((s+1).toString());

                // show legend
                if(toggle_legend.active.includes(0))
                    legend_colorbars[selected].visible = true;

            }else{
                cds_p_metadata.data[(s+1).toString()] = empty_y_values;
            }
        }
        metadata_heatmap.x_range.factors = x_factors;
        cds_p_metadata.change.emit();
        ''')

    if cds_d_metadata:
        ele["metadata"]["wid"]["metadata_multiselect"].js_on_change('value', metadata_multiselect_callback)
        ele["metadata"]["wid"]["toggle_legend"].js_on_click(metadata_multiselect_callback)


def link_obstable_filter(ele, cds_m_obstable, active_ranks):
    filter_callback = CustomJS(
        args=dict(cds_m_obstable=cds_m_obstable,
                  active_ranks=active_ranks,
                  filter=ele["obstable"]["filter"],
                  frequency_spinner=ele["obstable"]["wid"]["frequency_spinner"],
                  counts_perc_avg_spinner=ele["obstable"]["wid"]["counts_perc_avg_spinner"],
                  total_counts_spinner=ele["obstable"]["wid"]["total_counts_spinner"],
                  name_multichoice=ele["obstable"]["wid"]["name_multichoice"],
                  ),
        code='''
        const indices = [];
        for (var i = 0; i < cds_m_obstable.length; i++) {
            if (cds_m_obstable.data['col|frequency_perc'][i] < (frequency_spinner.value/100)){
                continue;
            }
            if (cds_m_obstable.data['col|counts_perc_avg'][i] < (counts_perc_avg_spinner.value/100)){
                continue;
            }
            if (cds_m_obstable.data['col|total_counts'][i] < (total_counts_spinner.value)){
                continue;
            }
            if (name_multichoice.value.length > 0 ){
                var found = false;
                for(let r = 0; r < active_ranks.length; r++){
                    // Compare all names on multichoice (array) against cell
                    if (name_multichoice.value.indexOf(cds_m_obstable.data["tax|"+active_ranks[r]][i]) >= 0){
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
        filter.indices = indices;
        cds_m_obstable.change.emit();
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
        var factors = new Set();
        for(let i = 0; i < cds_p_correlation.data["index"].length; i++){
            if(cds_p_correlation.data["rank"][i]==this.value){
                factors.add(cds_p_correlation.data["index"][i]);
                factors.add(cds_p_correlation.data["taxid"][i]);
            }
        }
        var sorted_factors = [...factors].sort();
        correlation.y_range.factors = sorted_factors;
        var rev_sorted_factors = [...sorted_factors].reverse();
        correlation.x_range.factors = rev_sorted_factors;
        ''')

    filter_callback = CustomJS(
        args=dict(filter=ele["correlation"]["filter"],
                  neg_slider=ele["correlation"]["wid"]["neg_slider"],
                  pos_slider=ele["correlation"]["wid"]["pos_slider"],
                  cds_p_correlation=cds_p_correlation),
        code='''
        const indices = [];
        for (var i = 0; i < cds_p_correlation.data["index"].length; i++) {
            const rho = cds_p_correlation.data["rho"][i];
            if ((rho >= neg_slider.value[0] && rho <= neg_slider.value[1]) ||
                (rho >= pos_slider.value[0] && rho <= pos_slider.value[1]))
            {
                indices.push(i)
            }
        }
        filter.indices = indices;
        cds_p_correlation.change.emit();
        ''')

    toggle_label_callback = CustomJS(
        args=dict(xaxis=ele["correlation"]["fig"].xaxis[0], yaxis=ele["correlation"]["fig"].yaxis[0]),
        code='''
        if(this.active.includes(0)){
            xaxis.major_label_text_font_size = "10px";
            xaxis.major_tick_line_color="black";
            yaxis.major_label_text_font_size = "10px";
            yaxis.major_tick_line_color="black";
        }else{
            xaxis.major_label_text_font_size = "0px";
            xaxis.major_tick_line_color=null;
            yaxis.major_label_text_font_size = "0px";
            yaxis.major_tick_line_color=null;
        }
        ''')

    ele["correlation"]["wid"]["toggle_label"].js_on_click(toggle_label_callback)
    ele["correlation"]["wid"]["pos_slider"].js_on_change('value', filter_callback)
    ele["correlation"]["wid"]["neg_slider"].js_on_change('value', filter_callback)
    ele["correlation"]["wid"]["rank_select"].js_on_change('value', rank_select_callback)


def link_obsbars_widgets(ele, cds_p_obsbars, dict_d_topobs, dict_d_sampleobs, cds_d_samples, top_obs_bars, dict_d_taxname, cds_d_metadata, cds_p_sampletable):
    rank_select_callback = CustomJS(
        args=dict(sort_select=ele["obsbars"]["wid"]["sort_select"],
                  legend=ele["obsbars"]["legend"],
                  cds_p_obsbars=cds_p_obsbars,
                  dict_d_sampleobs=dict_d_sampleobs,
                  cds_d_samples=cds_d_samples,
                  dict_d_topobs=dict_d_topobs,
                  dict_d_taxname=dict_d_taxname,
                  top_obs_bars=top_obs_bars),
        code='''
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
                    let sample = cds_p_obsbars.data["index"][s];
                    let val = 0;
                    if(dict_d_sampleobs[rank][taxid][sample]!=undefined){
                        val = dict_d_sampleobs[rank][taxid][sample];
                    }
                    cds_p_obsbars.data[i.toString()][s] = (val/total[s])*100;
                    // sum counts for sample
                    sum_assigned[s]+=val;
                    // update legend label
                    legend.items[i].label = (i+1).toString() + ") " + dict_d_taxname[taxid];
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
                  cds_p_sampletable=cds_p_sampletable,
                  cds_p_obsbars=cds_p_obsbars,
                  cds_d_samples=cds_d_samples,
                  cds_d_metadata=cds_d_metadata),
        code='''

        // Factors can be: index (sort_col|index), [md1, index] or [md1, md2, index]
        var factors = [];
        var sorted_factors = [];

        // get index of selected indices
        var selected_indices = cds_p_sampletable.selected.indices;

        if(selected_indices.length){
            // samples
            var annot_samples = cds_d_samples.data["index"];

            // Define value from Sort by select
            var sort_col;
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
                      return [groupby_col1[i], m, annot_samples[i]];
                    });

                    // only selected_indices
                    sorted_factors = grimer_sort(factors, sort_col, "numeric", false, groupby_col2, groupby_col1, selected_indices);
                }else{
                    sorted_factors = grimer_sort(factors, sort_col, "numeric", false, groupby_col1, [], selected_indices);
                }

            }else{
                // Single factors, just use the sample index
                factors = annot_samples;
                sorted_factors = grimer_sort(factors, sort_col, "numeric", false, [], [], selected_indices);
            }

            // Change value of the factors on the obsbars cds
            cds_p_obsbars.data["factors"] = factors; 

        }

        // Plot sorted factors
        obsbars.x_range.factors = sorted_factors;
        ''')

    toggle_label_callback = CustomJS(
        args=dict(xaxis=ele["obsbars"]["fig"].xaxis[0]),
        code='''
        if(this.active.includes(0)){
            xaxis.major_label_text_font_size = "10px";
            xaxis.major_tick_line_color="black";
        }else{
            xaxis.major_label_text_font_size = "0px";
            xaxis.major_tick_line_color=null;
        }
        ''')

    cds_p_sampletable.selected.js_on_change('indices', sort_groupby_callback)
    ele["obsbars"]["wid"]["toggle_label"].js_on_click(toggle_label_callback)
    ele["obsbars"]["wid"]["groupby1_select"].js_on_change('value', sort_groupby_callback)
    ele["obsbars"]["wid"]["groupby2_select"].js_on_change('value', sort_groupby_callback)
    ele["obsbars"]["wid"]["sort_select"].js_on_change('value', sort_groupby_callback)
    ele["obsbars"]["wid"]["rank_select"].js_on_change('value', rank_select_callback, sort_groupby_callback)


def link_sampletable_select(ele, cds_p_sampletable, cds_d_metadata):

    select_callback = CustomJS(
        args=dict(cds_p_sampletable=cds_p_sampletable,
                  cds_d_metadata=cds_d_metadata,
                  total_counts_spinner=ele["sampletable"]["wid"]["total_counts_spinner"],
                  assigned_spinner=ele["sampletable"]["wid"]["assigned_spinner"],
                  metadata_multichoice=ele["sampletable"]["wid"]["metadata_multichoice"],
                  ),
        code='''
        var selected_indices = [];
        for (var i = 0; i < cds_p_sampletable.length; i++) {
            if (cds_p_sampletable.data['col|total'][i] < total_counts_spinner.value){
                continue;
            }
            if (cds_p_sampletable.data['col|assigned_perc'][i] < (assigned_spinner.value/100)){
                continue;
            }
            if (metadata_multichoice.value.length > 0 ){
                var found = false;
                for (var m=0; m < metadata_multichoice.value.length; ++m){
                    const md = metadata_multichoice.value[m].split("|");
                    if(cds_d_metadata.data[md[0]][i]==md[1]){
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    continue;
                }
            }
            selected_indices.push(i);
        }
        cds_p_sampletable.selected.indices = selected_indices;
    ''')
    ele["sampletable"]["wid"]["total_counts_spinner"].js_on_change('value', select_callback)
    ele["sampletable"]["wid"]["assigned_spinner"].js_on_change('value', select_callback)
    if cds_d_metadata:
        ele["sampletable"]["wid"]["metadata_multichoice"].js_on_change('value', select_callback)
