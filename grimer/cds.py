#General
import pandas as pd
import numpy as np
from math import pi

#Internal
from grimer.func import print_df, transform_table, print_log, pairwise_rho, format_js_toString

#Bokeh
from bokeh.models import ColumnDataSource


def dict_taxname(tax, taxids):
    """
    mapping taxids to names
    (or names to names if taxid is not used)
    """
    id_name = {}
    for i in taxids:
        n = tax.name(i) if tax else i
        id_name[i] = n if n else i
    return id_name


def cds_plot_references(table, tax, references):
    # Stacked list of references, accounting for lineage matches
    # index -> observations (repeated)
    # columns -> "rank", "ref", "direct", "parent"
    clist = []
    if references is not None:
        for rank in table.ranks():
            for obs in table.observations(rank):
                for desc, ref in references.items():
                    direct = ref.get_refs_count(obs, direct=True)
                    parent = ref.get_refs_count(obs, parents=True)
                    if direct + parent > 0:
                        clist.append([obs, rank, desc, direct, parent])

    df_references = pd.DataFrame(clist, columns=["obs", "rank", "ref", "direct", "parent"])
    df_references.set_index('obs', inplace=True)

    print_df(df_references, "cds_p_references")
    return ColumnDataSource(df_references)


def cds_annotations(table, references, controls, decontam, control_samples):
    # Stacked matrix of true annotations (omit false)
    # index -> taxids
    # columns -> rank, annot

    df_annotations = pd.DataFrame(columns=["rank", "annot", "factors", "ov", "tv"])
    for i, rank in enumerate(table.ranks()):
        # Generate a DataFrame to use as source in tables
        df_rank = pd.DataFrame(index=table.observations(rank))

        if decontam:
            contaminants = decontam.get_contaminants(rank, df_rank.index).values
            if contaminants.any():
                df_rank["decontam"] = decontam.get_pscore(rank, df_rank.index)[contaminants]

        if references is not None:
            for desc, ref in references.items():
                df_rank[desc] = table.observations(rank).map(lambda x: ref.get_refs_count(x, direct=True))
                df_rank.loc[df_rank[desc] == 0, desc] = np.nan

        if controls is not None:
            for desc, ctrl in controls.items():
                control_table = table.get_subtable(samples=control_samples[desc], rank=rank)
                freq_perc_control = control_table.gt(0).sum(axis=0) / control_table.shape[0]
                df_rank[desc] = table.observations(rank).map(freq_perc_control).to_list()

        df_rank = pd.DataFrame(df_rank.stack(), columns=["ov"]).reset_index(1)
        df_rank.rename(columns={"level_1": "annot"}, inplace=True)

        # add transformed values to fit same scale on heatmap
        # Decontam reverse p-score normalized
        if not df_rank[df_rank["annot"] == "decontam"].empty:
            min_val = df_rank[df_rank["annot"] == "decontam"]["ov"].min()
            max_val = df_rank[df_rank["annot"] == "decontam"]["ov"].max()
            df_rank.loc[df_rank["annot"] == "decontam", "tv"] = 1 - ((df_rank[df_rank["annot"] == "decontam"]["ov"] - min_val) / (max_val - min_val))

        # max references divided by max
        if references is not None:
            for desc, ref in references.items():
                if not df_rank[df_rank["annot"] == desc].empty:
                    max_val = df_rank[df_rank["annot"] == desc]["ov"].max()
                    df_rank.loc[df_rank["annot"] == desc, "tv"] = df_rank.loc[df_rank["annot"] == desc, "ov"] / max_val

        # keep same percentage
        if controls is not None:
            for desc, ctrl in controls.items():
                if not df_rank.loc[df_rank["annot"] == desc].empty:
                    df_rank.loc[df_rank["annot"] == desc, "tv"] = df_rank.loc[df_rank["annot"] == desc, "ov"]

        df_rank["rank"] = rank  # set rank
        df_rank["factors"] = df_rank.index if i == 0 else ""  # initialize just for first rank (save space)

        # Concat in the main df
        df_annotations = pd.concat([df_annotations, df_rank], axis=0)

    print_df(df_annotations, "cds_p_annotations")
    return ColumnDataSource(df_annotations)


def cds_obstable(table, tax, references, controls, control_samples, decontam):
    # index unique taxids
    # col|...  values to plot to columns in the datatable
    # tax|...  auxiliary lineage of taxa entries
    # aux|ref  auxiliary references identifiers

    df_obstable = pd.DataFrame()
    # Create unified DataFrame with all ranks used
    for rank in table.ranks():
        # Generate a DataFrame to use as source in tables
        df_rank = pd.DataFrame(index=table.observations(rank))
        df_rank["col|rank"] = rank
        if tax:
            df_rank["col|name"] = table.observations(rank).map(lambda txid: tax.name(txid) if tax.name(txid) else txid).to_list()
        else:
            df_rank["col|name"] = table.observations(rank)

        # Frequency of taxa among all samples
        df_rank["col|frequency_perc"] = table.get_frequency_perc(rank)
        df_rank["col|counts_perc_avg"] = table.get_counts_perc_avg_samples(rank)
        # Average percentage of counts among all samples
        df_rank["col|total_counts"] = table.get_counts(rank)

        # If active - add decontam True/False results
        if decontam:
            df_rank["col|decontam"] = decontam.get_contaminants(rank, df_rank.index)

        # Add a column for each Annotation source
        if references is not None:
            for desc, ref in references.items():
                df_rank["col|" + desc] = table.observations(rank).map(lambda x: ref.get_refs_count(x, direct=True)).to_list()

        # Add a column for each Control source
        if controls is not None:
            # calculate frequency for each group of control provided
            for desc, ctrl in controls.items():
                control_table = table.get_subtable(samples=control_samples[desc], rank=rank)
                freq_perc_control = control_table.gt(0).sum(axis=0) / control_table.shape[0]
                df_rank["col|" + desc] = table.observations(rank).map(freq_perc_control).fillna(0).to_list()

        # Add col for each rank with parent taxid if exists, linking entries in their lineage for filtering and plotting
        for other_rank in table.ranks():
            if table.ranks().index(other_rank) > table.ranks().index(rank):
                df_rank["tax|" + other_rank] = ""
            elif other_rank != rank:
                df_rank["tax|" + other_rank] = table.observations(rank).map(lambda txid: table.get_lineage(txid, rank, other_rank)).fillna("")
            else:
                df_rank["tax|" + other_rank] = df_rank.index
        # Sort values by frequency to show on table
        df_rank.sort_values(by="col|frequency_perc", ascending=False, inplace=True)

        # Concat in the main df
        df_obstable = pd.concat([df_obstable, df_rank], axis=0)

    print_df(df_obstable, "cds_m_obstable")
    return ColumnDataSource(df_obstable)


def cds_sampletable(table):
    # index unique sample-ids
    # col|...  values to plot to columns in the datatable

    df_sampletable = pd.DataFrame(index=table.samples)
    df_sampletable["col|total"] = table.get_total() if not table.normalized else 0
    df_sampletable["col|assigned"] = table.get_assigned() if not table.normalized else 0
    df_sampletable["col|assigned_perc"] = table.get_assigned_perc()
    df_sampletable["col|unassigned"] = table.get_unassigned() if not table.normalized else 0
    df_sampletable["col|unassigned_perc"] = table.get_unassigned_perc()

    # assigned by rank
    for rank in table.ranks():
        df_sampletable["col|" + rank] = table.data[rank].sum(axis=1).divide(table.total, axis=0)

    df_sampletable.fillna(0, inplace=True)

    print_df(df_sampletable, "cds_p_sampletable")
    return ColumnDataSource(df_sampletable)


def cds_samplebars(table):
    # index unique sample-ids
    # aux| auxiliary values (not plotted)
    # bar| values plotted as bars (sample counts)
    # tax| values plotted as circles (taxa value)

    df_bars = pd.DataFrame(index=table.samples)
    # factors: set the x-axis reference for plotting, it can be dinamically changed (with groups)
    df_bars["aux|factors"] = df_bars.index
    df_bars["bar|unassigned"] = table.get_unassigned()
    # Initialized with "Assigned" of first rank
    df_bars["bar|selected"] = table.get_subtable(table.ranks()[0]).sum(axis=1)
    # Total assigned - assigned to rank
    df_bars["bar|others"] = (table.get_total() - table.get_unassigned()) - df_bars["bar|selected"]
    # Add empty cols for taxa values, to be dynamically inserted (None to avoid printing 0)
    for rank in table.ranks():
        df_bars["tax|" + rank] = None

    print_df(df_bars, "cds_p_samplebars")
    return ColumnDataSource(df_bars)


def cds_samples(table, references, controls, decontam):
    # index unique sample-ids
    # aux| auxiliary values (not plotted)
    # cnt| count values to be copied/traansformed to bars

    df_samples = pd.DataFrame(index=table.samples)
    # index to retrieve default input order
    df_samples["aux|input_order"] = range(df_samples.shape[0], 0, -1)
    df_samples["cnt|total"] = table.total
    df_samples["cnt|unassigned"] = table.unassigned

    # Keep total number of assignemnts for calculations
    df_samples["cnt|assigned"] = table.total - table.unassigned

    # Add specific rank assignements
    for rank in table.ranks():
        df_samples["cnt|" + rank + "|assigned"] = table.data[rank].sum(axis=1)

    # Add counts specific to sources
    source_list = []
    if references is not None:
        source_list.append(references.items())
    if controls is not None:
        source_list.append(controls.items())

    for sources in source_list:
        for desc, src in sources:
            for rank in table.ranks():
                idx = table.observations(rank).map(lambda x: src.get_refs_count(x, direct=True)) >= 1
                df_samples["cnt|" + rank + "|" + desc] = table.data[rank][table.observations(rank)[idx]].sum(axis=1)

    if decontam:
        contaminants = decontam.get_contaminant_list()
        for rank in table.ranks():
            idx = table.observations(rank).isin(contaminants)
            df_samples["cnt|" + rank + "|decontam"] = table.data[rank][table.observations(rank)[idx]].sum(axis=1)

    # fill NaN with zero so bars do not "dissapear" when plotting
    df_samples.fillna(0, inplace=True)

    print_df(df_samples, "cds_d_samples")
    return ColumnDataSource(df_samples)


def cds_metadata(metadata):
    # index -> sample-ids
    # columns -> metadata fields
    # values -> metadata values
    df_md = metadata.get_data()
    print_df(df_md, "cds_d_metadata")
    return ColumnDataSource(df_md)


def cds_plot_metadata(metadata, max_metadata_cols):
    # index (unique sample-ids)
    # md0, md1, ..., md(max_metadata_cols)
    # values (metadata field, metadata values)

    df_plot_md = pd.DataFrame(index=metadata.data.index, columns=["factors"] + [str(i) for i in range(1, max_metadata_cols + 1)])
    df_plot_md["factors"] = df_plot_md.index
    # Fill in only first metadata field
    first_field = metadata.get_col_headers()[0]

    df_plot_md["1"] = [(first_field, format_js_toString(md_value)) for md_value in metadata.get_col(first_field)]

    print_df(df_plot_md, "cds_p_metadata")
    return ColumnDataSource(df_plot_md)


def cds_plot_decontam(decontam):
    # index unique sample-ids
    # concentrations from decontam inputs
    # controls from decontam inputs
    # counts: field to be dynamically filled with click on obstable
    df_decontam = decontam.get_data()
    df_decontam["controls"] = df_decontam["controls"].map({True: 'Control', False: 'Sample'})
    df_decontam["counts"] = None
    print_df(df_decontam, "cds_p_decontam")
    return ColumnDataSource(df_decontam)


def cds_decontam(decontam, ranks):
    """
    cds based on a dict with valid values to plot model lines
    {taxid: (contam_y1, contam_y2, non_contam_y, pval)}
    """
    dict_coord_mod = {}
    for rank in ranks:
        df_valid_vals = decontam.rank[rank].dropna(subset=['contam'])
        pval = decontam.get_pscore(rank, df_valid_vals.index)
        vals = list(zip(df_valid_vals["contam"], df_valid_vals["contam_2"], df_valid_vals["non.contam"], pval))
        dict_coord_mod.update(dict(zip(df_valid_vals.index, vals)))

    print_df(dict_coord_mod, "cds_d_decontam_models")
    return ColumnDataSource(dict_coord_mod)


def cds_plot_decontam_models(decontam):
    """
    cds based on a dict with 3 pairs of values to plot. x is shared among y_cont and y_noncont
    # {x: [min,max], y_cont: [None,None], y_noncont: [None,None]}
    """
    dict_decontam_models = {}
    dict_decontam_models["x"] = [decontam.get_data()["concentration"].min(),
                                 decontam.get_data()["concentration"].max()]
    dict_decontam_models["y_cont"] = [None, None]
    dict_decontam_models["y_noncont"] = [None, None]
    print_df(dict_decontam_models, "cds_p_decontam_models")
    return ColumnDataSource(dict_decontam_models)


def dict_sampleobs(table):
    # dict with raw counts (not storing zeros)
    # dict_sampleobs[rank][obs][sample] = count
    dict_sampleobs = {}
    for rank in table.ranks():
        dict_sampleobs[rank] = {}
        for obs, sample_val in table.data[rank].to_dict().items():
            dict_sampleobs[rank][obs] = {}
            for sample, val in sample_val.items():
                if val > 0:
                    dict_sampleobs[rank][obs][sample] = val

    print_df(dict_sampleobs, "dict_d_sampleobs")
    return dict_sampleobs


def cds_heatmap(table, transformation, show_zeros):
    # Stacked matrix of raw counts + transformed value
    # index -> sample-ids (repeated)
    # obs
    # rank
    # ov -> original value (raw counts)
    # tv -> transformed values (user choice: log10, clr, ...)

    df_heatmap = pd.DataFrame(columns=["obs", "rank", "ov", "tv", "factors_sample", "factors_obs"])
    for i, rank in enumerate(table.ranks()):
        stacked_rank_df = pd.DataFrame(table.data[rank].stack(), columns=["ov"]).reset_index(1)
        # Rename first col to obs
        stacked_rank_df.rename(columns={stacked_rank_df.columns[0]: "obs"}, inplace=True)
        stacked_rank_df["rank"] = rank
        tv = transform_table(table.data[rank], table.total, transformation, table.zerorep)
        stacked_rank_df["tv"] = tv.stack().values
        #Drop zeros based on original counts
        if not show_zeros:
            stacked_rank_df = stacked_rank_df[stacked_rank_df["ov"] > 0]
        # initialize factors only for first rank
        #stacked_rank_df["factors_sample"] = stacked_rank_df.index
        #stacked_rank_df["factors_obs"] = stacked_rank_df["obs"]
        stacked_rank_df["factors_sample"] = stacked_rank_df.index if i==0 else ""
        stacked_rank_df["factors_obs"] = stacked_rank_df["obs"] if i==0 else ""

        df_heatmap = pd.concat([df_heatmap, stacked_rank_df], axis=0)

    df_heatmap.drop('ov', axis=1, inplace=True)
    print_df(df_heatmap, "cds_p_heatmap")
    return ColumnDataSource(df_heatmap)


def dict_hcluster(table, hcluster):
    # keys -> combination of hclusters
    # values -> sorted sample-ids

    leaves_x = {}
    # default order
    leaves_y = {"default": table.samples.to_list()}

    for rank in hcluster:
        # default order for each rank
        leaves_x["default|" + rank] = table.observations(rank).to_list()
        for method in hcluster[rank]:
            for metric in hcluster[rank][method]:
                # key
                key = rank + "|" + method + "|" + metric
                # samples
                leaves_y[key] = hcluster[rank][method][metric]["y"]["index"]
                # taxa
                leaves_x[key] = hcluster[rank][method][metric]["x"]["index"]

    print_df(leaves_x, "dict_d_hcluster_x")
    print_df(leaves_y, "dict_d_hcluster_y")
    return leaves_x, leaves_y


def cds_plot_dendro():
    # Empty CDS {"x": [], "y": [], "c": []}
    dendro_x = {"x": [], "y": [], "c": []}
    dendro_y = {"x": [], "y": [], "c": []}
    print_df(dendro_x, "cds_p_dendro_x")
    print_df(dendro_y, "cds_p_dendro_y")
    return ColumnDataSource(dendro_x), ColumnDataSource(dendro_y)


def dict_dendro(table, dendro):
    # dict_d_dedro_x and dict_d_dedro_y:
    # key -> key + "|x" , key + "|y" , key + "|c"
    # value -> list of lists (x and y) or list (c)
    dict_d_dedro_y = {}
    dict_d_dedro_x = {}

    for rank in dendro:
        for method in dendro[rank]:
            for metric in dendro[rank][method]:
                # key
                key = rank + "|" + method + "|" + metric
                # dendrogram values
                dict_d_dedro_y[key + "|x"] = dendro[rank][method][metric]["y"]["xs"]
                dict_d_dedro_y[key + "|y"] = dendro[rank][method][metric]["y"]["ys"]
                dict_d_dedro_y[key + "|c"] = dendro[rank][method][metric]["y"]["colors"]
                dict_d_dedro_x[key + "|x"] = dendro[rank][method][metric]["x"]["xs"]
                dict_d_dedro_x[key + "|y"] = dendro[rank][method][metric]["x"]["ys"]
                dict_d_dedro_x[key + "|c"] = dendro[rank][method][metric]["x"]["colors"]

    return dict_d_dedro_x, dict_d_dedro_y


def dict_topobs(table, top_obs_bars):
    dict_top_taxa = {}
    for rank in table.ranks():
        dict_top_taxa[rank] = table.get_top(rank, top_obs_bars)
    print_df(dict_top_taxa, "dict_d_topobs")
    return dict_top_taxa


def dict_refs(table, references):
    # dict with information about sources and references
    # references can be repeated among descriptions, sources and taxids
    # {taxid: {source: {desc: [refs]}}
    d_refs = {}
    # Get only valid taxids
    used_ids = set()
    for rank in table.ranks():
        used_ids.update(table.observations(rank))

    if references is not None:
        for i in used_ids:
            for sname, s in references.items():
                for ref, descs in s.get_refs_desc(i, direct=True).items():
                    for desc in descs:
                        # Only add items if they have a reference to it
                        if i not in d_refs:
                            d_refs[i] = {}
                        if sname not in d_refs[i]:
                            d_refs[i][sname] = {}
                        if desc not in d_refs[i][sname]:
                            d_refs[i][sname][desc] = []
                        d_refs[i][sname][desc].append(ref)

    print_df(d_refs, "dict_d_refs")
    return d_refs


def cds_correlation(table, corr):
    df_corr = pd.DataFrame(columns=["taxid", "rank", "rho"])
    for rank in table.ranks():
        stacked_rank_df = pd.DataFrame(corr[rank]["rho"], index=corr[rank]["observations"], columns=corr[rank]["observations"]).stack(dropna=False).reset_index(1)
        stacked_rank_df.rename(columns={"level_1": "taxid"}, inplace=True)
        stacked_rank_df.rename(columns={0: "rho"}, inplace=True)
        stacked_rank_df["rank"] = rank

        # Drop NA for rho (missing values and upper triangular matrix)
        stacked_rank_df.dropna(subset=['rho'], inplace=True)

        df_corr = pd.concat([df_corr, stacked_rank_df], axis=0)

    print_df(df_corr, "cds_p_correlation")
    return ColumnDataSource(df_corr)


def cds_obsbars(table, top_obs_bars):
    # index (unique sample-ids)
    # cols: 1, 2, ..., top_obs_bars, unassigned, others, factors

    #Load with data from first rank
    top_taxids = table.get_top(table.ranks()[0], top_obs_bars)
    df_obsbars = table.get_subtable(taxids=top_taxids, rank=table.ranks()[0], keep_shape=True)
    df_obsbars.rename(columns={c: str(i) for i, c in enumerate(df_obsbars.columns)}, inplace=True)
    # Complete table with None values
    ncol = len(df_obsbars.columns)
    while ncol < top_obs_bars:
        df_obsbars[str(ncol)] = 0
        ncol += 1
    # Other account for filtered taxa (not on top) and left over percentage for the rank without assignment
    df_obsbars["others"] = table.total - table.unassigned - df_obsbars.sum(axis=1)
    df_obsbars["unassigned"] = table.unassigned
    df_obsbars = transform_table(df_obsbars, table.total, "norm", 0) * 100
    df_obsbars["factors"] = df_obsbars.index.to_list()

    print_df(df_obsbars, "cds_p_obsbars")
    return ColumnDataSource(df_obsbars)


def cds_mgnify(mgnify, table, tax):
    # index (taxa, level, lineage)
    # count for each combination of index

    df_mgnify = pd.DataFrame(columns=["taxa", "level", "lineage", "count", "angle"])

    # Match uids (taxid or names) from input and keep only found elements
    uids = [txid for rank in table.ranks() for txid in table.observations(rank)]
    df_tmp = mgnify.data[mgnify.data['taxa'].isin(uids)]

    # reset index to properly concate later with biome lineages
    df_tmp.reset_index(drop=True, inplace=True)

    if df_tmp.empty:
        print_log("could not find matching entries on MGnify")
        return None

    # Split biome lineage
    biome_levels = df_tmp['biome'].str.split(':', expand=True)
    n_levels = biome_levels.shape[1]

    # Rename levels with full lineage, starting from second level
    biome_lineage = pd.DataFrame(biome_levels[1])
    for l in range(2, n_levels):
        biome_lineage[l] = pd.Series(biome_levels[[i for i in range(1, l + 1)]].values.tolist()).str.join(':')

    # Concat back
    df_tmp = pd.concat([biome_lineage, df_tmp], axis=1)

    # for each biome level (ignoring root 0)
    for l in range(1, n_levels):
        # group counts by biome, and fix fields
        df_biome = df_tmp.groupby(["taxa", l]).sum()
        df_biome["level"] = str(l)
        df_biome.reset_index(inplace=True)
        df_biome.rename(columns={l: "lineage"}, inplace=True)

        # Calculate angle for each taxa/level for wedges
        total_taxa_level = df_biome.groupby("taxa").sum().to_dict()["count"]
        df_biome["angle"] = (df_biome['count'] / df_biome['taxa'].map(total_taxa_level)) * (2*pi)

        # Group to the final df
        df_mgnify = pd.concat([df_mgnify, df_biome], axis=0, ignore_index=True)

    # set index
    df_mgnify.set_index('taxa', inplace=True)

    print_df(df_mgnify, "cds_p_mgnify")
    return ColumnDataSource(df_mgnify)
