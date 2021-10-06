#General
import pandas as pd
import numpy as np
from math import pi

#Internal
from grimer.utils import print_df, transform_table, print_log, pairwise_rho

#Bokeh
from bokeh.models import ColumnDataSource


def generate_dict_taxname(tax, taxids):
    id_name = {}
    for i in taxids:
        n = tax.name(i) if tax else i
        id_name[i] = n if n else i
    return id_name


def generate_cds_plot_references(table, tax, references):
    # Stacked list of references, accounting for lineage matches
    # index -> observations (repeated)
    # columns -> "rank", "ref", "direct", "parent"
    clist = []
    for rank in table.ranks():
        for obs in table.observations(rank):
            for desc, ref in references.items():
                direct = ref.get_refs_count(obs, direct=True)
                parent = ref.get_refs_count(obs, parents=True)
                if direct + parent > 0:
                    clist.append([obs, rank, desc, direct, parent])

    df_references = pd.DataFrame(clist, columns=["obs", "rank", "ref", "direct", "parent"])
    df_references.set_index('obs', inplace=True)

    print_df(df_references, "df_references -> cds_p_references")
    return ColumnDataSource(df_references)

def generate_cds_annotations(table, references, controls, decontam):
    # Stacked matrix of true annotations (omit false)
    # index -> taxids
    # columns -> rank, annot

    df_annotations = pd.DataFrame(columns=["rank", "annot"])
    for rank in table.ranks():
        # Generate a DataFrame to use as source in tables
        df_rank = pd.DataFrame(index=table.observations(rank))

        if decontam:
            df_rank["decontam"] = decontam.get_contaminants(rank, df_rank.index)

        for desc, ref in references.items():
            df_rank[desc] = table.observations(rank).map(lambda x: ref.get_refs_count(x, direct=True)) >= 1

        if controls:
            for desc, ctrl in controls.items():
                df_rank[desc] = table.observations(rank).map(lambda x: ctrl.get_refs_count(x, direct=True)) >= 1

        df_rank = pd.DataFrame(df_rank.stack(), columns=["val"]).reset_index(1)
        df_rank.rename(columns={"level_1": "annot"}, inplace=True)
        df_rank = df_rank[df_rank["val"]]  # keep only true entries
        df_rank["rank"] = rank  # set rank

        if "val" in df_rank.columns:
            df_rank.drop(columns="val", inplace=True)  # drop boolean col

        # Concat in the main df
        df_annotations = pd.concat([df_annotations, df_rank], axis=0)

    print_df(df_annotations, "df_annotations -> cds_p_annotations")
    return ColumnDataSource(df_annotations)


def generate_cds_obstable(table, tax, references, controls, control_samples, decontam):
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
        # Average percentage of counts among all samples
        df_rank["col|counts_perc_avg"] = table.get_counts_perc_avg(rank)
        df_rank["col|total_counts"] = table.get_total_counts(rank)

        # If active - add decontam True/False results
        if decontam:
            df_rank["col|decontam"] = decontam.get_contaminants(rank, df_rank.index)

        # Add a column for each Annotation source
        for desc, ref in references.items():
            df_rank["col|" + desc] = table.observations(rank).map(lambda x: ref.get_refs_count(x, direct=True)).to_list()

        # Add a column for each Control source
        if controls:
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

    print_df(df_obstable, "df_obstable -> cds_p_obstable")
    return ColumnDataSource(df_obstable)


def generate_cds_bars(table):
    # index unique sample-ids
    # aux| auxiliary values (not plotted)
    # bar| values plotted as bars (sample counts)
    # tax| values plotted as circles (taxa value)

    df_bars = pd.DataFrame(index=table.samples)
    # factors: set the x-axis reference for plotting, it can be dinamically changed (with groups)
    df_bars["aux|factors"] = df_bars.index
    df_bars["bar|unassigned"] = table.unassigned
    # Initialized with "Assigned" of first rank
    df_bars["bar|selected"] = table.data[table.ranks()[0]].sum(axis=1)
    # Total assigned - assigned to rank
    df_bars["bar|others"] = (table.total - table.unassigned) - df_bars["bar|selected"]
    # Add empty cols for taxa values, to be dynamically inserted (None to avoid printing 0)
    for rank in table.ranks():
        df_bars["tax|" + rank] = None

    print_df(df_bars, "df_bars -> cds_p_samplebars")
    return ColumnDataSource(df_bars)


def generate_cds_samples(table, references, controls, decontam):
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
    source_list = [references.items()]
    if controls:
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

    print_df(df_samples, "df_samples -> cds_d_samples")
    return ColumnDataSource(df_samples)


def generate_cds_metadata(metadata):
    # index -> sample-ids
    # columns -> metadata fields
    # values -> metadata values
    df_md = metadata.get_data()
    print_df(df_md, "df_md -> cds_d_metadata")
    return ColumnDataSource(df_md)


def generate_cds_plot_metadata(metadata, max_metadata_cols):
    # index (unique sample-ids)
    # md0, md1, ..., md(max_metadata_cols)
    # values (metadata field, metadata values)

    df_plot_md = pd.DataFrame(index=metadata.data.index)
    # Fill with first N metadata fields
    metadata_fields = metadata.get_col_headers().to_list()
    for i in range(max_metadata_cols):
        # Same transformation done in the colormap for numeric entries
        df_plot_md["md" + str(i)] = [(metadata_fields[i], '{:.16g}'.format(md_value) if not isinstance(md_value, str) else md_value) for md_value in metadata.get_col(metadata_fields[i])]

    print_df(df_plot_md, "df_plot_md -> cds_p_metadata")
    return ColumnDataSource(df_plot_md)


def generate_cds_plot_decontam(decontam):
    # index unique sample-ids
    # concentrations from decontam inputs
    # controls from decontam inputs
    # counts: field to be dynamically filled with click on obstable
    df_decontam = decontam.get_data()
    df_decontam["controls"] = df_decontam["controls"].map({True: 'Control', False: 'Sample'})
    df_decontam["counts"] = None
    print_df(df_decontam, "df_decontam -> cds_p_decontam")
    return ColumnDataSource(df_decontam)


def generate_cds_decontam(decontam, ranks):
    """
    cds based on a dict with valid values to plot model lines
    {taxid: (contam_y1, contam_y2, non_contam_y, pval)}
    """
    dict_coord_mod = {}
    for rank in ranks:
        df_valid_vals = decontam.rank[rank].dropna(subset=['contam'])
        pval = decontam.get_pvalues(rank, df_valid_vals.index)
        vals = list(zip(df_valid_vals["contam"], df_valid_vals["contam_2"], df_valid_vals["non.contam"], pval))
        dict_coord_mod.update(dict(zip(df_valid_vals.index, vals)))

    print_df(dict_coord_mod, "dict_coord_mod -> cds_d_decontam_models")
    return ColumnDataSource(dict_coord_mod)


def generate_cds_plot_decontam_models(decontam):
    """
    cds based on a dict with 3 pairs of values to plot. x is shared among y_cont and y_noncont
    # {x: [min,max], y_cont: [None,None], y_noncont: [None,None]}
    """
    dict_decontam_models = {}
    dict_decontam_models["x"] = [decontam.get_data()["concentration"].min(),
                                 decontam.get_data()["concentration"].max()]
    dict_decontam_models["y_cont"] = [None, None]
    dict_decontam_models["y_noncont"] = [None, None]
    print_df(dict_decontam_models, "dict_decontam_models -> cds_p_decontam_models")
    return ColumnDataSource(dict_decontam_models)


def generate_cds_sampleobs(table):
    # matrix-like cds with raw counts
    # index -> sample-ids
    # columns -> taxids (from all ranks)
    # values are observation raw counts
    df_sampleobs = pd.DataFrame(index=table.samples)
    for rank in table.ranks():
        df_sampleobs = pd.concat([df_sampleobs, table.data[rank]], axis=1)

    # fill NaN with zero so bars do not "dissapear" when plotting
    df_sampleobs.fillna(0, inplace=True)
    print_df(df_sampleobs, "df_sampleobs -> cds_d_sampleobs")
    return ColumnDataSource(df_sampleobs)


def generate_cds_heatmap(table, transformation, replace_zero_value, show_zeros):
    # Stacked matrix of raw counts + transformed value
    # index -> sample-ids (repeated)
    # obs
    # rank
    # ov -> original value (raw counts)
    # tv -> transformed values (user choice: log10, clr, ...)

    df_heatmap = pd.DataFrame(columns=["obs", "rank", "ov", "tv"])
    for rank in table.ranks():
        stacked_rank_df = pd.DataFrame(table.data[rank].stack(), columns=["ov"]).reset_index(1)
        # Rename first col to obs
        stacked_rank_df.rename(columns={stacked_rank_df.columns[0]: "obs"}, inplace=True)
        stacked_rank_df["rank"] = rank
        tv = transform_table(table.data[rank], table.total, transformation, replace_zero_value)
        stacked_rank_df["tv"] = tv.stack().values
        #Drop zeros based on original counts
        if not show_zeros:
            stacked_rank_df = stacked_rank_df[stacked_rank_df["ov"] > 0]

        df_heatmap = pd.concat([df_heatmap, stacked_rank_df], axis=0)

    print_df(df_heatmap, "df_heatmap -> cds_p_heatmap")
    return ColumnDataSource(df_heatmap)


def generate_dict_hcluster(table, hcluster):
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

    print_df(leaves_x, "leaves_x -> dict_d_hcluster_x")
    print_df(leaves_y, "leaves_y -> dict_d_hcluster_y")
    return leaves_x, leaves_y


def generate_cds_plot_dendro():
    # Empty CDS {"x": [], "y": [], "c": []}
    dendro_x = {"x": [], "y": [], "c": []}
    dendro_y = {"x": [], "y": [], "c": []}
    print_df(dendro_x, "dendro_x -> cds_p_dendro_x")
    print_df(dendro_y, "dendro_y -> cds_p_dendro_y")
    return ColumnDataSource(dendro_x), ColumnDataSource(dendro_y)


def generate_dict_dendro(table, dendro):
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


def generate_dict_topobs(table, top_obs_bars):
    dict_top_taxa = {}
    for rank in table.ranks():
        dict_top_taxa[rank] = table.get_top(rank, top_obs_bars)
    print_df(dict_top_taxa, "dict_top_taxa -> dict_d_topobs")
    return dict_top_taxa


def generate_dict_refs(table, references):
    # dict with information about sources and references
    # references can be repeated among descriptions, sources and taxids
    # {taxid: {source: {desc: [refs]}}
    d_refs = {}
    # Get only valid taxids
    used_ids = set()
    for rank in table.ranks():
        used_ids.update(table.observations(rank))

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

    print_df(d_refs, "d_refs -> dict_d_refs")
    return d_refs


def generate_cds_correlation(table, top_obs_corr, replace_zero_value):
    # index (repeated taxids)
    # other taxid
    # rank
    # rho

    df_corr = pd.DataFrame(columns=["taxid", "rank", "rho"])
    for rank in table.ranks():
        if top_obs_corr:
            top_taxids = table.get_top(rank, top_obs_corr)
            matrix = table.get_subtable(taxids=top_taxids, rank=rank)
        else:
            top_taxids = sorted(table.observations(rank))
            matrix = table.data[rank]

        # No correlation with just one observation
        if len(matrix.columns) >= 2:

            rho = pairwise_rho(transform_table(matrix, 0, "clr", replace_zero_value).values)

            if len(matrix.columns) == 2:
                # If there are only 2 observations, return in a float
                # re-format in a matrix shape
                rho = np.array([[np.nan, np.nan], [rho, np.nan]])
            else:
                # fill upper triangular matrix (mirrored values) with nan to be ignored by pandas
                # to save half of the space
                rho[np.triu_indices(rho.shape[0])] = np.nan

            stacked_rank_df = pd.DataFrame(rho, index=top_taxids, columns=top_taxids).stack(dropna=False).reset_index(1)
            stacked_rank_df.rename(columns={"level_1": "taxid"}, inplace=True)
            stacked_rank_df.rename(columns={0: "rho"}, inplace=True)
            stacked_rank_df["rank"] = rank

            # Drop NA for rho (missing values and upper triangular matrix)
            stacked_rank_df.dropna(subset=['rho'], inplace=True)

            df_corr = pd.concat([df_corr, stacked_rank_df], axis=0)

    print_df(df_corr, "df_corr -> cds_p_correlation")
    return ColumnDataSource(df_corr)


def generate_cds_obsbars(table, top_obs_bars):
    # index (unique sample-ids)
    # cols: 0, 1, ..., top_obs_bars, unassigned, others, factors

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

    print_df(df_obsbars, "df_obsbars -> cds_p_obsbars")
    return ColumnDataSource(df_obsbars)


def generate_cds_mgnify(mgnify, table, tax):
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

    print_df(df_mgnify, "df_mgnify -> cds_p_mgnify")
    return ColumnDataSource(df_mgnify)
