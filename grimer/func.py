#General
import numpy as np
import os
import sys
import subprocess
import shlex
import pandas as pd
from pandas.api.types import is_numeric_dtype
import yaml

#Internal
from grimer.config import Config
from grimer.decontam import Decontam
from grimer.metadata import Metadata
from grimer.reference import Reference
from grimer.table import Table

# Bokeh
from bokeh.palettes import Category10, Category20, Colorblind, linear_palette, Turbo256

# MultiTax
from multitax import *

#biom
import biom

# scikit-bio
from skbio.stats.composition import clr

# Scipy
import scipy.cluster.hierarchy as sch


def parse_config_file(config):
    cfg = None
    if config:
        try:
            with open(config, 'r') as file:
                cfg = yaml.safe_load(file)
        except Exception as e:
            print_log("Failed loading configuration file [" + config + "], skipping")
            print_log(str(e))
    else:
        print_log("Not provided, skipping")
    return cfg


def parse_taxonomy(taxonomy, tax_files):
    tax = None
    if taxonomy is not None:
        try:
            if not tax_files:
                print_log("Downloading taxonomy")
            if taxonomy == "ncbi":
                tax = NcbiTx(files=tax_files, extended_names=True)
            elif taxonomy == "gtdb":
                tax = GtdbTx(files=tax_files)
            elif taxonomy == "silva":
                tax = SilvaTx(files=tax_files)
            elif taxonomy == "greengenes":
                tax = GreengenesTx(files=tax_files)
            elif taxonomy == "ott":
                tax = OttTx(files=tax_files, extended_names=True)
            else:
                raise
        except Exception as e:
            print_log("Failed loading " + taxonomy + " taxonomy, skipping")
            print_log(str(e))
    else:
        print_log("Not provided, skipping")
    return tax


def parse_table(args, tax):
    # Specific default params if biom file is provided
    if args.input_file.endswith(".biom"):
        if not args.level_separator:
            args.level_separator = ";"
        args.transpose = True

    # Read and return full table with separated total and unassigned counts (sharing same index)
    table_df, total, unassigned = parse_input_file(args.input_file, args.unassigned_header, args.transpose, args.sample_replace, args.cumm_levels)

    if table_df.empty:
        raise Exception("Error parsing input file")

    # Define if table is already normalized (0-100) or has count data
    if args.values == "count":
        normalized = False
    elif args.values == "normalized":
        normalized = True
    elif (table_df.sum(axis=1).round() == 100).all() or (table_df % 1 != 0).any().any():
        normalized = True
    else:
        normalized = False

    # Zero replacement
    try:
        replace_zero_value = table_df[table_df.gt(0)].min().min() / int(args.replace_zeros)
    except:
        replace_zero_value = float(args.replace_zeros)
    if replace_zero_value == 1 and args.transformation == "log":
        replace_zero_value = 0.999999  # Do not allow value 1 using log

    # Split table into ranks. Ranks are either in the headers in multi level tables or will be created for a one level table
    if args.level_separator:
        ranked_tables, lineage = parse_multi_table(table_df, args.ranks, tax, args.level_separator, args.obs_replace, args.cumm_levels)
    else:
        ranked_tables, lineage = parse_single_table(table_df, args.ranks, tax, Config.default_rank_name)

    if not ranked_tables:
        raise Exception("Error parsing input file")

    table = Table(table_df.index, total, unassigned, lineage, normalized, replace_zero_value)

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
        filtered_trimmed_t = trim_table(filter_input_table(t, total, args.min_frequency, args.max_frequency, args.min_count, args.max_count, normalized))
        if t.empty:
            print_log("No valid entries, skipping")
        else:
            # Trim table for empty zeros rows/cols
            table.add_rank(r, filtered_trimmed_t)
            print_log("Total valid observations: " + str(len(table.observations(r))))

    print_log("")

    if not normalized:
        print_log("Total assigned (counts): " + str(table.get_total().sum() - table.get_unassigned().sum()))
        print_log("Total unassigned (counts): " + str(table.get_unassigned().sum()))
        print_log("")

    return table


def parse_metadata(args, samples):
    metadata = None

    # Parse metadata as DataFrame (md)
    md = pd.DataFrame()
    if args.metadata_file:
        # Parse table as dataframe
        md = pd.read_table(args.metadata_file, sep='\t', header=0, skiprows=0, index_col=0, dtype={0: str})
    elif args.input_file.endswith(".biom"):
        try:
            biom_in = biom.load_table(args.input_file)
            if biom_in.metadata() is not None:
                md = biom_in.metadata_to_dataframe(axis="sample")
        except:
            print_log("Error parsing metadata from BIOM file, skipping")
            return None

    if md.empty:
        print_log("No valid metadata, skipping")
        return None

    # Enforce string index
    md.index = md.index.astype('str')

    # Return type of columns, remove metadata row if present from metadata
    md_types = define_metadata_types(md)

    # types defined on file
    if str(md.index[0]).startswith("#"):
        # Drop row with types from main data
        md.drop(md_types.name, inplace=True)
        # Enforce column type on dataframe
        md[md_types[md_types == "categorical"].index] = md[md_types[md_types == "categorical"].index].astype(str)
        md[md_types[md_types == "numeric"].index] = md[md_types[md_types == "numeric"].index].apply(pd.to_numeric)

    # Convert datatypes to adequate numeric values (int, float)
    md = md.convert_dtypes(infer_objects=False, convert_string=False, convert_boolean=False)
    # Re-convert everything to object to standardize (int64 NA is not seriazable on bokeh)
    md = md.astype("object")

    # Remove empty fields
    null_cols = md.isna().all(axis=0)
    if any(null_cols):
        md = md.loc[:, ~null_cols]
        md_types = md_types[~null_cols]
        print_log(str(sum(null_cols)) + " metadata fields removed without valid values")

    # Convert NaN on categorical to ""
    md[md_types[md_types == "categorical"].index] = md[md_types[md_types == "categorical"].index].fillna('')
    # Convert boolean from categorical to String
    mask = md[md_types[md_types == "categorical"].index].applymap(type) != bool
    md[md_types[md_types == "categorical"].index] = md[md_types[md_types == "categorical"].index].where(mask, md[md_types[md_types == "categorical"].index].replace({True: 'True', False: 'False'}))

    # Remove names
    md.index.names = [None]
    md_types.name = None

    # sort and filter by given samples
    md = md.reindex(samples)

    # Check if matched metadata and samples
    null_rows = md.isna().all(axis=1)
    if any(null_rows):
        # Do not remove, just inform user
        #md = md.loc[~null_rows, :]
        print_log(str(sum(null_rows)) + " samples without valid metadata")

    if md.empty or sum(null_rows) == md.shape[0]:
        print_log("No valid metadata, skipping")
        return None

    metadata = Metadata(md, md_types)
    print_log("Samples: " + str(metadata.data.shape[0]))
    print_log("Numeric Fields: " + str(metadata.get_data("numeric").shape[1]))
    print_log("Categorical Fields: " + str(metadata.get_data("categorical").shape[1]))
    return metadata


def define_metadata_types(metadata):
    # Define all COLUMN TYPES as default
    types = pd.Series(Metadata.default_type, index=metadata.columns)
    # Set types
    if str(metadata.index[0]).startswith("#"):
        # types defined on file: get values defined on the first row
        types = metadata.iloc[0]
        # Validate declared types
        idx_valid = types.isin(Metadata.valid_types)
        if not idx_valid.all():
            print_log("Invalid metadata types replaced by: " + Metadata.default_type)
            types[~idx_valid] = default_type
    else:
        # guessed types from read_table
        types[metadata.dtypes.map(is_numeric_dtype)] = "numeric"

    return types


def parse_references(cfg, tax, taxonomy, ranks):
    references = None
    if cfg is not None and "references" in cfg:
        if taxonomy == "ncbi":
            references = {}
            for desc, sf in cfg["references"].items():
                references[desc] = Reference(file=sf)
                if tax:
                    # Update taxids / get taxid from name
                    references[desc].update_taxids(update_tax_nodes(references[desc].ids, tax))
                    for i in list(references[desc].ids.keys()):
                        # lineage of all parent nodes (without itself)
                        for l in tax.lineage(i)[:-1]:
                            references[desc].add_parent(l, i)
        else:
            print_log("References only possible with ncbi taxonomy, skipping")
    else:
        print_log("No references defined in the configuration file, skipping")
    return references


def parse_controls(cfg, table):
    controls = None
    control_samples = None
    if cfg is not None and "controls" in cfg:
        controls = {}
        control_samples = {}
        for desc, cf in cfg["controls"].items():
            with open(cf, "r") as file:
                samples = file.read().splitlines()
                obs = set()
                valid_samples = set()
                for rank in table.ranks():
                    # Retrieve sub-table for every rank
                    control_table = table.get_subtable(rank, samples=samples)
                    obs.update(control_table.columns.to_list())
                    valid_samples.update(control_table.index.to_list())

                # Add control observations as a reference
                controls[desc] = Reference(ids=obs)
                control_samples[desc] = list(valid_samples)
    else:
        print_log("No controls defined in the configuration file, skipping")

    return controls, control_samples


def parse_mgnify(run_mgnify, cfg, tax, ranks):
    mgnify = None
    if run_mgnify:
        if cfg is not None and "mgnify" in cfg["external"]:
            try:
                mgnify = pd.read_table(cfg["external"]["mgnify"], header=None, names=["rank", "taxa", "biome", "count"])
            except Exception as e:
                print_log("Failed parsing MGnify database file [" + cfg["external"]["mgnify"] + "], skipping")
                print_log(str(e))
                return None
            # Filter to keep only used ranks, if provided
            if ranks:
                mgnify = mgnify.loc[mgnify['rank'].isin(ranks)]
                mgnify.reset_index(drop=True, inplace=True)
            # Convert taxids if tax is provided
            if tax:
                updated_nodes = update_tax_nodes([tuple(x) for x in mgnify[["rank", "taxa"]].to_numpy()], tax)
                mgnify["taxa"] = mgnify[["rank", "taxa"]].apply(lambda rt: updated_nodes[(rt[0], rt[1])] if updated_nodes[(rt[0], rt[1])] is not None else rt[1], axis=1)
            if mgnify.empty:
                mgnify = None
                print_log("No matches with MGnify database, skipping")
        else:
            print_log("Not defined in the configuration file, skipping")
    else:
        print_log("Not activated, skipping")
    return mgnify


def run_correlation(table, top_obs_corr):
    corr = {}
    for rank in table.ranks():
        corr[rank] = {}
        if top_obs_corr:
            top_taxids = sorted(table.get_top(rank, top_obs_corr))
            matrix = table.get_subtable(taxids=top_taxids, rank=rank)
        else:
            top_taxids = sorted(table.observations(rank))
            matrix = table.data[rank]

        corr[rank]["observations"] = top_taxids
        corr[rank]["rho"] = []
        # No correlation with just one observation
        if len(matrix.columns) >= 2:
            rho = pairwise_rho(transform_table(matrix, 0, "clr", table.zerorep).values)
            if len(matrix.columns) == 2:
                # If there are only 2 observations, return in a float
                # re-format in a matrix shape
                rho = np.array([[np.nan, np.nan], [rho[1, 0], np.nan]])
            else:
                # fill upper triangular matrix (mirrored values) with nan to be ignored by pandas
                # to save half of the space
                rho[np.triu_indices(rho.shape[0])] = np.nan

            corr[rank]["rho"] = rho

    return corr


def parse_input_file(input_file, unassigned_header, transpose, sample_replace, cumm_levels):

    if input_file.endswith(".biom"):
        table_df = biom.load_table(input_file).to_dataframe(dense=True)
    else:
        # Default input_file: index=observations, columns=samples
        # table_df should have samples on indices and observations on columns
        table_df = pd.read_table(input_file, sep='\t', index_col=0, dtype={0: str}).transpose().fillna(0)
        # Enforce string observations
        table_df.columns = table_df.columns.astype(str)

    # If user is providing a reverse table, turn back
    if transpose:
        table_df = table_df.transpose()

    # Remove header on rows
    table_df.index.name = None

    # Replace text on sample labels
    if sample_replace:
        print_log("Replacing sample values:")
        before_replace = table_df.head(1).index
        #get index as series to use replace method
        new_index = table_df.reset_index()["index"].replace(regex=dict(zip(sample_replace[::2], sample_replace[1::2])))
        table_df.set_index(new_index, inplace=True)
        for b, a in zip(before_replace, table_df.head(1).index):
            print_log("  " + b + " -> " + a)
        print_log("  ...")

    # Sum total before split unassigned or filter
    if cumm_levels:
        total = table_df.max(axis=1)
    else:
        total = table_df.sum(axis=1)

    # unique unassigned/unclassified for table
    # Separate unassigned counts column from main data frame
    unassigned = pd.Series(0, index=table_df.index)
    if unassigned_header:
        for header in unassigned_header:
            if header in table_df.columns:
                if isinstance(table_df[header], pd.DataFrame):
                    # Sum in case there are several equally named headers
                    unassigned += table_df[header].sum(axis=1)
                else:
                    # return a pd.Series
                    unassigned += table_df[header]
                table_df.drop(columns=header, inplace=True)
            else:
                print_log("'" + header + "' header not found")

    if unassigned.sum() == 0:
        print_log("No unassigned entries defined")

    print_log("Trimming table")
    table_df = trim_table(table_df)

    # Filter based on the final table
    unassigned = unassigned.reindex(table_df.index)
    total = total.reindex(table_df.index)

    return table_df, total, unassigned


def filter_input_table(table_df, total, min_frequency, max_frequency, min_count, max_count, normalized):

    if min_count:
        cnt = table_df.sum().sum()
        if min_count < 1:
            table_df_norm = transform_table(table_df, total, "norm", 0) if not normalized else table_df
            table_df = table_df[table_df_norm >= min_count].fillna(0)
        elif min_count > 1:
            table_df = table_df[table_df >= min_count].fillna(0)
        print_log(str(int(cnt - table_df.sum().sum())) + " counts skipped with --min-count " + str(min_count))

    if max_count:
        cnt = table_df.sum().sum()
        if max_count < 1:
            table_df_norm = transform_table(table_df, total, "norm", 0) if not normalized else table_df
            table_df = table_df[table_df_norm <= max_count].fillna(0)
        elif max_count > 1:
            table_df = table_df[table_df <= max_count].fillna(0)
        print_log(str(int(cnt - table_df.sum().sum())) + " counts skipped with --max-count " + str(max_count))

    if min_frequency:
        cnt = table_df.shape[1]
        table_df_freq = table_df.gt(0).sum(axis=0)
        if min_frequency < 1:
            table_df_freq = table_df_freq / table_df.shape[0]
            table_df = table_df.loc[:, table_df_freq >= min_frequency]
        elif min_frequency > 1:
            table_df = table_df.loc[:, table_df_freq >= min_frequency]
        print_log(str(int(cnt - table_df.shape[1])) + " observations removed with --min-frequency " + str(min_frequency))

    if max_frequency:
        cnt = table_df.shape[1]
        table_df_freq = table_df.gt(0).sum(axis=0)
        if max_frequency < 1:
            table_df_freq = table_df_freq / table_df.shape[0]
            table_df = table_df.loc[:, table_df_freq <= max_frequency]
        elif max_frequency > 1:
            table_df = table_df.loc[:, table_df_freq <= max_frequency]
        print_log(str(int(cnt - table_df.shape[1])) + " observations removed with --max-frequency " + str(max_frequency))

    return table_df


def trim_table(table_df):
    # Check for cols/rows with sum zero
    zero_rows = table_df.sum(axis=1).eq(0)
    if any(zero_rows):
        table_df = table_df.loc[~zero_rows, :]
        print_log(str(sum(zero_rows)) + " samples with only zero removed")

    zero_cols = table_df.sum(axis=0).eq(0)
    if any(zero_cols):
        table_df = table_df.loc[:, ~zero_cols]
        print_log(str(sum(zero_cols)) + " observations with only zero removed")

    return table_df


def parse_multi_table(table_df, ranks, tax, level_separator, obs_replace, cumm_levels):
    from grimer.grimer import _debug

    # Transpose table (obseravations as index) and expand ranks in columns
    ranks_df = table_df.T.index.str.split(level_separator, expand=True).to_frame(index=False)

    # For every pair of replace arguments
    if obs_replace:
        print_log("Replacing observation values:")
        before_replace = ranks_df.dropna().head(1).values[0]
        ranks_df.replace(regex=dict(zip(obs_replace[::2], obs_replace[1::2])), inplace=True)
        for b, a in zip(before_replace, ranks_df.dropna().head(1).values[0]):
            print_log("  " + b + " -> " + a)
        print_log("  ...")

    # replace entirely space or empty with NaN
    ranks_df = ranks_df.replace(r'^\s*$', np.nan, regex=True)

    # Set rank names, matching user defined or default
    user_ranks = False
    if len(ranks) == ranks_df.shape[1]:
        parsed_ranks = {r: ranks[r] for r in range(ranks_df.shape[1])}
        user_ranks = True
    else:
        print_log("Ranks provided (" + str(len(ranks)) + ") do not match file (" + str(ranks_df.shape[1]) + " levels). Using default named ranks.")
        parsed_ranks = {r: "rank-" + str(r) for r in range(ranks_df.shape[1])}
    ranks_df.rename(columns=parsed_ranks, inplace=True)

    # Update taxids
    if tax is not None:
        unmatched_nodes = 0
        for i, r in parsed_ranks.items():
            rank_nodes = ranks_df[r].dropna().unique()

            # If there is at least one valid entry
            if rank_nodes.any():
                # If user-provided ranks are matching, update nodes with rank
                if user_ranks:
                    updated_nodes = {node: unode for (rank, node), unode in update_tax_nodes([(r, n) for n in rank_nodes], tax).items()}
                else:
                    updated_nodes = update_tax_nodes(rank_nodes, tax)

                # Add nan to keep missing ranks (different than tax.undefined_node [None] which will keep the name)
                updated_nodes[np.nan] = np.nan
                ranks_df[r] = ranks_df[r].map(lambda t: updated_nodes[t] if updated_nodes[t] is not None else t)
                del updated_nodes[np.nan]

                unmatched_nodes += list(updated_nodes.values()).count(tax.undefined_node)

        if unmatched_nodes:
            print_log(str(unmatched_nodes) + " observations not found in taxonomy (but kept)")

    # Check unique lineage
    for i, r in parsed_ranks.items():
        if i > 0:
            lin_count = ranks_df.iloc[:, :i+1].drop_duplicates().groupby(r).count()
            invalid = lin_count[(lin_count > 1).any(axis=1)].index.to_list()
            if invalid:
                print_log(str(len(invalid)) + " observations removed with invalid lineage at " + r)
                if _debug:
                    print_log(",".join(invalid) + " observations removed with invalid lineage at " + r)
                # Set to NaN to keep shape of ranks_df
                ranks_df.loc[ranks_df[r].isin(invalid), r] = np.nan

    ranked_tables = {}
    for i, r in parsed_ranks.items():
        # ranks_df and table_df.T have the same shape
        ranked_table_df = pd.concat([ranks_df[r], table_df.T.reset_index(drop=True)], axis=1)
        if cumm_levels:
            ranked_tables[r] = ranked_table_df.groupby([r], dropna=True).max().T
        else:
            ranked_tables[r] = ranked_table_df.groupby([r], dropna=True).sum().T
        ranked_tables[r].columns.name = None

    lineage = ranks_df
    return ranked_tables, lineage


def parse_single_table(table_df, ranks, tax, default_rank_name):

    # Update taxids
    if tax is not None:
        print(table_df)
        updated_nodes = update_tax_nodes(table_df.columns, tax)
        unmatched_nodes = list(updated_nodes.values()).count(tax.undefined_node)
        if unmatched_nodes:
            print_log(str(unmatched_nodes) + " observations not found in taxonomy")

        for node, upd_node in updated_nodes.items():
            if upd_node is not None and upd_node != node:
                # If updated node is a merge on an existing taxid, sum values
                if upd_node in table_df:
                    table_df[upd_node] += table_df[node]
                    table_df.drop(columns=node, inplace=True)
                    print_log("Updated and merged taxonomic nodes: " + node + " -> " + upd_node)
                else:
                    table_df.rename(columns={node: upd_node}, inplace=True)
                    print_log("Updated taxonomic node: " + node + " -> " + upd_node)


    # Generate ranks
    ranked_tables = {}
    for rank in ranks:
        # Special case for "default" rank
        if rank == default_rank_name:
            ranked_tables[rank] = table_df
        else:
            taxid_parent_rank = {i: tax.parent_rank(tax.latest(i), rank) for i in table_df.columns}
            rank_df = pd.DataFrame(index=table_df.index)
            for taxid, parent_rank_taxid in taxid_parent_rank.items():
                if parent_rank_taxid is None:
                    #no_rank += 1
                    continue
                if parent_rank_taxid not in rank_df:
                    rank_df[parent_rank_taxid] = 0
                rank_df[parent_rank_taxid] += table_df[taxid]

            if not rank_df.empty:
                ranked_tables[rank] = rank_df


    # Generate lineage
    if tax:
        lineage = pd.DataFrame(list(map(lambda t: tax.lineage(t, ranks=list(ranked_tables.keys())), table_df.columns)), columns=list(ranked_tables.keys()))
    else:
        lineage = pd.DataFrame()

    return ranked_tables, lineage


def transform_table(df, total_counts, transformation, replace_zero_value):
    # Special case clr with one observation (result in zeros)
    if transformation == "clr" and df.shape[1] == 1:
        print_log("WARNING: using log instead of clr with one observation")
        transformation = "log"

    if transformation == "log":
        transformed_df = (df + replace_zero_value).apply(np.log10)
    elif transformation == "clr":
        transformed_df = pd.DataFrame(clr(df + replace_zero_value), index=df.index, columns=df.columns)
    elif transformation == "norm":
        transformed_df = df.divide(total_counts, axis=0) + replace_zero_value
    else:
        transformed_df = df + replace_zero_value

    return transformed_df


def update_tax_nodes(nodes, tax):
    """
    nodes can be a list of strings: taxids or names or a list of tuples with (rank, taxid/name)
    Return a dictionary mapping nodes and updated nodes (or None)
    First look for id, if nothing found, lookup by unique name
    """

    updated_nodes = {}
    for node in nodes:
        if isinstance(node, tuple):
            r = node[0]
            n = node[1]
        else:
            r = None
            n = node

        # Either returns same node, updated or tax.undefined_node (None)
        updated_taxid = tax.latest(n)
        if updated_taxid:
            # Assign updated or same taxid
            updated_nodes[node] = updated_taxid
        else:
            names = tax.search_name(n, rank=r, exact=True)
            # Assign taxid if found unique name only
            if names and len(names) == 1:
                updated_nodes[node] = names[0]
            else:
                updated_nodes[node] = tax.undefined_node

    return updated_nodes


def run_decontam(run_decontam, cfg, table, metadata, control_samples, script_dir):

    if not run_decontam:
        print_log("Not activated, skipping")
        return None

    if cfg is None:
        print_log("Not defined in the configuration file, skipping")
        return None

    df_decontam = pd.DataFrame(index=table.samples, columns=["concentration", "controls"])
    cfg_decontam = cfg["external"]["decontam"]
    tmp_output_prefix = "tmp_"

    # Collect metadata for DECONTAM (concentrations to use frequency and control for prevalence)
    out_table = tmp_output_prefix + "table_counts.tsv"
    out_concentration = tmp_output_prefix + "concentration_counts.tsv"
    out_controls = tmp_output_prefix + "control_samples_list.txt"
    if cfg_decontam["method"] in ["frequency", "combined"]:
        out_concentration = tmp_output_prefix + "concentration_counts.tsv"
        # Load frequency file, if provided
        if "frequency_file" in cfg_decontam:
            if os.path.isfile(cfg_decontam["frequency_file"]):
                # Load concentrations from file and sort (reindex) based on table inputs
                df_decontam["concentration"] = pd.read_table(cfg_decontam["frequency_file"], sep='\t', header=None, skiprows=0, index_col=0).reindex(table.samples)
                # If any entry is unknown, input is incomplete
                if df_decontam["concentration"].isnull().values.any():
                    print_log("File " + cfg_decontam["frequency_file"] + " is incomplete (Missing: " + ",".join(df_decontam[df_decontam.isnull().any(axis=1)].index.to_list()) + "), skipping")
                    return None
            else:
                print_log("File " + cfg_decontam["frequency_file"] + " not found, skipping")
                return None
        elif "frequency_metadata" in cfg_decontam:
            if cfg_decontam["frequency_metadata"] in metadata.get_col_headers():
                # Get concentrations from metadata
                df_decontam["concentration"] = metadata.get_col(cfg_decontam["frequency_metadata"])
            else:
                print_log("Could not find " + cfg_decontam["frequency_metadata"] + " in the metadata, skipping.")
                return None
        elif not table.normalized:
            # Use total from table
            print_log("No concentration provided, using total counts as concentration (frequency for DECONTAM)")
            df_decontam["concentration"] = table.get_total()
        else:
            print_log("Cannot run DECONTAM without defined concentration and normalized input values, skipping")
            return None
        # Print concentrations to file
        df_decontam["concentration"].to_csv(out_concentration, sep="\t", header=False, index=True)

    if cfg_decontam["method"] in ["prevalence", "combined"]:
        control_list = set()
        if "prevalence_file" in cfg_decontam:
            for file in cfg_decontam["prevalence_file"]:
                if os.path.isfile(file):
                    # Load controls from file
                    control_list.update([line.rstrip() for line in open(file)])
                else:
                    print_log("File not found " + file)
        elif "prevalence_metadata" in cfg_decontam:
            for field, value in cfg_decontam["prevalence_metadata"].items():
                if field in metadata.get_col_headers():
                    control_list.update(metadata.get_subset(field, value).index)
                else:
                    print_log("Could not find " + field + " in the metadata.")
        else:
            # Use all samples passed as controls
            for cs in control_samples.values():
                control_list.update(cs)

        # Select valid controls
        df_decontam["controls"] = table.samples.isin(control_list)

        if df_decontam["controls"].any():
            print_log(str(df_decontam["controls"].sum()) + " valid control samples to be used by DECONTAM")
            outf = open(out_controls, "w")
            print("\n".join(df_decontam.index[df_decontam["controls"]]), file=outf)
            outf.close()
        else:
            print_log("Could not find valid control entries, skipping")
            return None

    decontam = Decontam(df_decontam)
    # Run DECONTAM for each for each
    for rank in table.ranks():
        if len(table.observations(rank)) == 1:
            decontam.add_rank_empty(rank, table.observations(rank))
        else:
            # normalize and write temporary table for each rank
            if not table.normalized:
                transform_table(table.data[rank], table.get_total()[table.data[rank].index], "norm", 0).to_csv(out_table, sep="\t", header=True, index=True)
            else:
                table.data[rank].to_csv(out_table, sep="\t", header=True, index=True)

            cmd = " ".join([os.path.join(script_dir, "scripts", "run_decontam.R"),
                            "--resout " + tmp_output_prefix + "decontam_out.tsv",
                            "--modout " + tmp_output_prefix + "decontam_mod.tsv",
                            "--counts " + out_table,
                            "--concentrations " + out_concentration if cfg_decontam["method"] in ["frequency", "combined"] else "",
                            "--controls " + out_controls if cfg_decontam["method"] in ["prevalence", "combined"] else "",
                            "--method " + cfg_decontam["method"],
                            "--threshold " + str(cfg_decontam["threshold"])])
            stdout, stderr = run_cmd(cmd)

            decontam.add_rank_results(rank, tmp_output_prefix + "decontam_out.tsv", tmp_output_prefix + "decontam_mod.tsv")

    for file in [out_table, out_concentration, out_controls, tmp_output_prefix + "decontam_out.tsv", tmp_output_prefix + "decontam_mod.tsv"]:
        if os.path.isfile(file):
            os.remove(file)

    return decontam


def run_hclustering(table, linkage_methods, linkage_metrics, transformation, skip_dendrogram, optimal_ordering):
    hcluster = {}
    dendro = {}

    for rank in table.ranks():

        # Get .values of transform, numpy array
        matrix = transform_table(table.data[rank], table.get_total(), transformation, table.zerorep).values

        hcluster[rank] = {}
        dendro[rank] = {}
        for method in linkage_methods:
            hcluster[rank][method] = {}
            dendro[rank][method] = {}
            for metric in linkage_metrics:
                hcluster[rank][method][metric] = {}
                hcluster[rank][method][metric]["x"] = {}
                hcluster[rank][method][metric]["y"] = {}

                #H.clustering, returning dendrogram
                # Only one observation does not cluster
                if matrix.shape[1] > 1:
                    x = sch.dendrogram(sch.linkage(matrix.transpose(), method=method, metric=metric, optimal_ordering=optimal_ordering), no_plot=True)
                    hcluster[rank][method][metric]["x"]["index"] = table.observations(rank)[x["leaves"]].to_list()
                else:
                    hcluster[rank][method][metric]["x"]["index"] = table.observations(rank).to_list()

                # Only one samples does not cluster
                if matrix.shape[0] > 1:
                    y = sch.dendrogram(sch.linkage(matrix, method=method, metric=metric, optimal_ordering=optimal_ordering), no_plot=True)
                    hcluster[rank][method][metric]["y"]["index"] = table.samples[y["leaves"]].to_list()
                else:
                    hcluster[rank][method][metric]["y"]["index"] = table.samples.to_list()

                if not skip_dendrogram:
                    dendro[rank][method][metric] = {}
                    dendro[rank][method][metric]["y"] = {}
                    dendro[rank][method][metric]["x"] = {}

                    # Save dendrogram values and colors
                    xs, ys, colors = [[]] * 3
                    if matrix.shape[1] > 1:
                        xs, ys, colors = dendro_lines_color(x, "x")
                    dendro[rank][method][metric]["x"]["xs"] = xs
                    dendro[rank][method][metric]["x"]["ys"] = ys
                    dendro[rank][method][metric]["x"]["colors"] = colors
                    if matrix.shape[0] > 1:
                        xs, ys, colors = dendro_lines_color(y, "y")
                    dendro[rank][method][metric]["y"]["xs"] = xs
                    dendro[rank][method][metric]["y"]["ys"] = ys
                    dendro[rank][method][metric]["y"]["colors"] = colors

    return hcluster, dendro


def dendro_lines_color(dendro, axis):
    icoord = pd.DataFrame(dendro["icoord"])
    icoord = icoord * ((len(dendro["icoord"]) + 0.5) / icoord.max().max())
    icoord = icoord.values.tolist()
    if axis == "y":
        dcoord = dendro["dcoord"]
    else:
        dcoord = [[-j for j in i] for i in dendro['dcoord']]

    color_list = dendro["color_list"]
    unique_colors = sorted(set(color_list))
    cp = make_color_palette(len(unique_colors))
    colors = [cp[unique_colors.index(colorid)] for colorid in color_list]

    if axis == "y":
        return dcoord, icoord, colors
    else:
        return icoord, dcoord, colors


def pairwise_vlr(mat):
    cov = np.cov(mat.T, ddof=1)
    diagonal = np.diagonal(cov)
    return -2 * cov + diagonal[:, np.newaxis] + diagonal


def pairwise_rho(mat):
    variances = np.var(mat, axis=0, ddof=1)
    return 1 - (pairwise_vlr(mat) / np.add.outer(variances, variances))


def include_scripts(scripts):
    # Insert global js functions and css and return template
    template = "{% block postamble %}"
    for file, t in scripts.items():
        with open(file, 'r') as file:
            template += "<" + t + ">"
            template += "".join(file.readlines())
            template += "</" + t + ">"
    template += "{% endblock %}"
    return template


def format_js_toString(val):
    # Transform numeric value to float and string to match toString
    return str(float(val)) if isinstance(val, (int, float)) else str(val)


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

def run_cmd(cmd, print_stderr: bool=False, exit_on_error: bool=True):
    errcode = 0
    stdout = ""
    stderr = ""
    try:
        process = subprocess.Popen(shlex.split(cmd),
                                   universal_newlines=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        # wait for the process to terminate
        stdout, stderr = process.communicate()
        errcode = process.returncode
        if exit_on_error and errcode != 0:
            raise Exception()
        if print_stderr and stderr:
            print_log(stderr)

    except Exception as e:
        print_log('The following command failed to run:\n' + cmd)
        print_log(str(e))
        print_log("Error code: " + str(errcode))
        print_log("Out: ")
        if stdout:
            print_log(stdout)
        print_log("Error: ")
        if stderr:
            print_log(stderr)
        sys.exit(errcode)

    return stdout, stderr


def print_log(text):
    sys.stderr.write(text + "\n")
    sys.stderr.flush()


def print_df(df, name: str=None):
    from grimer.grimer import _debug
    if _debug:
        print(name)
        if isinstance(df, dict):
            if df:
                print(len(df.keys()), "keys:", list(df.keys())[0], "...", list(df.keys())[-1])
                #print(list(df.values())[0], "...", list(df.values())[-1])
        else:
            #print(df.columns)
            print(df.head())
            print(df.shape)
        print("size:", sys.getsizeof(df))
        print("-----------------------------------------------")


def print_logo_cli(version):
    print_log("==================")
    print_log(" ╔═╗╦═╗╦╔╦╗╔═╗╦═╗ ")
    print_log(" ║ ╦╠╦╝║║║║║╣ ╠╦╝ ")
    print_log(" ╚═╝╩╚═╩╩ ╩╚═╝╩╚═ ")
    print_log("           v" + version)
    print_log("==================")
