import pandas as pd


class Decontam:
    cols_rank = ["freq", "prev", "p.freq", "p.prev", "p", "contaminant"]

    def __init__(self, df_concentration_controls):
        self.data = df_concentration_controls
        self.rank = {}

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k, v) in vars(self).items()]
        return 'Decontam({})'.format(', '.join(args))

    def add_rank_results(self, rank, decontam_out_file, decontam_mod_file):
        self.rank[rank] = pd.read_table(decontam_out_file, sep='\t', header=0, skiprows=0, index_col=0, names=self.cols_rank, dtype={0: str})

        # Parse models enforcing index as string
        mod = pd.read_table(decontam_mod_file, sep='\t', header=0, skiprows=0, index_col=0, dtype={0: str})

        # Remove point counter at end (.1 or .1000)
        mod.index = mod.index.map(lambda txid: txid[:-5] if txid.endswith(".1000") else txid[:-2]).to_list()

        # Merge first point of model
        self.rank[rank] = self.rank[rank].merge(mod.iloc[0::2, 0], left_index=True, right_index=True)

        # Merge second point of model and non-contant line
        self.rank[rank] = self.rank[rank].merge(mod.iloc[1::2, :], suffixes=["", "_2"], left_index=True, right_index=True)

    def add_rank_empty(self, rank, idx):
        self.rank[rank] = pd.DataFrame(index=idx, columns=self.cols_rank + ["contam", "contam_2", "non.contam"])
        self.rank[rank]["contaminant"] = False

    def get_data(self):
        return self.data.fillna(False)

    def get_contaminants(self, rank, idx):
        return self.rank[rank].reindex(idx)["contaminant"]

    def get_pvalues(self, rank, idx):
        return self.rank[rank].reindex(idx)["p"]

    def get_contaminant_list(self):
        clist = []
        for r in self.rank:
            clist.extend(self.rank[r].index[self.rank[r]["contaminant"]==True].to_list())
        return clist
