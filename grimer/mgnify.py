import pandas as pd


class MGnify:

    def __init__(self, mgnify_file, ranks: list=[]):
        self.data = self.parse(mgnify_file, ranks)

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k, v) in vars(self).items()]
        return 'MGnify({})'.format(', '.join(args))

    def parse(self, file, ranks):
        mgnify_df = pd.read_table(file, header=None, names=["rank", "taxa", "biome", "count"])
        # Filter by ranks if provided
        if ranks:
            mgnify_df = mgnify_df.loc[mgnify_df['rank'].isin(ranks)]
            mgnify_df.reset_index(drop=True, inplace=True)

        #mgnify_df.drop(columns="rank", inplace=True)
        return mgnify_df

    def update_taxids(self, taxid_updated):
        # Update taxonomy to taxid or keep name if not available
        self.data["taxa"] = self.data[["rank", "taxa"]].apply(lambda rt: taxid_updated[(rt[0], rt[1])] if taxid_updated[(rt[0], rt[1])] is not None else rt[1], axis=1)


