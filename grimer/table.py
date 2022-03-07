from collections import OrderedDict


class Table:
    def __init__(self, samples, total, unassigned, lineage, normalized, zerorep):
        # Ordered dict to keep rank insert order
        self.data = OrderedDict()
        self.lineage = lineage
        self.samples = samples
        self.total = total
        self.unassigned = unassigned
        self.normalized = normalized
        self.zerorep = zerorep

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k, v) in vars(self).items()]
        return 'Table({})'.format(', '.join(args))

    def add_rank(self, rank, table):
        self.data[rank] = table

    def observations(self, rank):
        return self.data[rank].columns

    def ranks(self):
        return list(self.data.keys())

    def get_min_valid_count_perc(self):
        return min([self.get_counts_perc(rank)[self.get_counts_perc(rank) > 0].min().min() for rank in self.ranks()])

    def get_total(self):
        return self.total

    def get_unassigned(self):
        return self.unassigned

    def get_assigned(self):
        return self.get_total() - self.get_unassigned()

    def get_unassigned_perc(self):
        return self.get_unassigned().divide(self.get_total(), axis=0) if not self.normalized else self.get_unassigned().divide(100, axis=0)

    def get_assigned_perc(self):
        return self.get_assigned().divide(self.get_total(), axis=0) if not self.normalized else self.get_assigned().divide(100, axis=0)

    def get_lineage(self, taxid, rank, other_rank):
        # get lineage up-to requested rank
        return self.lineage[self.lineage[rank] == taxid][other_rank].values[0]

    def get_frequency(self, rank):
        return self.data[rank].gt(0).sum(axis=0)

    def get_frequency_perc(self, rank):
        return self.get_frequency(rank) / len(self.samples)

    def get_counts(self, rank):
        return self.data[rank].sum(axis=0) if not self.normalized else 0

    def get_counts_perc(self, rank):
        return self.data[rank].divide(self.get_total(), axis=0) if not self.normalized else self.data[rank].divide(100, axis=0)

    def get_counts_perc_avg_samples(self, rank):
        return self.get_counts_perc(rank).sum(axis=0) / len(self.samples)

    def get_top(self, rank, n):
        return sorted(self.get_counts_perc_avg_samples(rank).sort_values(ascending=False).index[:n].to_list())

    def get_subtable(self, rank, samples: list=[], taxids: list=[], keep_shape: bool=False):
        subtable = self.data[rank]

        if samples:
            valid_samples = []
            for s in samples:
                if s in self.samples:
                    valid_samples.append(s)

            if valid_samples:
                subtable = subtable.loc[subtable.index.intersection(valid_samples)]
                if not keep_shape:
                    subtable = subtable.loc[:, subtable.sum(axis=0) > 0]
            else:
                return None

        if taxids:
            subtable = subtable[taxids].copy()
            if not keep_shape:
                subtable = subtable.loc[subtable[taxids].sum(axis=1) > 0, :]

        return subtable
