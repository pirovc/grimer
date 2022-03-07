class Metadata:
    valid_types = ["categorical", "numeric"]
    default_type = "categorical"

    def __init__(self, md, types):
        self.data = md
        self.types = types

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k, v) in vars(self).items()]
        return 'Metadata({})'.format(', '.join(args))

    def get_col_headers(self):
        return self.data.columns

    def get_data(self, metadata_type: str=None):
        if metadata_type is not None:
            return self.data[self.types[self.types == metadata_type].index]
        else:
            return self.data

    def get_col(self, col):
        return self.data[col]

    def get_unique_values(self, col):
        return self.get_col(col).dropna().unique()

    def get_type(self, col):
        return self.types[col]

    def get_subset(self, column, value):
        return self.data[self.data[column] == value]
