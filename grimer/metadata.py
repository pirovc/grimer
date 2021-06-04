import pandas as pd
from pandas.api.types import is_numeric_dtype
from grimer.utils import print_log


class Metadata:
    valid_types = ["categorical", "numeric"]
    default_type = "categorical"

    def __init__(self, metadata_file, samples: list=[]):
        # Read metadata and let pandas guess dtypes, index as str
        self.data = pd.read_table(metadata_file, sep='\t', header=0, skiprows=0, index_col=0, dtype={0:str})

        # Enforce string index
        self.data.index = self.data.index.astype('str')

        # Define all COLUMN TYPES as default
        self.types = pd.Series(self.default_type, index=self.data.columns)

        # Set types
        if str(self.data.index[0]).startswith("#"):
            # types defined on file
            self.set_hard_types()
        else:
            # guessed types from read_table
            self.types[self.data.dtypes.map(is_numeric_dtype)] = "numeric"

        # Convert datatypes to adequate numeric values (int, float)
        self.data = self.data.convert_dtypes(infer_objects=False, convert_string=False)
        # Re-convert everython to object to standardize (int64 NA is not seriazable on bokeh)
        self.data = self.data.astype("object")

        # Remove empty fields
        null_cols = self.data.isna().all(axis=0)
        if any(null_cols):
            self.data = self.data.loc[:, ~null_cols]
            self.types = self.types[~null_cols]
            print_log(str(sum(null_cols)) + " fields removed without valid values")

        # Convert NaN on categorical to ""
        self.data[self.types[self.types == "categorical"].index] = self.data[self.types[self.types == "categorical"].index].fillna('')

        # Remove names
        self.data.index.names = [None]
        self.types.name = None

        # sort and filter by given samples
        if samples:
            self.data = self.data.reindex(samples)

        # Check if matched metadata and samples
        null_rows = self.data.isna().all(axis=1)
        if any(null_rows):
            #self.data = self.data.loc[~null_rows, :]
            print_log(str(sum(null_rows)) + " samples without valid metadata")

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k, v) in vars(self).items()]
        return 'Metadata({})'.format(', '.join(args))

    def set_hard_types(self):
        # Get values defined on the first row
        self.types = self.data.iloc[0]
        # Drop row with types from main data
        self.data.drop(self.types.name, inplace=True)
        # Validate declared types
        idx_valid = self.types.isin(self.valid_types)
        if not idx_valid.all():
            print_log("Invalid metadata types replaced by: " + self.default_type)
            self.types[~idx_valid] = self.default_type
        # Enforce column type on dataframe
        self.data[self.types[self.types == "categorical"].index] = self.data[self.types[self.types == "categorical"].index].astype(str)
        self.data[self.types[self.types == "numeric"].index] = self.data[self.types[self.types == "numeric"].index].apply(pd.to_numeric)

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
        return sorted(self.get_col(col).dropna().unique())

    def get_formatted_unique_values(self, col):
        if self.types[col] == "categorical":
            return self.get_unique_values(col)
        else:
            return list(map('{:.16g}'.format, self.get_unique_values(col)))

    def get_type(self, col):
        return self.types[col]

    def get_subset(self, column, value):
        return self.data[self.data[column] == value]
