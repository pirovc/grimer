import yaml


class SourceOld:
    def __init__(self, file: str=None, ids: list=[]):
        # Only leaf ids/nodes
        self.ids = set()
        self.lineage = set()

        # {id: {ref1: set(desc1, desc2,...), ref2: set(desc3,...)}
        self.refs = {}

        if file is not None:
            self.parse(file)
        elif ids:
            self.ids.update(ids)
            self.lineage.update(ids)

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k, v) in vars(self).items()]
        return 'Source({})'.format(', '.join(args))

    def parse(self, file):
        with open(file, 'r') as fh:
            if file.endswith(".yml") or file.endswith(".yaml"):
                src = yaml.safe_load(fh)
                for desc, val in src.items():
                    for ref, v in val.items():
                        str_ids = list(map(str, v["ids"]))
                        self.ids.update(str_ids)
                        self.lineage.update(str_ids)
                        for i in str_ids:
                            self.add_refs_desc(i, (ref, v["url"]), desc)
            else:
                for line in fh:
                    main_id = line.rstrip()
                    self.ids.add(main_id)
                    self.lineage.add(main_id)

    def update_lineage(self, ids):
        self.lineage.update(ids)

    def update_taxids(self, taxid_updated):
        # Update taxonomy entries or convert names to taxid
        for node, upd_node in taxid_updated.items():
            if upd_node is not None and upd_node != node:
                print("Updated taxonomic node: " + node + " -> " + upd_node)
                self.ids.discard(node)
                self.ids.add(upd_node)
                self.lineage.discard(node)
                self.lineage.add(upd_node)
                if node in self.refs:
                    self.refs[upd_node] = self.refs.pop(node)

    def add_refs_desc(self, i, ref, desc):
        if i not in self.refs:
            self.refs[i] = {}
        if ref not in self.refs[i]:
            self.refs[i][ref] = set()
        if desc is not None:
            self.refs[i][ref].add(desc)

    def get_refs_desc(self, i):
        return self.refs[i] if i in self.refs else {}

    def get_refs(self, i):
        return list(self.refs[i].keys()) if i in self.refs else ()

    def get_refs_count(self, i):
        return len(self.refs[i]) if i in self.refs else 0

    def update_refs(self, taxid_parent_rank):
        for taxid, parent_taxid in taxid_parent_rank.items():
            if parent_taxid is not None and taxid in self.refs:
                for i in self.refs[taxid]:
                    for r in self.refs[taxid][i]:
                        self.add_refs_desc(parent_taxid, i, r)
