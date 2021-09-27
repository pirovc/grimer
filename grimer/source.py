import yaml


class Source:
    def __init__(self, file: str=None, ids: list=[]):
        self.ids = {}  # {refid: {ref1: set(desc1, desc2,...), ref2: set(desc3,...)}}
        self.children = {}  # {child_id: set(refids)}
        self.parents = {}  # {parent_id: set(refids)}

        if file is not None:
            self.parse(file)
        elif ids:
            for i in ids:
                self.add(i)

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k, v) in vars(self).items()]
        return 'Source({})'.format(', '.join(args))

    def add(self, i, ref: str=None, desc: str=None):
        if i not in self.ids:
            self.ids[i] = {}
        if ref is not None:
            if ref not in self.ids[i]:
                self.ids[i][ref] = set()
            if desc is not None:
                self.ids[i][ref].add(desc)

    def add_child(self, child, refid):
        if child not in self.children:
            self.children[child] = set()
        self.children[child].add(refid)

    def add_parent(self, parent, refid):
        if parent not in self.parents:
            self.parents[parent] = set()
        self.parents[parent].add(refid)

    def parse(self, file):
        with open(file, 'r') as fh:
            if file.endswith(".yml") or file.endswith(".yaml"):
                src = yaml.safe_load(fh)
                for desc, val in src.items():
                    for ref, v in val.items():
                        for i in map(str, v["ids"]):
                            self.add(i, (ref, v["url"]), desc)
            else:
                for line in fh:
                    main_id = line.rstrip()
                    self.add(main_id)

    def update_taxids(self, taxid_updated):
        for node, upd_node in taxid_updated.items():
            if upd_node is not None and upd_node != node:
                print("Updated taxonomic node: " + node + " -> " + upd_node)
                self.add(upd_node)
                self.ids[upd_node].update(self.ids[node])
                self.ids.discard(node)

    def get_refs_desc(self, i, direct: bool=False, children: bool=False, parents: bool=False):
        refs_desc = {}
        if direct and i in self.ids:
            refs_desc.update(self.ids[i])
        if children and i in self.children:
            for refid in self.children[i]:
                refs_desc.update(self.ids[refid])
        if parents and i in self.parents:
            for refid in self.parents[i]:
                refs_desc.update(self.ids[refid])
        return refs_desc

    def get_refs_count(self, i, direct: bool=False, children: bool=False, parents: bool=False):
        return len(self.get_refs_desc(i, direct, children, parents))
