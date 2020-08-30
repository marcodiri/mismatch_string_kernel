import re


class MismatchTree:
    def __init__(self, alphabet, k, m):
        self.alphabet = alphabet
        self.k = k
        self.m = m

    # the only purpose is to compute the mismatch vector,
    # so it's useless to save the whole tree.
    # No nodes hold reference to other nodes, so they get
    # garbage collected as soon as they are overridden
    class Node:
        def __init__(self):
            self.table = {}
            self.height = 0

        # this function only uses the parent table to compute
        # the child table. No reference to parent is kept.
        def add_child(self, node, edge, m):
            node.height = self.height + 1
            for mer in self.table.keys():
                count = self.table[mer]
                mer = list(mer.replace(".", ""))  # remove previously attached dots
                first_letter = mer.pop(0)
                remaining_mer = "".join(mer)  # rebuild string from chars array
                while remaining_mer in node.table.keys():  # avoid override in the dictionary
                    remaining_mer += "."
                node.table[remaining_mer] = count + 1 if first_letter != edge else count
                if node.table[remaining_mer] > m:
                    node.table.pop(remaining_mer)

    def normalize_input(self, x):
        letters = "".join(self.alphabet)
        return re.sub(r"[^{}]".format(letters), "", x)

    def vectorize(self, x):
        root = self.Node()
        # if running into memory issues for large K, remove mismatch_vector
        # since it is not necessary to keep it in memory since we use the dok
        mismatch_vector = []
        # for K > 2 mismatch_vector are sparse, so compress it
        # using Dictionary of keys (DOK) because when computing
        # the kernel we need fast access to check whether the two
        # of these vectors have a non zero element at a given index
        dok = {}

        if len(x) < self.k:
            raise RuntimeError("Cannot do {}-mer of a {} length string".format(self.k, len(x)))

        x_norm = self.normalize_input(x)

        def k_mers(string):
            a = {}
            for i in range(len(string)-self.k+1):
                a[string[i:i+self.k]] = 0  # duplicate k-mers get overridden
            return a

        root.table = k_mers(x_norm)

        def add_child(node):
            # build the tree in a depth-first fashion
            for a in self.alphabet:
                # each cycle child_node gets overridden and since
                # no references are kept it gets garbage collected.
                # The only nodes kept in memory are the nodes on
                # the path from the root to the leaf along the edges
                # that form the current k-mer (one of the
                # len(alphabet)^k possible k-mers).
                # So at any time max k+1 = O(k) nodes are kept in memory.
                child_node = self.Node()
                node.add_child(child_node, a, self.m)
                if child_node.height < self.k:
                    add_child(child_node)
                else:  # max depth reached, save number of survived mers
                    survived_mers = len(child_node.table)
                    mismatch_vector.append(survived_mers)
                    if survived_mers > 0:
                        dok[len(mismatch_vector)-1] = survived_mers

        add_child(root)

        return x_norm, dok
