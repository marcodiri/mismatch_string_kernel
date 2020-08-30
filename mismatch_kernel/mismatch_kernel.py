from math import sqrt
from mismatch_kernel.mismatch_tree import MismatchTree


class MismatchKernel:
    def __init__(self, a, k, m, mismatch_vectors=None, kernel_matrix=None):
        """
        (k-m)-mismatch kernel

        :param a: alphabet which will compose the k-mers
        :param k: length of the k-mers (substrings of length k)
        :param m: maximum number of mismatches between k-mers
        :param mismatch_vectors: (optional) a dictionary with some strings as keys
            and the corresponding vector in DOK (Dictionary of keys) format as values.
            When calling get_kernel(x1, x2) if x1 and/or x2 are present as a key in this
            dictionary, the corresponding vector will be used to calculate the kernel;
            if not present the vector will be calculated.
        :param kernel_matrix: (optional) a dictionary with some strings as keys
            and as values a dictionary with the same strings as keys and the kernel
            as values.
            Example: suppose we have strings {"a", "b", "c"}, a kernel matrix can have
            all or only some of the kernel values between these strings;
            i.e. a full kernel_matrix could be
            {"a": {"a": 1, "b": 2, "c": -1}, "b": {"b": 1, "c": 3}, "c": {"c": 1}}
            meaning that kernel("a", "a")=1, kernel("a", "b")=2, ecc..
            Note that "b" dictionary has no "a" key because kernel operation is commutative
            so kernel("b", "a")=kernel("a", "b") that is already present in the "a"
            dictionary (the actual matrix would be symmetric).
            When calling get_kernel(x1, x2), if kernel_matrix[x1][x2] or
            kernel_matrix[x2][x1] exists, the corresponding value will be immediately
            returned; otherwise it will be calculated.
        """
        self.mismatch_tree = MismatchTree(a, k, m)
        if mismatch_vectors is None:
            mismatch_vectors = {}
        if kernel_matrix is None:
            kernel_matrix = {}
        self.MISMATCH_VECTORS = mismatch_vectors
        self.KERNEL_MATRIX = kernel_matrix

    def vectorize(self, x):
        # if mismatch vector not already computed for an input
        # compute it and save it
        x = self.mismatch_tree.normalize_input(x)
        if x in self.MISMATCH_VECTORS:
            return x, self.MISMATCH_VECTORS[x]
        else:
            x, v = self.mismatch_tree.vectorize(x)
            self.MISMATCH_VECTORS[x] = v
            return x, v

    def get_kernel(self, x1, x2):
        x1 = self.mismatch_tree.normalize_input(x1)
        x2 = self.mismatch_tree.normalize_input(x2)

        # look up matrix for already computed kernel
        if x1 in self.KERNEL_MATRIX and x2 in self.KERNEL_MATRIX[x1]:
            return self.KERNEL_MATRIX[x1][x2]
        elif x2 in self.KERNEL_MATRIX and x1 in self.KERNEL_MATRIX[x2]:
            return self.KERNEL_MATRIX[x2][x1]

        x1, v1 = self.vectorize(x1)
        x2, v2 = self.vectorize(x2)

        # since mismatch vectors are in DOK format we need to compute
        # the dot product accordingly
        kernel = 0
        # iterate on the shortest dok
        dok1, dok2 = (v1, v2) if len(v1) < len(v2) else (v2, v1)
        # normalizer_dok1 = 0  # to compute normalized kernel
        # normalizer_dok2 = sum(e**2 for e in dok2.values())
        for i in dok1:
            # normalizer_dok1 += dok1[i]**2
            # if both indexes contains non zero elements, sum their product
            if i in dok2:
                kernel += dok1[i]*dok2[i]

        # normalized_kernel = kernel / (sqrt(normalizer_dok1) * sqrt(normalizer_dok2)) if kernel != 0 else 0
        # d = np.sqrt(2-2*np.around(normalized_kernel, 3))  # distance between x1 and x2

        if x1 not in self.KERNEL_MATRIX:
            self.KERNEL_MATRIX[x1] = {}
        self.KERNEL_MATRIX[x1][x2] = kernel

        return kernel
