# Mismatch string kernel
This module can be used to vectorize strings and compute kernel between them.

A Python3 implementation of the mismatch kernel described in the publication below:

    %0 Journal Article
    %T Mismatch string kernels for discriminative protein classification
    %A Leslie, Christina S
    %A Eskin, Eleazar
    %A Cohen, Adiel
    %A Weston, Jason
    %A Noble, William Stafford
    %J Bioinformatics
    %V 20
    %N 4
    %P 467-476
    %@ 1460-2059
    %D 2004
    %I Oxford University Press
    %U https://doi.org/10.1093/bioinformatics/btg431

## Usage
To understand the technicalities of what this kernel does please refer to the article above.
### Initializing the kernel
First you have to define an alphabet from which the k-mers will be generated, 
the length k of the k-mers and m the maximum number of mismatches between mers, for example:
```python
ALPHABET = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
            'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', ' ']
k = 3
m = 1
```

Then you can create a `MismatchKernel` object with such parameters:
```python
from mismatch_kernel import MismatchKernel

mk = MismatchKernel(ALPHABET, k, m)
```

### Mapping a string to (k-m) feature space
You can use the `vectorize(x)` function to map a string `x` to the (k-m) feature space.

Note that the alphabet is in general case sensitive, so if your strings needs to be 
case sensitive (i.e. "string" != "StRiNg"), your alphabet should contain both uppercase
and lowercase letters. Also this will much increase computational time because the 
k-mer feature space has dimension #(ALPHABET)^k; same thing goes if you need to 
distinguish punctuation, for example in the alphabet above the strings will be different
based on the spaces they contain (i.e. "space" != "spa ce"). 
In general the strings you pass to this module functions will be normalized, i.e. every
character not in the alphabet will be removed. For example if you call `vectorize("String")`
after defining the above alphabet you are actually vectorizing "tring", so you should have
called `vectorize("String".lower())` instead.

The `vectorize(x)` function returns a tuple `(x_norm, dok)` where `x_norm` is the actual string
that has been vectorized (i.e. `x` normalized) so you can check if that's what you actually 
wanted to vectorize, and `dok` is the vector in DOK (dictionary of keys)
format (because the vectors are generally sparse), so it will be a dictionary like
`{2: 1, 3: 1, 14: 1, 17: 2, 30: 1, 41: 1, ...}` meaning that the vector has non-zero values 
only at the position of the dictionary keys, i.e. `[0, 0, 1, 1, 0, ..., 0, 1, 0, 0, 2, ...]`.
You can push `x_norm` in a dictionary along with the vector so you don't have to 
vectorize it again, this is what the `get_kernel()` function actually does.
#### Example
```python
x_norm, vect = mk.mismatch_tree.vectorize("doc. Frankenstein".lower())
print("{} -> {}".format(x_norm, vect))
```

    > doc frankenstein -> {10: 1, 13: 1, 37: 1, 64: 1, ...}

### Calculating the kernel between two strings
You can use the `get_kernel(x1, x2)` function to get the kernel between `x1` and `x2`,
the kernel varies between 0 and 1, the more similar the two strings the greater it will be 
(1 if the strings are equal).
The function will automatically normalize and vectorize the two strings to compute the 
kernel.

#### Example
```python
ker = mk.get_kernel("doc. Frankenstein".lower(), "doc. Drunkenstein".lower())
print(ker)
```
    
    > 0.7500011542039571

### Using or supplying already calculated mismatch vectors and kernels
The `get_kernel` function will save in the `MismatchKernel` object the mismatch vectors of every
string it vectorizes in the `MISMATCH_VECTORS` attribute, that is a dictionary that stores
strings as keys and the corresponding vector as values 
(i.e. `{'doc frankenstein': {10: 1, 13: 1, 37: 1, 64: 1, ...}, 
'doc drunkenstein': {80: 1, 98: 1, 116: 1, 121: 1, ...}}`) so if you call next
`mk.get_kernel("doc drunkenstein", "doc nykterstein"` it won't vectorize again `"doc drunkenstein"`.

Likewise every calculated kernel will be stored in the `KERNEL_MATRIX` attribute, that is a 
dictionary that stores strings as keys and another dictionary with strings as keys and the
kernel value between the two keys as values 
(i.e. `{'doc frankenstein': {'doc drunkenstein': 0.7500011542039571, 'doc nykterstein': 0.5041614599291009}}`).
If you have to calculate the kernel for a batch of strings you can call `get_kernel` from the
same `MismatchKernel` object so the strings for which the mismatch vector or the kernel have
already been calculated won't be calculated again.

If you already have one or both of these dictionaries, for example if you pickled the `MISMATCH_VECTORS` 
and `KERNEL_MATRIX` attributes from a previous run, you can pass them to the `MismatchKernel`
constructor:
```python
mk = MismatchKernel(ALPHABET, k, m, vectors_dict, kernels_dict)
```