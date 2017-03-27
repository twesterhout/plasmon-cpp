#!/usr/bin/env python3


import random
import subprocess

import numpy as np



def random_matrix(element_type, n, m):
    """
    Generates a random n x m matrix

    element_type - type of elements in the matrix. Should be one
    of 'float', 'double', 'complex-float', 'complex-double'. There
    is no difference in the behavior for 'float' and 'double', as
    python does not distinguish these type. It, however, matters
    when passing element_type to C++.
    n - height of the matrix (number of rows)
    m - width of the matrix (number of columns)

    Returns an np.ndarray of dimensions (n, m)
    """

    if (element_type == 'float') \
       or (element_type == 'double'):
        return np.random.rand(n, m)
    elif (element_type == 'complex-float') \
         or (element_type == 'complex-double'):
        return np.random.rand(n, m) + 1.j * np.random.rand(n, m)
    else:
        raise Exception('Invalid type: ' + element_type)



def to_number(x, element_type):
    """
    Converts a string to a number.

    x - string to convert
    element_type - type of number. It should be one of
    'float', 'double', 'complex-float', 'complex-double'.

    Returns a number of type element_type
    """
    def to_complex(x):
        x = x.strip().strip('(').strip(')')
        real, imag = x.split(',')
        return float(real) + 1j * float(imag)

    if (element_type == 'float') \
       or (element_type == 'double'):
        return float(x)
    elif (element_type == 'complex-float') \
         or (element_type == 'complex-double'):
        return to_complex(x)
    else:
        raise Exception('Invalid type: ' + element_type)


def to_matrix(input, element_type):
    """
    Converts string to matrix of element_type

    input - input string
    element_type - type of element

    Returns an np.ndarray of element_type
    """
    return np.array(
        [ [ to_number(x, element_type) 
            for x in l.split('\t') if x != '' ]
          for l in input.strip('\n').split('\n') ]
        )


def to_cxx_input(array, element_type):
    ls = [x for x in np.nditer(array, order='C')]

    to_str = None
    if (element_type == 'float') \
       or (element_type == 'double'):
        to_str = str
    elif (element_type == 'complex-float') \
         or (element_type == 'complex-double'):
        to_str = lambda x: '({}, {})'.format(x.real, x.imag)
    else:
        raise Exception('Invalid type!')

    return ''.join((to_str(e) + '\t') for e in ls) + '\n'
