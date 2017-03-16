"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

Porechop makes use of C++ functions which are compiled in cpp_functions.so. This module uses ctypes
to wrap them in similarly named Python functions.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
from ctypes import CDLL, cast, c_char_p, c_int, c_void_p

SO_FILE = 'cpp_functions.so'
SO_FILE_FULL = os.path.join(os.path.dirname(os.path.realpath(__file__)), SO_FILE)
if not os.path.isfile(SO_FILE_FULL):
    sys.exit('could not find ' + SO_FILE + ' - please reinstall')
C_LIB = CDLL(SO_FILE_FULL)

C_LIB.adapterAlignment.argtypes = [c_char_p,  # Read sequence
                                   c_char_p,  # Adapter sequence
                                   c_int,     # Match score
                                   c_int,     # Mismatch score
                                   c_int,     # Gap open score
                                   c_int]     # Gap extension score
C_LIB.adapterAlignment.restype = c_void_p     # String describing alignment


# This function cleans up the heap memory for the C strings returned by the other C functions. It
# must be called after them.
C_LIB.freeCString.argtypes = [c_void_p]
C_LIB.freeCString.restype = None


def adapter_alignment(read_sequence, adapter_sequence, scoring_scheme_vals):
    """
    Python wrapper for adapterAlignment C++ function.
    """
    match_score = scoring_scheme_vals[0]
    mismatch_score = scoring_scheme_vals[1]
    gap_open_score = scoring_scheme_vals[2]
    gap_extend_score = scoring_scheme_vals[3]
    ptr = C_LIB.adapterAlignment(read_sequence.encode('utf-8'), adapter_sequence.encode('utf-8'),
                                 match_score, mismatch_score, gap_open_score, gap_extend_score)
    result_string = c_string_to_python_string(ptr)
    return result_string


def c_string_to_python_string(c_string):
    """
    This function casts a C string to a Python string and then calls a function to delete the C
    string from the heap.
    """
    python_string = cast(c_string, c_char_p).value.decode()
    C_LIB.freeCString(c_string)
    return python_string
