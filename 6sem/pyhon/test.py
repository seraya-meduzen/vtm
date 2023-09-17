import numba
from numba import njit, prange, cfunc, carray, types
import numpy as np
from itertools import repeat
from multiprocessing import Pool, cpu_count
import textwrap
import os
import time
from ctypes import cdll
import ctypes
from ctypes import *

lib = ctypes.CDLL("./libfoo.so")


lib.pyadd.argtypes = [ctypes.POINTER(
    c_char), ctypes.POINTER(c_char), ctypes.POINTER(c_char)]
lib.pyadd.restype = None


lhs = "11666"
rhs = "123132"

lhs = np.array([*lhs])
rhs = np.array([*rhs])
res = np.array(["0"] * (lhs.size + rhs.size))

lhs_ptr = lhs.ctypes.data_as(ctypes.POINTER(ctypes.c_char))
rhs_ptr = rhs.ctypes.data_as(ctypes.POINTER(ctypes.c_char))
res_ptr = res.ctypes.data_as(ctypes.POINTER(ctypes.c_char))

lib.pyadd(lhs_ptr, rhs_ptr, res_ptr)
print(res)
