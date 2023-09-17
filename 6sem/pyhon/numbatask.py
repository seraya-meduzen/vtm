import numba
from numba import njit, prange
import numpy as np
from itertools import repeat
from multiprocessing import Pool, cpu_count
import textwrap
import os
import time


K = 3

matrix = [[("1", "1"), ("0", "1"), ("0", "1"), ("0", "1"), ("0", "1")],
          [("1", "2"), ("1", "3"), ("-1", "1"), ("1", "6"), ("-2", "1")],
          [("-1", "1"), ("1", "2"), ("1", "2"), ("0", "1"), ("-1", "1")],
          [("-1", "2"), ("1", "6"), ("1", "2"), ("-1", "6"), ("2", "1")],
          [("0", "1"), ("0", "1"), ("0", "1"), ("0", "1"), ("1", "1")]]


# parent_pid = os.getpid()


@numba.jit(nopython=True)
def eval_i(m, n):
    return max(((len(m) // 4) // 3), ((len(n) // 4) // 3)) + 1


@numba.jit(nopython=True)
def divide(m, n):
    return str(int(m) // int(n))


# @njit
# def row_multiplication(row, col):
#     results_private = ["0", "1"]

#     for i in range(len(col)):
#         tmp = (toom_cook(
#             col[i], row[i][0]), row[i][1])
#         num = str(int(toom_cook(
#             tmp[0], results_private[1])) + int(toom_cook(tmp[1], results_private[0])))
#         denom = toom_cook(tmp[1], results_private[1])

#         results_private[0] = num
#         results_private[1] = denom

#     # return str(int(results_private[0]) //
#     #            int(results_private[1]))

#     return results_private


# @njit(parallel=True)
# def matmul(matrix, vector_):
#     result = ["0"] * len(matrix)
#     for row in prange(len(matrix)):
#         result[row] = row_multiplication(matrix[row], vector_)

#     return result

# @njit
# def toom_cook(m_, n_):
#     m = m_
#     n = n_

#     if (len(m) <= 9 and len(n) <= 9) or (m[0] == '-' and n[0] == '-' and len(m) <= 10 and len(n) <= 10):
#         # return str(int(m) * int(n))
#         return str(n)

    # m_minus = False
    # n_minus = False

    # if m[0] == '-':
    #     m = m[1:]
    #     m_minus = True

    # if n[0] == '-':
    #     n = n[1:]
    #     n_minus = True

    # B = 4 * eval_i(m, n)

    # p = [i[::-1] for i in textwrap.wrap(m[::-1], B)]
    # q = [i[::-1] for i in textwrap.wrap(n[::-1], B)]

    # p = p + ["0"] * (3 - len(p))
    # q = q + ["0"] * (3 - len(q))

    # p_p = [p[0], str(int(p[0]) + int(p[1]) + int(p[2])), str(int(p[0]) - int(p[1]) + int(
    #     p[2])), str(int(p[0]) - int(toom_cook("2", p[1])) + int(toom_cook("4", p[2]))), p[2]]
    # q_q = [q[0], str(int(q[0]) + int(q[1]) + int(q[2])), str(int(q[0]) - int(q[1]) + int(
    #     q[2])), str(int(q[0]) - int(toom_cook("2", q[1])) + int(toom_cook("4", q[2]))), q[2]]

    # r_in_points = [toom_cook(p_p[0], q_q[0]), toom_cook(p_p[1], q_q[1]), toom_cook(
    #     p_p[2], q_q[2]), toom_cook(p_p[3], q_q[3]), toom_cook(p_p[4], q_q[4])]

    # res = list(map(str, map(int, matmul(matrix, np.array(
    #     list(map(int, r_in_points)), dtype=np.float64)))))

    # recompos = [res[0], res[1] + "0" * B, res[2] + "0" *
    #             B * 2, res[3] + "0" * B * 3, res[4] + "0" * B * 4]

    # if (m_minus and n_minus) or (not m_minus and not n_minus):
    #     return str(int(recompos[4]) + int(recompos[3]) + int(recompos[2]) + int(recompos[1]) + int(recompos[0]))

    # return "-" + str(int(recompos[4]) + int(recompos[3]) + int(recompos[2]) + int(recompos[1]) + int(recompos[0]))

@numba.jit(nopython=True)
def multiply_strings(num1, num2):
    if num1 == '0' or num2 == '0':
        return '0'

    result = [0] * (len(num1) + len(num2))

    num1 = num1[::-1]
    num2 = num2[::-1]

    for i in range(len(num1)):
        for j in range(len(num2)):
            digit1 = ord(num1[i]) - ord("0")
            digit2 = ord(num2[j]) - ord("0")
            product = digit1 * digit2

            result[i + j] += product % 10
            result[i + j + 1] += product // 10

            if result[i + j] >= 10:
                result[i + j + 1] += result[i + j] // 10
                result[i + j] %= 10

    while len(result) > 1 and result[-1] == 0:
        result.pop()

    return ''.join(map(str, result[::-1]))


start = time.time()
tmp = multiply_strings("11564",
                       "12644")
end = time.time()
print(end - start, tmp)
