import numba
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


parent_pid = os.getpid()


@numba.jit(nopython=True)
def eval_i(m, n):
    return max(((len(m) // 4) // 3), ((len(n) // 4) // 3)) + 1


@numba.jit(nopython=True)
def divide(m, n):
    return str(int(m) // int(n))


def row_multiplication(row, col):
    results_private = ["0", "1"]

    for i in range(len(col)):
        tmp = (toom_cook(
            col[i], row[i][0]), row[i][1])
        num = str(int(toom_cook(
            tmp[0], results_private[1])) + int(toom_cook(tmp[1], results_private[0])))
        denom = toom_cook(tmp[1], results_private[1])

        results_private[0] = num
        results_private[1] = denom

    return str(int(results_private[0]) //
               int(results_private[1]))


def matmul(matrix, vector_, target=None):
    if target == "Parallel":
        with Pool(cpu_count()) as pool:
            result = list(pool.starmap(row_multiplication,
                                       zip(matrix, repeat(vector_))))
    else:
        result = [row_multiplication(row, vector_) for row in matrix]

    return result
    # return [row_multiplication(row, vector_) for row in matrix]


def toom_cook(m_, n_):
    m = m_
    n = n_

    if (len(m) <= 9 and len(n) <= 9) or (m[0] == '-' and n[0] == '-' and len(m) <= 10 and len(n) <= 10):
        return str(int(m) * int(n))

    m_minus = False
    n_minus = False

    if m[0] == '-':
        m = m[1:]
        m_minus = True

    if n[0] == '-':
        n = n[1:]
        n_minus = True

    B = 4 * eval_i(m, n)

    p = [i[::-1] for i in textwrap.wrap(m[::-1], B)]
    q = [i[::-1] for i in textwrap.wrap(n[::-1], B)]

    p = p + ["0"] * (3 - len(p))
    q = q + ["0"] * (3 - len(q))

    p_p = [p[0], str(int(p[0]) + int(p[1]) + int(p[2])), str(int(p[0]) - int(p[1]) + int(
        p[2])), str(int(p[0]) - int(toom_cook("2", p[1])) + int(toom_cook("4", p[2]))), p[2]]
    q_q = [q[0], str(int(q[0]) + int(q[1]) + int(q[2])), str(int(q[0]) - int(q[1]) + int(
        q[2])), str(int(q[0]) - int(toom_cook("2", q[1])) + int(toom_cook("4", q[2]))), q[2]]

    if parent_pid == os.getpid():
        with Pool(cpu_count()) as pool:
            r_in_points = list(pool.starmap(toom_cook, zip(p_p, q_q)))
    else:
        r_in_points = [toom_cook(p_p[0], q_q[0]), toom_cook(p_p[1], q_q[1]), toom_cook(
            p_p[2], q_q[2]), toom_cook(p_p[3], q_q[3]), toom_cook(p_p[4], q_q[4])]
    # r_in_points = [toom_cook(p_p[0], q_q[0]), toom_cook(p_p[1], q_q[1]), toom_cook(
    #         p_p[2], q_q[2]), toom_cook(p_p[3], q_q[3]), toom_cook(p_p[4], q_q[4])]

    if parent_pid == os.getpid():
        res = matmul(matrix, r_in_points, target="Parallel")
    else:
        res = matmul(matrix, r_in_points)

    recompos = [res[0], res[1] + "0" * B, res[2] + "0" *
                B * 2, res[3] + "0" * B * 3, res[4] + "0" * B * 4]

    if (m_minus and n_minus) or (not m_minus and not n_minus):
        return str(int(recompos[4]) + int(recompos[3]) + int(recompos[2]) + int(recompos[1]) + int(recompos[0]))

    return "-" + str(int(recompos[4]) + int(recompos[3]) + int(recompos[2]) + int(recompos[1]) + int(recompos[0]))


start = time.time()
tmp = toom_cook("1234567890123456782342323423442349023423423434534534512123131231233423234234423492342341231321231231231231377977895645",
                "9876543219823423424765432234534535342341023234423498123123131312334232342344234952342341231231231231231313178456478784")
end = time.time()
print(end - start)
