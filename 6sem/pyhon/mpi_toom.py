from mpi4py import MPI
import textwrap
import time

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

K = 3

matrix = [[("1", "1"), ("0", "1"), ("0", "1"), ("0", "1"), ("0", "1")],
          [("1", "2"), ("1", "3"), ("-1", "1"), ("1", "6"), ("-2", "1")],
          [("-1", "1"), ("1", "2"), ("1", "2"), ("0", "1"), ("-1", "1")],
          [("-1", "2"), ("1", "6"), ("1", "2"), ("-1", "6"), ("2", "1")],
          [("0", "1"), ("0", "1"), ("0", "1"), ("0", "1"), ("1", "1")]]


def eval_i(m, n):
    return max(((len(m) // 4) // 3), ((len(n) // 4) // 3)) + 1


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


def matmul_nonparallel(matrix, vector_):
    return [row_multiplication(row, vector_) for row in matrix]


def matmul_parallel(matrix, vector_):
    chunk_size = len(matrix) // size
    start_row = rank * chunk_size
    end_row = start_row + chunk_size

    if (rank == size - 1):
        end_row += len(matrix) % size

    local_result = []
    for j in range(start_row, end_row):

        results_private = ["0", "1"]

        for i in range(len(vector_)):
            tmp = (toom_cook(vector_[i], matrix[j][i][0]), matrix[j][i][1])
            num = str(int(toom_cook(tmp[0], results_private[1])) +
                      int(toom_cook(tmp[1], results_private[0])))
            denom = toom_cook(tmp[1], results_private[1])

            results_private[0] = num
            results_private[1] = denom

        local_result.append(
            str(int(results_private[0]) // int(results_private[1])))

    all_results = comm.gather(local_result, root=0)

    if rank == 0:
        result = []
        for sublist in all_results:
            result.extend(sublist)
    else:
        result = None

    result = comm.bcast(result, root=0)
    return result


def toom_cook(m_, n_, matmul=matmul_nonparallel):
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

    r_in_points = [toom_cook(p_p[0], q_q[0]), toom_cook(p_p[1], q_q[1]), toom_cook(
        p_p[2], q_q[2]), toom_cook(p_p[3], q_q[3]), toom_cook(p_p[4], q_q[4])]

    res = matmul(matrix, r_in_points)

    recompos = [res[0], res[1] + "0" * B, res[2] + "0" *
                B * 2, res[3] + "0" * B * 3, res[4] + "0" * B * 4]

    if (m_minus and n_minus) or (not m_minus and not n_minus):
        return str(int(recompos[4]) + int(recompos[3]) + int(recompos[2]) + int(recompos[1]) + int(recompos[0]))

    return "-" + str(int(recompos[4]) + int(recompos[3]) + int(recompos[2]) + int(recompos[1]) + int(recompos[0]))


if __name__ == '__main__':
    start = time.time()
    tests = [toom_cook("1234567890123456782342323423442349023423423434534534512123131231233423234234423492342341231321231231231231377977895645",
                       "9876543219823423424765432234534535342341023234423498123123131312334232342344234952342341231231231231231313178456478784", matmul_parallel)]
    end = time.time()

    if rank == 0:
        # print(*tests, sep="\n")
        print(end - start)
