program fcn2
    implicit none

    integer, parameter :: block_size = 32
    integer, parameter :: matrix_size = 4096
    double precision, dimension(:, :), allocatable :: x, y, res_1, res_2
    real :: start, finish
    integer :: i, j
    allocate (x(matrix_size, matrix_size), y(matrix_size, matrix_size), res_2(matrix_size, matrix_size))

    ! x = reshape((/ 1, 5, -4, -8, -10, 4, 5, 6, 1, 3, 5, 8, 7, 9, -4, 2 /), shape(x))
    ! y = reshape((/ 4, 5, -4, -8, -10, 4, 6, 6, 1, 3, 2, 8, 7, 10, -4, 2 /), shape(y))
    !
    call random_number(x)
    call random_number(y)

    call cpu_time(start)
    res_1 = mtx_mult(x, y, matrix_size, block_size)
    call cpu_time(finish)

    print '("Time = ",f6.3," seconds.")',finish-start
    !write (*, *) ((res_1(i, j), i = 1, matrix_size), j = 1, matrix_size)

contains

function mtx_sum(mtx_A, mtx_B, split_index, multiplier) result(res)
    implicit none

    integer, optional,value :: multiplier
    integer, intent(in) :: split_index

    double precision, allocatable, dimension(:, :), intent(in) :: mtx_B, mtx_A

    double precision, allocatable, dimension(:, :) :: res
    integer :: i, j

    allocate(res(split_index, split_index))

    if (.NOT. PRESENT(multiplier)) multiplier = 1

    forall(i = 1 : block_size, j = 1 : block_size)
            res(i, j) = mtx_A(i, j) + (multiplier * mtx_B(i, j))
    end forall

end function mtx_sum


recursive function mtx_mult(mtx_A, mtx_B, n, block_size) result(result_matrix)
    implicit none

    integer, intent(in) :: block_size, n
    double precision :: ALPHA, BETA

    double precision, allocatable, dimension(:, :), intent(in) :: mtx_A, mtx_B

    double precision, allocatable, dimension(:, :) :: result_matrix

    integer :: split_index

    integer :: i, j

    double precision,allocatable, dimension(:, :) :: a_11, a_12, a_21, a_22, b_11, b_12, b_21, b_22, D, D_1, D_2, H_2, H_1, V_1, V_2
    double precision,allocatable, dimension(:, :) :: result_matrix_00, result_matrix_01, result_matrix_10, result_matrix_11

    allocate(result_matrix(n, n), result_matrix_11(n / 2, n / 2))

    allocate(result_matrix_00(n / 2, n / 2), result_matrix_01(n / 2, n / 2), result_matrix_10(n / 2, n / 2))

    ALPHA = 1.0
    BETA = 0.0

    if (n <= 1024) then
        call dgemm('N', 'N', n, n, n, ALPHA, mtx_A, n, mtx_B, n, BETA, result_matrix, n);
        return

    else
        split_index = n / 2

        a_11 = mtx_A(1:split_index, 1:split_index)
        a_12 = mtx_A(1:split_index, 1 + split_index : 2*split_index)
        a_21 = mtx_A(split_index + 1 : 2 * split_index, 1:split_index)
        a_22 = mtx_A(1 + split_index : 2 * split_index, 1 + split_index : 2*split_index)
        b_11 = mtx_B(1:split_index, 1:split_index)
        b_12 = mtx_B(1:split_index, 1 + split_index : 2*split_index)
        b_21 = mtx_B(split_index + 1 : 2 * split_index, 1:split_index)
        b_22 = mtx_B(1 + split_index : 2 * split_index, 1 + split_index : 2*split_index)


        D = mtx_mult(mtx_sum(a_11, a_22, split_index), mtx_sum(b_11, b_22, split_index), split_index, block_size)
        D_1 = mtx_mult(mtx_sum(a_12, a_22, split_index, -1), mtx_sum(b_21, b_22, split_index), split_index, block_size)
        D_2 = mtx_mult(mtx_sum(a_21, a_11, split_index, -1), mtx_sum(b_11, b_12, split_index), split_index, block_size)
        H_2 = mtx_mult(mtx_sum(a_21, a_22, split_index), b_11, split_index, block_size)
        H_1 = mtx_mult(mtx_sum(a_11, a_12, split_index), b_22, split_index, block_size)
        V_1 = mtx_mult(a_22, mtx_sum(b_21, b_11, split_index, -1),  split_index, block_size)
        V_2 = mtx_mult(a_11, mtx_sum(b_12, b_22, split_index, -1), split_index, block_size)

        result_matrix_00 = mtx_sum(mtx_sum(mtx_sum(D, D_1, split_index), V_1, split_index), H_1, split_index, -1)
        result_matrix_01 = mtx_sum(V_2, H_1, split_index)
        result_matrix_10 = mtx_sum(V_1, H_2, split_index)
        result_matrix_11 = mtx_sum(mtx_sum(mtx_sum(D, D_2, split_index), V_2, split_index), H_2, split_index, -1)

        result_matrix(1:block_size, 1:block_size) = result_matrix_00(1:block_size, 1:block_size)
        result_matrix(1:block_size, 1 + split_index:2 * split_index) = result_matrix_01(1:block_size, 1:block_size)
        result_matrix(split_index + 1:2 * split_index, 1:block_size) = result_matrix_10(1:block_size, 1:block_size)
        result_matrix(split_index + 1:2 * split_index, split_index + 1:2 * split_index) =&
        result_matrix_11(1:block_size, 1:block_size)


    end if

end function mtx_mult

end program fcn2
