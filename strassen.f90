program fcn2
    implicit none

    integer, parameter :: block_size = 2
    integer, parameter :: matrix_size = 4

    double precision, dimension(:, :), allocatable :: x, y, res_1, res_2
    ! double precision, dimension(:), allocatable :: z

    integer :: i, j
    allocate (x(matrix_size, matrix_size), y(matrix_size, matrix_size))

    ! x = reshape((/ 1, 5, -4, -8, -10, 4, 5, 6, 1, 3, 5, 8, 7, 9, -4, 2 /), shape(x))
    ! y = reshape((/ 4, 5, -4, -8, -10, 4, 6, 6, 1, 3, 2, 8, 7, 10, -4, 2 /), shape(y))

    call random_number(x)
    call random_number(y)

    res_1 = mtx_mult(x, y, matrix_size, block_size)
    res_2 = matmul(x, y, matrix_size, block_size)

    write (*, *) ((res_1(i, j), i = 1, matrix_size), j = 1, matrix_size)
    write (*, *) ((res_2(i, j), i = 1, matrix_size), j = 1, matrix_size)

contains

subroutine multiply_block(X, Y, res, n, block_size)
    implicit none
    integer, intent(in) :: n, block_size

    double precision, allocatable, intent(inout) :: res(:, :)
    double precision, allocatable, intent(in) :: X(:, :)
    double precision, intent(in) :: Y(:, :)
    ! double precision, allocatable, intent(in) :: Y(:)

    integer :: i, j, s
    double precision :: tmp

    do i = 1, block_size
        do s = 1, block_size
            do j = 1, block_size
                tmp = X(s, i) * Y(j, s)
                res(j, i) = res(j, i) + tmp
            end do
        end do
    end do

end subroutine multiply_block


subroutine copy_a(block_a, A, n, block_size)
    implicit none

    integer, intent(in) :: n, block_size
    double precision, allocatable, dimension(:, :), intent(inout) :: block_a
    double precision, dimension(:, :), intent(in) :: A
    ! double precision, allocatable, dimension(:), intent(in) :: A

    integer :: i, j

    forall(i = 1 : block_size, j = 1 : block_size)
            block_a(j, i) = A(j, i);
    end forall

end subroutine copy_a
!
subroutine copy_c(block_c, C, n, block_size)
    implicit none

    integer, intent(in) :: n, block_size

    double precision, dimension(:, :), intent(inout) :: C
    ! double precision, allocatable, dimension(:), intent(inout) :: C
    double precision, allocatable, dimension(:, :), intent(in) :: block_c

    integer :: i, j

    forall(i = 1 : block_size, j = 1 : block_size)
            C(j, i) = block_c(j, i);
    end forall

end subroutine copy_c


function matmul(A, B, n, block_size) result(C)
    implicit none
    integer, intent(in) :: n, block_size

    double precision, allocatable, dimension(:, :), intent(in) :: A, B

    double precision, allocatable, dimension(:, :) :: C, block_c, block_a

    integer :: filler
    integer :: i, j, k

    allocate(C(n, n), block_c(block_size, block_size), block_a(block_size, block_size))

    do i = 1, n, block_size
        do j = 1, n, block_size
            block_c = reshape((/ (0, filler = 1, size(block_c)) /), shape(block_c))

            do k = 1, n, block_size
                call copy_a(block_a, A(k:, i:), n, block_size)
                call multiply_block(block_a, B(j:, k:), block_c, n, block_size)
            end do
            call copy_c(block_c, C(j:, i:), n, block_size);
        end do
    end do

end function matmul
!
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

    double precision, allocatable, dimension(:, :), intent(in) :: mtx_A, mtx_B

    double precision, allocatable, dimension(:, :) :: result_matrix

    integer :: split_index

    integer :: i, j

    double precision,allocatable, dimension(:, :) :: a_11, a_12, a_21, a_22, b_11, b_12, b_21, b_22, D, D_1, D_2, H_2, H_1, V_1, V_2
    double precision,allocatable, dimension(:, :) :: result_matrix_00, result_matrix_01, result_matrix_10, result_matrix_11

    allocate(result_matrix(n, n))

    allocate(a_11(n / 2, n / 2), a_12(n / 2, n / 2), a_21(n / 2, n / 2), a_22(n / 2, n / 2), b_11(n / 2, n / 2), b_12(n / 2, n / 2))
    allocate(b_21(n / 2, n / 2), b_22(n / 2, n / 2), D(n / 2, n / 2), D_1(n / 2, n / 2), D_2(n / 2, n / 2), H_2(n / 2, n / 2))
    allocate(H_1(n / 2, n / 2), V_1(n / 2, n / 2), V_2(n / 2, n / 2), result_matrix_11(n / 2, n / 2))

    allocate(result_matrix_00(n / 2, n / 2), result_matrix_01(n / 2, n / 2), result_matrix_10(n / 2, n / 2))


    if (n <= 1024) then
        result_matrix = matmul(mtx_B, mtx_A, n, block_size)
        return
    end if

    if (n == 1) then
        result_matrix(1, 1) = mtx_A(1, 1) * mtx_B(1, 1)
    else
        split_index = n / 2

        forall(i = 1 : block_size, j = 1 : block_size)
                a_11(i, j) = mtx_B(i, j)
				a_12(i, j) = mtx_B(i, j + split_index)
				a_21(i, j) = mtx_B(split_index + i, j)
				a_22(i, j) = mtx_B(i + split_index, j + split_index)
				b_11(i, j) = mtx_A(i, j)
				b_12(i, j) = mtx_A(i, j + split_index)
				b_21(i, j) = mtx_A(split_index + i, j)
				b_22(i, j) = mtx_A(i + split_index, j + split_index)
        end forall

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

        forall(i = 1 : block_size, j = 1 : block_size)
                result_matrix(i, j) = result_matrix_00(i, j)
				result_matrix(i, j + split_index) = result_matrix_01(i, j)
				result_matrix(split_index + i, j) = result_matrix_10(i, j)
				result_matrix(i + split_index, j + split_index) = result_matrix_11(i, j)
        end forall

    end if

end function mtx_mult

end program fcn2
