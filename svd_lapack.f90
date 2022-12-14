program fcn2
    implicit None

    integer :: i, info, lda, ldu, ldvt, lwork, m, n

    double precision, allocatable :: a(:, :), a_copy(:, :), b(:), s(:), u(:, :), vt(:, :), work(:)
    double precision :: dummy(1, 1)

    double precision, allocatable :: iwork(:)

    m = 4
    n = 4

    lda = m
    ldu = m
    ldvt = n

    allocate (a(lda,n), s(n), vt(ldvt,n), u(ldu,m), iwork(8 * m))

    call random_number(a)

    lwork = -1

    call dgesdd('A', m, n, a, lda, s, u, ldu, vt, ldvt, dummy, lwork, iwork, info)

    lwork = dummy(1,1)

    allocate (work(lwork))

    call dgesdd('A', m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)

    print *, u
    print *, s
    print *, vt

end program fcn2
