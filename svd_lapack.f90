program fcn2
    implicit None
    INTEGER :: t1, t2, count_rate, count_max
    integer :: i, info, lda, ldu, ldvt, lwork, m, n

    double precision, allocatable :: a(:, :), a_copy(:, :), b(:), s(:), u(:, :), vt(:, :), work(:)
    double precision :: dummy(1, 1)

    double precision, allocatable :: iwork(:)

    m = 2048
    n = m

    lda = m
    ldu = m
    ldvt = n

    allocate (a(lda,n), s(n), vt(ldvt,n), u(ldu,m), iwork(8 * m))

    call random_number(a)

    lwork = -1
    CALL SYSTEM_CLOCK(t1, count_rate, count_max)
    call dgesdd('A', m, n, a, lda, s, u, ldu, vt, ldvt, dummy, lwork, iwork, info)

    lwork = dummy(1,1)

    allocate (work(lwork))

    call dgesdd('A', m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
    call system_clock (t2, count_rate, count_max )
    write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real ( count_rate )
    !print *, u
    !print *, s
    !print *, vt

end program fcn2
