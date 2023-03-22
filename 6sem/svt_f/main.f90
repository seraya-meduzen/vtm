program main
    implicit none

    real*8, dimension (:), allocatable :: a, rhs, alu, w, sol, vv
    integer, dimension (:), allocatable :: ia, ja, jlu, ju, levs, jw
    integer :: nmax, nzmax , n, nnz, ierr, iounit, i, lfil, iwk, maxits, im, iout
    real :: start, finish
    real*8 :: eps

    nmax = 264752
    nzmax = 3541545
    iounit = 1

    allocate(ia(nmax+1), ja(nzmax), a(nzmax))
    open(iounit, file="./../matr/4.dat", status ="old")

    call readsk(nmax, nzmax, n, nnz, a, ja, ia, 1, ierr)

    print '("Matrix Size: ",i0)',n
    if(ierr /= 0) then
        print '("read problems")'
        stop ierr
    end if

    allocate(rhs(n))
    do i = 1, n
        rhs(i) = sin(real(i))
    end do

    iwk = nnz * 100
    allocate(jw(3 * n), w(n), alu(iwk), jlu(iwk), levs(iwk), ju(n))
    lfil = 0

    call cpu_time(start)
    call iluk(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr)
    call cpu_time(finish)
    print '("Time ILU = ",f6.3," seconds.")',finish-start
    if(ierr /= 0) then
        print '("iluk failed")'
        stop ierr
    end if

    eps = 1e-8
    im = 10
    maxits = n * n
    iout = 2

    open(iout, file="../log.txt", status ="old")
    allocate(sol(n), vv(n * (im + 1)))

    call cpu_time(start)
    call pgmres(n, im, rhs, sol, vv, eps, maxits, iout, a, ja, ia, alu, jlu, ju, ierr)
    call cpu_time(finish)

    print '("Time for pgmres = ",f9.3," seconds.")',finish-start
    if(ierr /= 0) then
        print '("pgmres failed")'
        stop ierr
    end if


    deallocate(ia, ja, a, jw, w, alu, jlu, levs, ju, sol, vv, rhs)

end program main
