program hello

      use, intrinsic :: iso_c_binding
      interface
      integer(C_INT) function SSC_MARK_START ( i,j ) bind ( c )
      use, intrinsic :: iso_c_binding
      integer(C_INT), value :: i,j
      !!!DEC$ ATTRIBUTES C:: SSC_MARK_START
      end function
      integer(C_INT) function SSC_MARK_STOP ( i,j ) bind ( c )
      use, intrinsic :: iso_c_binding
      integer(C_INT),value :: i,j
      !!!DEC$ ATTRIBUTES C:: SSC_MARK_STOP
      end function
      end interface

      include 'mpif.h'
      real, allocatable :: a(:), b(:)
      integer :: i, rum=100000000, rank, ierror

      real(8) :: mkrts, mkrte
      integer(C_INT) :: sscfake = 1

      call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      sscfake = SSC_MARK_STOP(0,rank)
      allocate(a(1:rum))
      allocate(b(1:rum))
      !$OMP PARALLEL DO
      do i=1, rum
        a(i) = 2
        b(i) = 3
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      do i=1, rum
        a(i) = a(i)*b(i)
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      do i=1, rum
        b(i) = 4
      end do
      !$OMP END PARALLEL DO

      mkrts = MPI_WTIME()
      sscfake = SSC_MARK_START(1,rank)

      !$OMP PARALLEL DO
      do i=1, rum
        a(i) = a(i)*b(i)
      end do
      !$OMP END PARALLEL DO

      sscfake = SSC_MARK_STOP(1,rank)
      mkrte = MPI_WTIME()
      write(*,"(A,F10.6,A,I)") "Walltime of the main kernel: ", mkrte - mkrts, " sec", rank

      write(*,*) "Hello World!", a(1), a(rum)
      deallocate(a, b)
      call MPI_FINALIZE(ierror)

      sscfake = SSC_MARK_START(0,rank)
      write(*,*) sscfake
end program hello
