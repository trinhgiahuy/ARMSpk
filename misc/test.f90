program hello

      use, intrinsic :: iso_c_binding
      interface
      integer(4) function SSC_MARK_START ( i ) bind ( c )
      integer(4) :: i
      !!!DEC$ ATTRIBUTES C:: SSC_MARK_START
      end function
      integer(4) function SSC_MARK_STOP ( j ) bind ( c )
      integer(4) :: j
      !!!DEC$ ATTRIBUTES C:: SSC_MARK_STOP
      end function
      end interface

      include 'mpif.h'
      real, allocatable :: a(:), b(:)
      integer :: i, ierror

      real(8) :: mkrts, mkrte
      integer(4) :: sscfake = 1
      sscfake = SSC_MARK_STOP(sscfake)

      call MPI_INIT(ierror)
      allocate(a(1:1024))
      allocate(b(1:1024))
      !$OMP PARALLEL DO
      do i=1, 1024
        a(i) = 2
        b(i) = 3
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      do i=1, 1024
        a(i) = a(i)*b(i)
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      do i=1, 1024
        b(i) = 4
      end do
      !$OMP END PARALLEL DO

      mkrts = MPI_WTIME()
      sscfake = SSC_MARK_START(sscfake)

      !$OMP PARALLEL DO
      do i=1, 1024
        a(i) = a(i)*b(i)
      end do
      !$OMP END PARALLEL DO

      sscfake = SSC_MARK_STOP(sscfake)
      mkrte = MPI_WTIME()
      write(*,"(A,F10.6,A)") "Walltime of the main kernel: ", mkrte - mkrts, " sec"

      write(*,*) "Hello World!", a(1), a(1024)
      deallocate(a, b)
      call MPI_FINALIZE(ierror)

      sscfake = SSC_MARK_START(sscfake)
      write(*,*) sscfake
end program hello
