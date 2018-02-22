program hello

      use, intrinsic :: iso_c_binding
      interface
      subroutine SSC_MARK_START ( ) bind ( c )
      !!!DEC$ ATTRIBUTES C:: SSC_MARK_START
      end subroutine
      subroutine SSC_MARK_STOP ( ) bind ( c )
      !!!DEC$ ATTRIBUTES C:: SSC_MARK_STOP
      end subroutine
      end interface

      include 'mpif.h'
      real, allocatable :: a(:), b(:)
      integer :: i, ierror

      real(8) :: mkrts, mkrte
      call SSC_MARK_STOP()

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
      call SSC_MARK_START()

      !$OMP PARALLEL DO
      do i=1, 1024
        a(i) = a(i)*b(i)
      end do
      !$OMP END PARALLEL DO

      call SSC_MARK_STOP()
      mkrte = MPI_WTIME()
      write(*,"(A,F10.6,A)") "Walltime of the main kernel: ", mkrte - mkrts, " sec"

      write(*,*) "Hello World!", a(1), a(1024)
      deallocate(a, b)
      call MPI_FINALIZE(ierror)

      call SSC_MARK_START()
end program hello
