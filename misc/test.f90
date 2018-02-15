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

      real, allocatable :: a(:), b(:)
      integer :: i
      call SSC_MARK_STOP();
      allocate(a(1:1024))
      allocate(b(1:1024))
      do i=1, 1024
        a(i) = 2
        b(i) = 3
      end do
      do i=1, 1024
        a(i) = a(i)*b(i)
      end do
      do i=1, 1024
        b(i) = 4
      end do
      call SSC_MARK_START();
      do i=1, 1024
        a(i) = a(i)*b(i)
      end do
      call SSC_MARK_STOP();
      print *, "Hello World!", a(1),a(1024)
      deallocate(a, b)
      call SSC_MARK_START();
end program hello
