        program dimtest
        !real a(1577,27), b(27,32), c(1577,32)
        !real a(3200,3100), b(3100,3200), c(3200,3200)
        real a(32,31), b(31,32), c(32,32)
        real s, f
        integer i
        a=1
        b=1
        c=1
        call cpu_time(s)
        !do i=1, 1000
        !do i=1, 7
        do i=1, 7007702
        c = c + matmul(a,b)
        end do
        call cpu_time(f)
        print '("matmul Time = ",f12.3," seconds.")', (f-s)
        print '("matmul C[0,0] = ",f12.3," .")', (c(1,1))
        stop
        end
