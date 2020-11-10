        program dimtest
        real a(1577,27), b(27,32), c(1577,32)
        real s, f
        a=1
        b=1
        call cpu_time(s)
        c = matmul(a,b)
        call cpu_time(f)
        print '("matmul Time = ",f6.3," seconds.")', (f-s)
        print '("matmul C[0,0] = ",f6.3," .")', (c(1,1))
        stop
        end
