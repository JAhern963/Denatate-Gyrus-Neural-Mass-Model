include("single_pulse.jl")

function high_freq_protocol(t,t1,f,amp,width,n)

    if n == 1
        f = 0
        #throw(ArgumentError("Can't have n = 1, use 'single_pulse.jl'"))
        if (t>t1-width/2) & (t<t1+width/2)
            f = amp
        end
    else

        #make pulse times
        ti = zeros(1,n)
        ti[1] = t1
        for i = 2:n
            ti[i] = t1 + (i-1)*(1/f)
        end

        #loop to prevent pulse overlap and pulses outside time window
        if (ti[2]-width/2) <= (ti[1]+width/2)
            throw(ArgumentError("Pulses overlap: Reconsider arguments WIDTH and F."))
        end

        #make the pulses
        f = 0
        for i = 1:n
            if (t>ti[i]-width/2) & (t<ti[i]+width/2)
                f = amp
            end
        end
        tend = ti[n] + 3*(1/f)
    end 

    return f
end
