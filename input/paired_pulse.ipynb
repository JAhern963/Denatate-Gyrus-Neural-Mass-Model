function paired_pulse(t,t1,T,amp,width,tmax)
    # function that creates a paired pulse current input as a time dependent
    # function.

    # t1    = start time of first pulse (s)
    # T     = time interval between pulses, 1/f, (s)
    # amp   = amplitude of a pulse
    # width = width of a pulse
    # tmax  = end time of the ODE solution

    t2 = t1 + T
    ti = [t1 t2]
    #loop to prevent pulse overlap and pulses outside time window
    if (ti[2]-width/2) <= (ti[1]+width/2)
        throw(ArgumentError("Pulses overlap: Reconsider arguments."))
    elseif ti[2]+width/2 >= tmax
        throw(ArgumentError("Input is not within time window."))
    end
    # loop to make the pulses
    f = 0
    for i = 1:2
        if (t>ti[i]-width/2) & (t<ti[i]+width/2)
            f = amp
        end
    end
    return f
end
