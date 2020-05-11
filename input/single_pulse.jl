function single_pulse(t,t1,amp,width)
    f = 0
    if (t>t1-width/2) & (t<t1+width/2)
        f = amp
    end
    return f
end
