function Sigmoid(x,e0,r,v0)
    #Sigmoid function for translating the postsynaptic activity of a
    #neuron (x) into a firing rate
    S = 2*e0 ./ (1 .+ exp.(r.*(v0 .- x)))
    return S
end
