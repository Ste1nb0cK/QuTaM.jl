function samplechannel(Vs::Array{ComplexF64}, P::Vector{Float64},
                                    psi::Vector{ComplexF64}, sys::System, tau_index::Int64)
        # 3. Sample the channel
        channel = StatsBase.sample(1:sys.NCHANNELS, StatsBase.weights(P))
               return channel
end
