############# Precomputation routine ######################
"""

    precompute!(sys::System, nsamples::Int64,
         ts::Vector{Float64}, Qs::Vector{Matrix{ComplexF64}})

Does the precomputation routine for the Gillipsie algorithm.
The result is stored in the `Qs`.
# Arguments

`sys::System`: system
`nsamples::Int64`: number of samples
`ts::Vector{Float64}`: fine grid vector. *IT IS ASSUMED HOMOGENOUS*
`Qs::Vector{Matrix{ComplexF64}}`: array of matrices to store the precomputation,
                                  dimensions must be (sys.NLEVELS, sys.NLEVELS, nsamples)
"""
function precompute!(sys::System, nsamples::Int64,
         ts::Vector{Float64}, Qs::Array{ComplexF64}, Vs::Array{ComplexF64}; progbar::Bool=true)
        @inbounds @simd for k in 1:nsamples
            # expm = exp(-1im*ts[k]*sys.Heff)
            Vs[:,:,k] = exp(-1im*ts[k]*sys.Heff)
            Qs[:, :, k] = Vs[:, :, k] * sys.J * adjoint(Vs[:, :, k])
        end
end
