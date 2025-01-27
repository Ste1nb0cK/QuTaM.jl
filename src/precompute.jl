############# Precomputation routine ######################
"""

```
setVs!(sys::System, nsamples::Int64, ts::Vector{Float64}, Vs::Array{ComplexF64})
````

Calculate the matrix exponentials ``\\mathrm{exp}(-iH_e t_s)`` for each ``t_s`` in the vector `ts`,
where ``H_e`` is the effective hamiltonian of the system `sys`, and the results are written
in `Vs`, which is an array of dimensions `(sys.NLEVELS, sys.NLEVELS, nsamples)`. To access
the exponential corresponding to `ts[k]` you would do `Vs[:, ;, k]`.
"""
function setVs!(sys::System, nsamples::Int64, ts::Vector{Float64}, Vs::Array{ComplexF64})
    tmp = copy(sys.Heff)
    @inbounds @simd for k in 1:nsamples
            Vs[:,:,k] = exp(-1im*ts[k]*sys.Heff)
    end
end


"""

```
setQs!(sys::System, nsamples::Int64,ts::Vector{Float64}, Qs::Array{ComplexF64}, Vs::Array{ComplexF64})
```

Calculate the matrix producs ``VJV^\\dagger`` for each ``V``
in  `Vs`, where  ``J=\\sum_k L_{k}^\\dagger L_k`` is  `sys.J`, and the results are written in `Qs` which is an array of dimensions
 `(sys.NLEVELS, sys.NLEVELS, nsamples)`. To access the product corresponding to `ts[k]` you would do `Qs[:, ;, k]`.
"""
function setQs!(sys::System, nsamples::Int64,
         ts::Vector{Float64}, Qs::Array{ComplexF64}, Vs::Array{ComplexF64})
        @inbounds @simd for k in 1:nsamples
            Qs[:, :, k] = Vs[:, :, k] * sys.J * adjoint(Vs[:, :, k])
        end
end


"""

```
precompute!(sys::System, nsamples::Int64, ts::Vector{Float64}, Qs::Array{ComplexF64}, Vs::Array{ComplexF64})
```

Precompute the ``Q(t_s)`` and ``V(t_s)`` necessary for running the *Quantum Gillipsie Algorithm*
 [radaelli2024gillespie](@cite) with the time grid `ts`. The result is written in `Qs` and `Vs`.
Under the hood, this simply calls `setVs!` and `setQs!`.
"""
function precompute!(sys::System, nsamples::Int64,
         ts::Vector{Float64}, Qs::Array{ComplexF64}, Vs::Array{ComplexF64})

    setVs!(sys, nsamples, ts, Vs)
    setQs!(sys, nsamples, ts, Qs, Vs)
end
