# BackAction.jl
Yet another library for Open Quantum Systems and Quantum Information
[![](https://img.shields.io/badge/docs-maker?style=flat&color=blue&link=https%3A%2F%2Fste1nb0ck.github.io%2FQuTaM.jl%2Fdev%2Findex.html)](https://ste1nb0ck.github.io/QuTaM.jl/dev/index.html)

See the [tutorial](https://ste1nb0ck.github.io/QuTaM.jl/dev/tutorial.html) 
and the full [documentation](https://ste1nb0ck.github.io/QuTaM.jl/dev/index.html).

Recommended: in case you run into very long precompilation times each time you want to use the library [see](https://www.youtube.com/watch?v=_3vJSBk0Bls&t=15s).

For a more elaborated example look at `resonance_fluorescene.ipynb`.
## Installation
## Installation
For the moment the library is not available through Pkg, so you will need to clone the repository 

```console
$ git clone https://github.com/Ste1nb0cK/BackAction.jl
```

and build it with `] develop path_to_clonedrepo` in the Julia REPL, this automatically adds it to your
current project. 

```julia
import Pkg; Pkg.develop(path="./BackAction.jl")
```

To check the installation do

```julia
Pkg.status
```

and you should see `BackAction` in the output e.g.

```console
  [32f8aca8] BackAction v0.1.0 `..`
  [e30172f5] Documenter v1.8.0
  [daee34ce] DocumenterCitations v1.3.5
  [b964fa9f] LaTeXStrings v1.4.0
  [91a5bcdd] Plots v1.40.9
  [10745b16] Statistics v1.10.0
  [8dfed614] Test

```

## List of examples
You can find a list of examples in the `notebooks` directory. It contains for the moment:

- _Resonance Fluorescene_(RF): at zero temperature 
- _Driven Qubit_(dq): Same model as RF, but at non zero temperature. 

## Tests

Automated tests to check the physics can be found in the `test` directory, these check the correct WTD statistics and the convergence to the 
Lindblad equation (among other things) for different models.
