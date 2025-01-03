# Quantum Trajectories and Metrology
[![](https://img.shields.io/badge/docs-maker?style=flat&color=blue&link=https%3A%2F%2Fste1nb0ck.github.io%2FQuTaM.jl%2Fdev%2Findex.html)](https://ste1nb0ck.github.io/QuTaM.jl/dev/index.html)

See the [tutorial](https://ste1nb0ck.github.io/QuTaM.jl/dev/tutorial.html) 
and the full [documentation](https://ste1nb0ck.github.io/QuTaM.jl/dev/index.html).

Recommended: in case you run into very long precompilation times each time you want to use the library [see](https://www.youtube.com/watch?v=_3vJSBk0Bls&t=15s).

For a more elaborated example look at `resonance_fluorescene.ipynb`.

## List of examples
You can find a list of examples in the `notebooks` directory. For the moment it contains
- _Resonance Fluorescene_(RF): at zero temperature 
- _Driven Qubit_(dq): Same model as RF, but at non zero temperature. At the moment it is incomplete, the comparison with the Lindblad equation is pending

## Tests
Automated tests to check the physics can be found in the `test` directory, the test the correct WTD statistics and the convergence to the 
Lindblad equation (among other things) for different models.
