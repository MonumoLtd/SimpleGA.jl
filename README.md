# SimpleGA

A simple, fast implementation of Geometric Algebra in Julia. The emphasis is on the geometric product. Inner and outer products are defined in terms of the geometric product.


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://monumoltd.github.io/SimpleGA.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://monumoltd.github.io/SimpleGA.jl/dev/)

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/MonumoLtd/SimpleGA.jl/blob/main/LICENSE)
[![Build Status](https://github.com/MonumoLtd/SimpleGA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/monumoltd/SimpleGA.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MonumoLtd/SimpleGA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MonumoLtd/SimpleGA.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

## Installation

Install:

```julia
julia>]
pkg> add SimpleGA
```

Start using the package. SimpleGA builds on some functionality in LinearAlgebra, so make sure you are using this.

```julia
using LinearAlgebra
using SimpleGA
```

## Example

Check out the documentation for more examples [check out the documentation.](https://monumoltd.github.io/SimpleGA.jl/dev/)

```julia
e = GA20.basis  # Creates a 2D basis, named 'e'.
a = e[1] + e[2]
b = 2.0*e[1] + 3*e[2] 
a*b  # Evaluates the geometric product of a and b.
```


