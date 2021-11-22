# Cartpole.jl

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

Written using Julia v1.6.3.

## Installing

### VSCode

(Install the Julia VSCode extension if you haven't already.)

Clone this repository and then open it in VSCode. To download the required dependencies, start a REPL (Shift-Ctrl p, Julia: Start REPL) and type in

```julia
using Pkg
Pkg.instantiate()
```

(See <https://youtu.be/CCfd8Y49xcU> for a quick introduction to using Julia projects in VSCode.)

### Julia REPL

Clone this repository and then start Julia. Type in

```julia
using Pkg
Pkg.activate("Path\To\Cloned\Project")
Pkg.instantiate()
```

## Usage

Afterward installing, simply type

```julia
using Cartpole
cartpole_example()
```

(The first run will be slow since everything needs to be compiled.)
