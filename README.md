# Polyhedron.jl

A Julia library for the computation and analysis of invariant sets in discrete-time linear dynamical systems.

`Polyhedron.jl` provides tools to define dynamical systems, represent state and input constraints, verify positive invariance, and compute invariant polyhedra using mathematical optimization and polyhedral computation.

The library is built around the Julia ecosystem, integrating optimization modeling and polyhedral manipulation tools.


## Installation

Since this package has not yet been published to the official Julia General Registry, it must be installed manually using Julia's `Pkg` mode.

### 1. Clone the repository

```bash
git clone https://github.com/Nicholas-Wallace/Polyhedron.git
```

### 2. Navigate to the cloned directory

```bash
cd Polyhedron
```

### 3. Start Julia

```bash
julia
```

### 4. Enter Pkg mode

In the Julia REPL, press `]` to enter package management mode. The prompt will change from `julia>` to `(@v1.x) pkg>`.

### 5. Add the package in development mode

```julia
dev path/to/Polyhedron
```

> Replace `path/to/Polyhedron` with the absolute or relative path to the cloned repository directory.

# Examples

The following examples illustrate how to define dynamical systems and use the library to investigate invariant polyhedra.

## Example 1.1 — Linear System with State Feedback

This example demonstrates how to define a system using its state-space matrices, specify state and input constraints, and investigate the invariance of a polyhedral set.

### Julia implementation

```julia
# ============================================================
# Example 1: Linear System with State Feedback
# ============================================================

# 1. Define the state-space matrices

A = [
    1.0  1.0
    0.0  1.0
]

B = [
    2.0
    1.0
]

C = [1.0  0.0]


# 2. Define the dynamical system
sys = SLIT(A, B, C, SSF)


# 3. Define state and input constraints
X = [0.8 0;
          0  1;
         -1  0;
          0  -1]
    
U = [1.2;
    -1.5;;]


# 4. Verify or compute an invariant polyhedron
ip = is_pinvariant(sys, U, X)
```

### Expected workflow

1. Define the system matrices.
2. Instantiate the corresponding system object.
3. Specify the state and input constraints.
4. Run the invariant set verification or computation procedure.
5. Inspect the resulting polyhedron and optimization status.

## Example 1.2 — Finding a Invariant Polyhedron

Use the same system, but with SOF, right now

```Julia
set_feedback!(sys, SOF)
```

We can verify that Positive Invariance cannot be acchived by the set defined by $X$, but we can find a polyhedron inside $X$ that can be Positive Invariant

```Julia
ip = finding_L_pinvariant(sys, U, X)
```

---

## Example 2 — Discrete-Time System with State and Input Constraints

Another example

### Julia implementation

```julia
# ============================================================
# Example 2: Constrained Linear System
# ============================================================

# 1. Define the state-space matrices

A = [
    0.9666  0.0
    0.0328  0.9666
]

B = [
    0.1980
    0.0033
]

C = [0.0  1.0]


# 2. Define the input constraints

U = [
     1/4
    -1/4
]


# 3. Define the state constraints

X = [
     1/15   0
    -1/15   0
     0      1/15
     0     -1/15
]


# 4. Define the dynamical system
system = SLIT(A, B, C, SSF)


# 5. Verify or compute an invariant polyhedron
ip = finding_L_pinvariant(system, U, X)
```

### Expected workflow

1. Define the system matrices.
2. Specify the state and input constraint matrices.
3. Instantiate the system object.
4. Execute the desired invariant set computation or verification.
5. Analyze the resulting polyhedron and, if applicable, visualize the results.

---

## Documentation and Further Examples

For more information about the library, its implementation, and its supported functionality, explore the source code and examples available in this repository.

Additional examples can be added here as the library evolves.

## Dependencies

The library builds on the Julia ecosystem for mathematical optimization and polyhedral computation. Its implementation uses tools including:

- [JuMP](https://jump.dev/) — mathematical optimization modeling.
- [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl) — polyhedral representation and manipulation.
- [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl) — polyhedral computation through CDDlib.

## Contributing

Contributions, bug reports, and suggestions are welcome. Feel free to open an issue or submit a pull request.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
