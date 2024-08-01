# Analysis of an Incompressible Flow Inside a Backward-Facing Step Using Artificial Compressibility Method

**Written By:** Seyed MohammadAmin Taleghani  
**Date:** July 2021

## Table of Contents
1. [Description of the Problem](#description-of-the-problem)
2. [Governing Equations](#governing-equations)
   - [Non-dimensionalization](#non-dimensionalization)
   - [Artificial Compressibility](#artificial-compressibility)
3. [Implicit Discretization](#implicit-discretization)
4. [Boundary Conditions](#boundary-conditions)
5. [Solution Algorithm](#solution-algorithm)
6. [Generated Grids](#generated-grids)
7. [Results](#results)
   - [Convergence History](#convergence-history)
   - [Wall Shear](#wall-shear)
   - [Pressure Contours](#pressure-contours)
   - [Streamlines and Velocity Vectors](#streamlines-and-velocity-vectors)
   - [Wake Length Separation and Reattachment Locations](#wake-length-separation-and-reattachment-locations)

## Description of the Problem
The flow past a backward-facing step is a good test case to examine the accuracy and performance of a numerical method. The geometry of this test case is simple yet contains complex flow features associated with separation and reattachment. The main difficulty is the high dependence of the recirculation regions on grid size and resolution.

## Governing Equations
Two-dimensional Cartesian coordinate dimensional nonconservative form of the Incompressible Navier-Stokes equations:

### Non-dimensionalization
The variables in the equations above are non-dimensionalized for simplicity.

### Artificial Compressibility
The continuity equation is modified by including a time-dependent term:

\[ \frac{\partial p}{\partial t} + \nabla \cdot \mathbf{v} = 0 \]

Where \( \beta \) is the artificial compressibility of the fluid. This can be related to a pseudo-speed of sound and an artificial density.

## Implicit Discretization
For the two-dimensional problems, the implicit formulation results in a block penta-diagonal system of equations. To make the solution process less expensive, approximate factorization and artificial viscosity are implemented.

## Boundary Conditions
### Left Boundary (i=2)
The left boundary condition controls the equation for \( U \). Since the velocity components are invariant at the left boundary, specific conditions for pressure are set using a second-order forward FDE.

### Lower and Upper Walls
Derivations for the lower and upper walls are provided, along with conditions at the outlet boundary (i = imax-1).

## Solution Algorithm
The solution algorithm involves debugging and final grids.

## Generated Grids
- **Debugging Grid 51x21**  
  ![Debugging Grid 51x21](images/debugging_grid.png)
- **Final Grid 141x51**  
  ![Final Grid 141x51](images/final_grid.png)

## Results

### Convergence History
- Convergence history of the solution based on the divergence of the velocity vector.
- ![Convergence History](images/convergence_history.png)

### Wall Shear
- Wall shear for both upper and lower walls.
- ![Wall Shear](images/wall_shear.png)

### Pressure Contours
Pressure contours at various Reynold numbers:
- Re = 800  
  ![Pressure Contours Re800](images/pressure_contours_re800.png)
- Re = 600  
  ![Pressure Contours Re600](images/pressure_contours_re600.png)
- Re = 400  
  ![Pressure Contours Re400](images/pressure_contours_re400.png)
- Re = 200  
  ![Pressure Contours Re200](images/pressure_contours_re200.png)
- Re = 80  
  ![Pressure Contours Re80](images/pressure_contours_re80.png)
- Re = 40  
  ![Pressure Contours Re40](images/pressure_contours_re40.png)
- Re = 1  
  ![Pressure Contours Re1](images/pressure_contours_re1.png)

### Streamlines and Velocity Vectors
Streamlines and velocity vectors at various Reynold numbers follow the pattern shown in the pressure contours section.

### Wake Length Separation and Reattachment Locations
Comparison of wake length, separation location, and reattachment location at different Reynold numbers against Hejranfar and Khajeh-Saeed's results.

---

This README file provides a structured overview of the project, making it easier for users to understand the problem, methodology, and results. For detailed equations, figures, and more in-depth explanations, refer to the full project report.
