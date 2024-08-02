# Analysis of an Incompressible Flow Inside a Backward-Facing Step Using Artificial Compressibility Method

(Based on a Problem in “Computational Fluid Dynamics for Engineers Vol. 1” by Klaus A. Hoffmann)

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
   - [Streamlines and Velocity Vectors](#streamlines-and-velocity-vectors)

## Description of the Problem
The flow past a backward-facing step is a good test case to examine the accuracy and performance of a numerical method. The geometry of this test case is simple yet contains complex flow features associated with separation and reattachment. The main difficulty is the high dependence of the recirculation regions on grid size and resolution.

![image](https://github.com/user-attachments/assets/611ed275-3927-45bf-9f1d-2607e7e3b055)


## Governing Equations
Two-dimensional Cartesian coordinate dimensional nonconservative form of the Incompressible Navier-Stokes equations:

$$ \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0 $$

$$ \frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} + \frac{1}{\rho} \frac{\partial p}{\partial x} = \nu \left\( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right\) $$

$$ \frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} + \frac{1}{\rho} \frac{\partial p}{\partial y} = \nu \left\( \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right\) $$

### Non-dimensionalization
The variables in the equations above are non-dimensionalized for simplicity.

$$ \frac{\partial u^\*}{\partial x^\*} + \frac{\partial v^\*}{\partial y^\*} = 0 $$

$$ \frac{\partial u^\*}{\partial t^\*} + u^\* \frac{\partial u^\*}{\partial x^\*} + v^\* \frac{\partial u^\*}{\partial y^\*} + \frac{\partial p^\*}{\partial x^\*} = \frac{1}{Re} \left( \frac{\partial^2 u^\*}{\partial x^{\*^2}} + \frac{\partial^2 u^\*}{\partial y^{\*^2}} \right) $$

$$ \frac{\partial v^\*}{\partial t^\*} + u^\* \frac{\partial v^\*}{\partial x^\*} + v^\* \frac{\partial v^\*}{\partial y^\*} + \frac{\partial p^\*}{\partial y^\*} = \frac{1}{Re} \left( \frac{\partial^2 v^\*}{\partial x^{\*^2}} + \frac{\partial^2 v^\*}{\partial y^{\*^2}} \right) $$

The variables in the equations above are nondimensionalized as follows:

$$ t^{\*} = \frac{t u_\infty}{L} ~~~~ x^{\*} = \frac{x}{L} ~~~~ y^{\*} = \frac{y}{L}$$

$$ u^{\*} = \frac{u}{u_\infty} ~~~~ v^{\*} = \frac{v}{u_\infty}  ~~~~ p^{\*} = \frac{p}{\rho_\infty u_\infty^2}$$

From now on, we shall drop the * from non-dimensional parameters for the sake of simplicity.


### Artificial Compressibility
The continuity equation is modified by including a time-dependent term:

$$ \frac{\partial p}{\partial t} + \nabla \cdot \mathbf{v} = 0 $$

Where $\beta$ can be interpreted as the “artificial compressibility” of the fluid. Following the equation of state, the compressibility can be related to a pseudo-speed of sound and to an artificial density by the following relations:

$$ p = \beta \tilde{\rho} = \tilde{a}^2 \tilde{\rho} $$

$$ \left\{ \begin{array}{l}
\beta=\text { artificial compressibility, } \\
\tilde{a}=\text { artificial sound speed }=\beta^{\frac{1}{2}}
\end{array} \right. $$

Where $0.1<\beta<10$.

Thus, the steady incompressible Navier-Stokes equations (two-dimensional Cartesian coordinates) are expressed in a pseudo-transient form as:

$$
\frac{1}{\beta} \frac{\partial p}{\partial \tau}+\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0 $$

$$ \frac{\partial u}{\partial t}+u \frac{\partial u}{\partial x}+v \frac{\partial u}{\partial y}=-\frac{\partial p}{\partial x}+\frac{1}{Re}\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right) $$

$$ \frac{\partial v}{\partial t}+u \frac{\partial v}{\partial x}+v \frac{\partial v}{\partial y}=-\frac{\partial p}{\partial y}+\frac{1}{Re}\left(\frac{\partial^2 v}{\partial x^2}+\frac{\partial^2 v}{\partial y^2}\right)
$$

Above equations are written ion a flux vector form as:

$$
\frac{\partial Q}{\partial t}+\frac{\partial F}{\partial x}+\frac{\partial G}{\partial y}=\frac{1}{{Re}} N \nabla^2 Q
$$

$$
Q=\left(\begin{array}{l}
p \\
u \\
v
\end{array}\right), F=\left(\begin{array}{l}
\beta u \\
u^2+p \\
u v
\end{array}\right), G=\left(\begin{array}{l}
\beta v \\
u v \\
v^2+p
\end{array}\right), N=\left(\begin{array}{lll}
0 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{array}\right)
$$

## Implicit Discretization

$$
\begin{gathered}
\Delta Q+\Delta t\left[\frac{\partial}{\partial x}(A \Delta Q)+\frac{\partial}{\partial y}(B \Delta Q)-\frac{N}{{Re}}\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right) \Delta Q\right] \\
=\Delta t\left[ \frac{\partial E^n}{\partial x}-\frac{\partial F^n}{\partial y}+\frac{N}{\{Re}}\left(\frac{\partial^2 Q}{\partial x^2}+\frac{\partial^2 Q}{\partial y^2}\right)\right]
\end{gathered}
$$

or

![image](https://github.com/user-attachments/assets/30bb62a8-5694-4ef3-a47a-ccc0288d48ec)


For the two-dimensional problems under consideration, the implicit formulation will result in a block penta-diagonal system of equations. As we are aware, the solution of such a system is expensive. To overcome this problem, approximate factorization will be implemented. Second, to overcome any possible instability in the solution, artificial viscosity (damping terms) will be added to the equation. Thus, with the stated considerations, the FDE is formulated as:

$$
\begin{aligned}
& {\left[I+\Delta t\left(\frac{\partial A}{\partial x}-\frac{N}{{Re}} \frac{\partial^2}{\partial x^2}\right)+\varepsilon_i(\Delta x)^2 \frac{\partial^2}{\partial x^2}\right] \Delta Q^*} \\
& \quad=R H S-\varepsilon_e\left[(\Delta x)^4 \frac{\partial^4}{\partial x^4}+(\Delta y)^4 \frac{\partial^4}{\partial y^4}\right] Q
\end{aligned}
$$

A second-order central difference approximation is applied spatially to the above to provide:

![image](https://github.com/user-attachments/assets/fdf54c3b-24b6-4061-b902-7e31cb52c137)

Equation above is now rearranged to provide a block tridiagonal system as follows:

![image](https://github.com/user-attachments/assets/6adec44a-0911-4cba-984c-efc05378d992)

where

![image](https://github.com/user-attachments/assets/107c7569-9c04-452a-acc4-9cbb55dfa32f)

![image](https://github.com/user-attachments/assets/1af85ffd-d3cc-4d68-9eb6-9c5c5c0c9dd7)


## Boundary Conditions
### Left Boundary (i=2)
The left boundary condition controls the equation for $U$. Since the velocity components are invariant at the left boundary, specific conditions for pressure are set using a second-order forward FDE.

### Lower and Upper Walls
Derivations for the lower and upper walls are provided, along with conditions at the outlet boundary (i = imax-1).

## Solution Algorithm

![image](https://github.com/user-attachments/assets/f21e7b79-ff14-4b3b-add5-877006ea2bd5)


## Generated Grids
- **Debugging Grid 51x21**  
  ![image](https://github.com/user-attachments/assets/efa351f1-5fe3-4225-bbd0-ba92f2898169)

- **Final Grid 141x51**  
 ![image](https://github.com/user-attachments/assets/592de4cb-01f0-49b8-89f7-93a408876543)


## Results

### Convergence History
- Convergence history of the solution based on the divergence of the velocity vector.
- ![image](https://github.com/user-attachments/assets/e91e0fe3-5b17-4e2c-a71e-2477716b720b)


### Wall Shear
- Wall shear for both upper and lower walls.
- ![image](https://github.com/user-attachments/assets/a99a72e0-9d1a-4a5b-8ba8-b59e206ef732)

- ![image](https://github.com/user-attachments/assets/5c46194e-91ab-44f0-be5f-74065fb94a39)


### Streamlines and Velocity Vectors
Streamlines and velocity vectors at various Reynold numbers follow the pattern shown in the pressure contours section.

- ** Re = 800 **
- ![image](https://github.com/user-attachments/assets/2fff5a82-18a7-45c1-9ee9-ac526ad05abf)

- ** Re = 600 **
- ![image](https://github.com/user-attachments/assets/a5d00e74-e108-4913-ab6a-46b77570a64b)


- ** Re = 400 **
- ![image](https://github.com/user-attachments/assets/d1d0aa04-79b8-4db8-9621-be7c91b14ce4)


- ** Re = 200 **
- ![image](https://github.com/user-attachments/assets/25ce042a-3d11-4b5f-aab4-b386b4349773)
