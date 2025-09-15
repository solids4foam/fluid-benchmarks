# Flow Past a Horizontally Oscillating Sphere in a Cube Cavity

---

## Tutorial Aims

- Benchmark the accuracy of the `newtonIcoFluid` coupled Newton-Raphson fluid
  solver for flow around a horizontally oscillating sphere in a cube cavity.

---

## Case Overview

Gilmanov and Sotiropoulos[^gilmanov2005] provided a numerical solution for the
 3-D, incompressible, laminar, viscous flow past a rigid sphere that oscillates
 sinusoidally in a cubic cavity.

### Solution Domain

The cubic cavity is defined over the domain $`[-1, 1] \times [-1, 1] \times [-1,
 1] \, \mathrm{m}`$, and a sphere of diameter of $`D = 1 \, m`$ is initially
 located at the origin, $`(0, 0, 0) \, \mathrm{m}`$, (see figure below).

![Domain and dimensions.](./images/domain.png)

### Governing Equations

To account for the mesh motion resulting from the oscillating traverse of the
 `sphere` patch in the cavity, the ALE approach is adopted. In this formulation,
 the continuity and momentum equations are expressed as follows:

```math
\begin{gather}
%
% Continuity equation
%
  \int_{\Omega(t)} \nabla \cdot \mathbf{U} \mathrm{d}V = 0, \\
%
% Momentum equation
%
  \frac{\partial}{\partial t} \int_{\Omega(t)} \mathbf{U} \mathrm{d}V
+ \int_{\Omega(t)} \nabla
  \cdot
  \big[
      \mathbf{U} \otimes (\mathbf{U}
    - \mathbf{U}_{\mathrm{g}})
  \big] \mathrm{d}V
- \int_{\Omega(t)} \nu \nabla^2 \mathbf{U} \mathrm{d}V
- \int_{\Omega(t)} \nabla p \mathrm{d}V
= 0,
\end{gather}
```

where $`\mathbf{U}`$ denotes the fluid velocity, $`p`$ is the kinematic
 pressure, $`\nu`$ is the kinematic viscosity of the fluid, $`\rho`$ represents
 the fluid density, and $`\mathbf{U}_{\mathrm{g}}`$ is the grid (mesh) velocity.
 The grid velocity $`\mathbf{U}_g`$ is obtained by solving the **Space
 Conservation Law (SCL)**, also known as the **Grid Conservation Law (GCL)**:

```math
\begin{equation}
    \frac{\partial}{\partial t} \int_{\Omega(t)} \mathrm{d}V
  - \int_{\Omega(t)} \nabla \cdot \mathbf{U}_{\mathrm g} \mathrm{d}V = 0.
\end{equation}
```

> [!NOTE]
> In the Arbitrary Lagrangian–Eulerian formulation, if $`\mathbf{U}_{\mathrm{g}}
>  = 0`$, we recover the Eulerian description (static mesh). If
> $`\mathbf{U}_{\mathrm{g}} = \mathbf{U}`$, we obtain the Lagrangian
> description (mesh moves with fluid particles). Any other case corresponds to
> the ALE description.

### Sphere (Boundary) Motion

The displacement of points on the `sphere` patch, with an initial sphere center
 at $`(x_0, y_0, z_0) = (0, 0, 0)`$ are given by:

```math
\begin{equation}
  \mathbf{x}_{\text{sph}}
  = \left(
        \mathbf{x}_{\text{sph}_0}^x
      + h (1 - \cos(2 \pi t))
    \right) \mathbf{i}
  + \mathbf{x}_{{\text{sph}}_0}^y \mathbf{j}
  + \mathbf{x}_{{\text{sph}}_0}^z \mathbf{k},
\end{equation}
```

where $`t`$ is time, $`h = D/8 \, \mathrm{m}`$ is the oscillation amplitude, and
 $`\mathbf{x}_{\text{sph}_0}`$ is the initial location points on the `sphere`
 patch, where $`\mathbf{x}_{\text{sph}_0}^i`$, $`(i = x, y ,z)`$, are its $`x`$,
 $`y`$, and $`z`$ components respectively. This motion is specified in the
 `0/pointDisplacement` file.

> [!NOTE]
> If $`h = D/8 = 0.125 \, \mathrm{m}`$, and that $`\cos(2 \pi t) \in [-1, 1]`$,
> then the sphere's centre oscillates over $`[0, 2h] = [0, 0.25] \,
> \mathrm{m}`$. Thus, the centre remains in the positive half of the cavity, in
> the $x$-direction, while the sphere extends from $`x = -0.5 \, \mathrm{m}`$ to
> $`x = 0.75 \, \mathrm{m}`$, remaining at least $`0.25 \, \mathrm{m}`$ away
> from the `xmax` patch and $`0.5 \, \mathrm{m}`$ away from the `xmin` patch.

> [!NOTE]
> Times of extrema:
>
> - Minima at
>
> ```math
> \begin{gather}
>   \cos(2 \pi t) = 1, \\
>   \begin{aligned}
>     &\rightarrow 2 \pi t = 2 \pi k, \\
>     &\rightarrow t = k = 0, 1, 2, \ldots \, \mathrm{s}.
>   \end{aligned}
> \end{gather}
> ```
>
> - Maxima at
>
> ```math
> \begin{gather}
>   \cos(2 \pi t) = -1, \\
>   \begin{aligned}
>     &\rightarrow 2 \pi t = 2 \pi (k + 1/2), \\
>     &\rightarrow t = k + \tfrac{1}{2} = 0.5, 1.5, 2.5, \ldots \, \mathrm{s}.
>   \end{aligned}
> \end{gather}
> ```

### Boundary condition

A no-slip wall condition ($`\mathbf{U} = \mathbf{\bar{U}}`$) and a zero normal
 pressure gradient ($`{\partial p}/{\partial n} = 0`$) are used on all solid
 surfaces. On the outer boundaries, $`\mathbf{\bar{U}} = \mathbf{0}`$, while on
 the `sphere` patch, $`\mathbf{\bar{U}}`$ equals the patch velocity.

The fluid is initially at rest $`\mathbf{U}_0 = \mathbf{0} \, \mathrm{m/s}`$ and
 the pressure in the domain is set to $`p = 0 \, \mathrm{Pa}`$. As all pressure
 boundaries are of gradient type, the pressure field is defined up to a
 constant. To deal with this in the current case, the pressure in one cell of
 the domain is fixed to $`0 \, \mathrm{Pa}`$. For each model, this is specified
 by setting the `pRefPoint` entry in the following dictionaries:

- `newtonIcoFluid`: `constant/fluidProperties.newtonIcoFluidCoeffs`,
- `pimpleFluid`: `system/fvSolution.PIMPLE`.

The `xmin` and `xmax` patches are assigned a zero displacement condition, while
 `ymin`, `ymax`, `zmin` and `zmax` use a _slip_ boundary condition, which
 constrains motion in the normal direction, $`\mathbf{x}_{\mathrm{g}} \cdot
 \mathbf{n} = 0 \, \mathrm{m}`$, while allowing tangential displacement with the
 same displacement as the nearest interior cell, satisfying $`{\partial
 (\mathbf{x}_g \cdot \mathbf{t})}/{\partial n} = 0`$, where $`\mathbf{x}_g`$
 represents the coordinate of a point on the patch.

![Initial and boundary conditions](./images/initial_and_boundary_condition.png)

Also, note that since the `sphere` patch is moving, we need to use the motion
 boundary condition, `newMovingWallVelocity`, in the `0/U` file.

```foam
boundaryField
{
    sphere
    {
        type            newMovingWallVelocity;
        value           uniform (0 0 0);
    }
    ...
}
```

> [!NOTE]
> `newMovingWallVelocity` is a version of `movingWallVelocity`, where boundary
> non-orthogonal corrections are included.

### Constant Properties

The fluid is Newtonian with the density, $`\rho = 1 \mathrm{kg}/\mathrm{m}^3`$,
 (needed by `newtonIcoFluid` fluid model). In this case, the kinematic
 viscosity, `nu`, is calculated so that the Reynolds number is 20.

``` math
\begin{equation}
  \mathrm{Re} = \frac{|\mathbf{U}|_{\mathrm{max}} D}{\nu} = 20,
\end{equation}
```

where

```math
\begin{equation}
    |\mathbf{U}|_{\text{max}}
  = |\mathbf{U}|_{\text{sph}_\text{max}}
  = 2 \pi h = 2 \pi \frac{D}{8} \, \frac{\mathrm{m}}{\mathrm{s}}.
\end{equation}
```

Therefore, with $`D = 1 \, \mathrm{m}`$ and $`\mathrm{Re} = 20`$,

```math
\begin{align}
  \nu =& \frac{2 \pi \frac{D}{8} D}{\mathrm{Re}} \\
      =& \frac{2 \pi}{8 \times 20} = \frac{\pi}{80}
      =  0.03926990817 \, \frac{\mathrm{m}^2}{\mathrm{s}}.
\end{align}
```

These properties are specified in `constant/transportProperties`

```foam
transportModel  Newtonian;

rho             rho [1 -3 0 0 0 0 0] 1;
nu              nu  [0 2 -1 0 0 0 0] 0.03926990817;
```

Since the flow is laminar, `constant/turbulenceProperties` contains:

```foam
simulationType  laminar;
```

Gravity has no effect in this problem; therefore, the `constant/g` file
 contains:

```foam
dimensions      [0 1 -2 0 0 0 0];
value           (0 0 0);
```

In the `constant/physicsProperties` file, we select the fluid analysis as the
 physics type.

```foam
type fluid;
```

Next, we select a fluid model and specify its coefficients in the
`constant/fluidProperties` file.

```foam
fluidModel ${FLUID_MODEL:-pimpleFluid};

pimpleFluidCoeffs
{}

newtonIcoFluidCoeffs
{
    // Pressure stabilisation
    stabilisation
    {
        type        RhieChow;
        scaleFactor 1.0;
    }

    // Set the pressure to zero at the centre cell
    pRefPoint       (0.5 0.2 0);
    pRefValue       0.0;

    // PETSc options file used by PETSc SNES
    optionsFile     petscOptions.pcsplit3D;
}
```

> [!NOTE]
> Coefficients for `pimpleFluid` should be specified in the
> `system/fvSolution.PIMPLE` dictionary.

> [!TIP]
> `${FLUID_MODEL:-newtonIcoFluid}` uses shell parameter expansion, meaning it
> will insert the value of the environment variable `FLUID_MODEL` if it's
> non-empty, or fall back to `newtonIcoFluid` if empty. This allows us to set
> `fluidModel` dynamically without editing the file and by simply setting the
> environment variable before running the case, e.g., `FLUID_MODEL=pimpleFluid
>  solids4Foam` or `FLUID_MODEL=pimpleFluid ./Allrun`

---

## Initial Mesh

Five types of meshes are availabe in the case:

- Structured hexagonal uniform mesh (`MESH=HEX`),
- Structured hexagonal graded mesh (`MESH=HEX_GRADED`),
- Unstructured hexagonal mesh (`MESH=CARTESIAN`),
- Unstructured tetrahedral mesh (`MESH=TRI`),
- Unstructured polyhedral mesh (`MESH=POLY`).

Examples of the above mesh types are shown in the following figure:

![Various mesh types: (a) `MESH=HEX`, (b) `MESH=CARTESIAN`, (c) `MESH=TET`, (d)
 `MESH=POLY`.](./images/meshes.png)

---

## Mesh Motion

The dynamic mesh details are specified in `constant/dynamicMeshDict`, where the
 type of dynamic mesh is set to `dynamicMotionSolverFvMesh`, which uses a mesh
 motion solver to compute the mesh motion based on a prescribed motion.

```foam
dynamicFvMesh       dynamicMotionSolverFvMesh;
motionSolverLibs    (fvMotionSolvers);
```

In this case, among other mesh motion solvers implemented in OpenFOAM, we use
 `displacementLaplacian`, which is a diffusion-based solver that smoothly
 propagates motion throughout the domain.

```foam
motionSolver        displacementLaplacian;
```

It solves a Laplace equation for cell displacement field, `cellDisplacement`,
 which is stored at cell centers,

```math
\begin{equation}
    \nabla \cdot
    \left(
      \Gamma \nabla (\Delta \mathbf{X}_{\text{cell}})
    \right)
  = 0,
\end{equation}
```

where $`\Delta \mathbf{X}_{\text{cell}}`$ is the cell displacement field and
 $`\Gamma`$ is the motion _diffusivity_ (scalar) field.

The diffusivity field determines how the motion of the boundary patch is
 distributed throughout the domain. The value for diffusivity is determined by
 specifying a diffusivity model via the
 `displacementLaplacianCoeffs.diffusivity` dictionary entry:

```foam
displacementLaplacianCoeffs
{
    diffusivity     inverseDistance (sphere);
}
```

In this case we use the `inverseDistance` model which calculates diffusivity
 based on inverse distance to given patches, which in this case is the `sphere`
 patch. This means that mesh points farther away from the sphere patch are
 deformed more, while the cells near to the sphere act stiffer and deform
 less, thus preserving their quality.

As we have selected the `displacementLaplacian` motion solver, the boundary
 conditions for the mesh displacement field should be prescribed in
 `0/pointDisplacement`. The `sphere` displacement is specified as a rigid body
 motion using the `solidBodyMotionDisplacement` boundary condition:

```foam
boundaryField
{
    sphere
    {
        type                    solidBodyMotionDisplacement;
        solidBodyMotionFunction oscillatingLinearMotion;
        /*
            Reference displacement function written in terms of cosine
            x(t) = -h*\cos(2*\pi*t) + h
                 = -h*\sin(2*\pi*t + \pi/2) + h
                 = -h*\sin(2*\pi*(t + (\pi/2)/(2\pi))) + h
                 = -h*\sin(2*\pi*(t + 1/4)) + h

            oscillatingLinerMotion function:
            x =  A*\sin(B*(t + C)) + D
        */
        oscillatingLinearMotionCoeffs
        {
            h               0.125;          // [m]

            amplitude       (-$h 0 0);      // [m] vector, x-direction motion
            omega           ${{ 2*pi() }};  // [rad/s] (period = 1s)
            phaseShift      0.25;           // [s]
            verticalShift   ($h 0 0);       // [m] vector shift
        }
    }
    ...
}
```

Note that the cosine function $`x(t) = h \left(1 - \cos(2\pi t)\right)`$ was
converted to its equivalent sine form, $`-h \sin\left(2\pi(t + 1/4)\right) +
h`$, as follows:

```math
\begin{align}
  x(t) &= -h \cos(2\pi t) + h \\
       &= -h \sin(2\pi t + \pi/2) + h \\
       &= -h \sin\left(2\pi\left(t + \tfrac{\pi/2}{2\pi}\right)\right) + h \\
       &= -h \sin\left(2\pi(t + 1/4)\right) + h \, \mathrm{m}, \\
\end{align}
```

Then a pattern matching was performed against the `oscillatingLinearMotion`
function, with an optional _phase-_ and a _vertical-shift_.

```math
\begin{equation}
  x = A \sin(B (t + C)) + D,
\end{equation}
```

where

- $`x`$: Displacement $`[\mathrm{m}]`$ (vector)
- $`A = (-h, 0, 0)`$: Amplitude $`[\mathrm{m}]`$ (vector)
- $`B = 2\pi`$: Radial velocity $`[\mathrm{rad/s}]`$ (scalar)
- $`C = 1/4`$: Phase-shift to left $`[\mathrm{s}]`$ (scalar)
- $`D = (h, 0, 0)`$: Vertical shift $`[\mathrm{m}]`$ (vector)

> [!NOTE]
> `h` is a user-defined dictionary entry here; It is not a mandatory entry
> required by the `oscillatingLinear` function.

For the `xmin` and `xmax` patches, a no-slip boundary condition is applied,
 while the `ymin`, `ymax`, `zmin` and `zmax` patches use a slip boundary
 condition.

```foam
boundaryField
{
    ...
    "xmin|xmax"
    {
        type    fixedDisplacement;
        value   (0 0 0);
    }
    "ymin|ymax|zmin|zmax"
    {
        type    slip;
    }
    ...
}
```

---

## Running the Case

The case can be run using the included `Allrun` script located in the template
 case (`oscillatingSphereInCube/base`) directory:

```bash
# Serial run with default values of environment variables
# Please refer to the top of `./Allrun` to inspect or change the default values
./Allrun

# Specify the mesh density level
# options: $START_MESH <= MESH_LEVEL <= $END_MESH
MESH_LEVEL=4 ./Allrun

# Specify the mesh type
# HEX, HEX_GRADED, CARTESIAN, TET, POLY
MESH_LEVEL=8 MESH=POLY ./Allrun

# Specify the physics model to use
# options: newtonIcoFluid, pimpleFluid
FLUID_MODEL=pimpleFluid ./Allrun

# Parallel run (currently using `scotch` method)
NUMBER_OF_SUBDOMAINS=8 ./Allrun -p

# Example:
MESH_LEVEL=8 MESH=POLY \
NUMBER_OF_SUBDOMAINS=32 \
FLUID_MODEL=pimpleFluid \
    ./Allrun -p
```

or by executing the `Allrun` script in the parent `oscillatingSphereInCube`
 directory.

```bash
# Serial run
./Allrun

# Parallel run (currently using `scotch` method)
PARALLEL=8 ./Allrun
PARALLEL=8 NUMBER_OF_SUBDOMAINS=32 ./Allrun
```

The `oscillatingSphereInCube/Allrun` script runs multiple cases based on the
 configurations and a range of mesh density level given at the top of the
 script, e.g.

```bash
# Define configurations as space-separated strings
configs=(
    "BASE=base NAME=cartesian.pimpleFluid MESH=CARTESIAN PIMPLEFLUID=1"
)
readonly configs
# Description
#     BASE:            name of the base template case
#     NAME:            base name given to each case that is created
#     MESH:            type of the initial mesh
#                      - structured hexahedral mesh (HEX or HEX_GRADED)
#                      - unstructured hexagonal mesh (CARTESIAN)
#                      - unstructured tetrahedral mesh (TET)
#                      - unstructured polyhedral mesh (POLY)
#     PIMPLEFLUID:     use the pimpleFluid (1) or the newtonIcoFluid (0) models
#
# Example
# configs=(
#     "BASE=base NAME=poly.newtonIcoFluid MESH=POLY PIMPLEFLUID=0"
#     "BASE=base NAME=poly.pimpleFluid MESH=POLY PIMPLEFLUID=1"
# )

# Mesh level (density) range (integer).
declare -i START_MESH=1
declare -i END_MESH=1

# ...
```

The `Allrun` script automates the setup, execution, and post-processing of the
 cases. It creates a unique run directory, iterates over predefined case
 configurations and mesh refinement levels, generates the mesh using one of
 `blockMesh`, cfMesh, or Gmsh, and selects the fluid model (`newtonIcoFluid` or
 `pimpleFluid`). Simulations are then run in serial (`PARALLEL=0`) or parallel
 (`PARALLEL=1`), while performance metrics are logged using some version of
 `time` utility. Finally, plots of the drag coefficient, $`C_d`$, is generated
 using Gnuplot.

---

## Expected Results

The 2-D contour plots ($x$-$y$) at four different time instances of $`t = 0 \,
 s`$, $`T / 4 \, s`$, $`T / 2 \, s`$ and $`3 T / 4 \, s`$, ($`T = 1 \,
 \mathrm{s}`$ being the period of the oscillation), are illustrated below for
 $`Re = 20`$.. The results show that the $x$-component of the velocity reaches
 higher magnitudes when the sphere passes $`x = 0.125 \, \mathrm{m}`$).

<!-- TODO(abzrg): a note about the position of the sphere and streamlines at
t=T/4s and t=3T/4s -->

<!-- ![TODO(abzrg): ADD FIGURE: CONTOUR PLOT AND STREAM TRACES](placeholder) -->

To quantitatively compare the results, in the following figure, the drag
 coefficient, $C_{\mathrm d}$, is plotted along with the results reported by
 Erzincanli and Sahin [^erzincanli2013]. The comparison shows that the present
 numerical results are visually indistinguishable from those of the reference.

<!-- ![TODO(abzrg): ADD FIGURE: DRAG COEFFICIENT PLOT](placeholder) -->

<!-- ![TODO(abzrg): ADD FIGURE: LIFT COEFFICIENT PLOT](placeholder) -->

### Function Objects for Calculating Forces and Force Coefficients

The $`x`$-axis is defined as the drag direction, with the drag coefficient given
 by

```math
\begin{equation}
    C_{\mathrm{d}}
  = \frac{2 F_x}{\rho_\infty \mathbf{U}_\infty^2 A_{\text{ref}}},
\end{equation}
```

Here, $F_x$ is the $`x`$-components of the fluid force acting on the sphere.
 $`\rho`$ is the freestream fluid density ($`1 \mathrm{kg/m^3}`$), $`U_\infty =
 2 \pi h`$ is the magnitude of the freestream flow velocity relative to the
 sphere, which in this case is considered to be the maximum velocity of the
 sphere (while passing through the channel center). $`A_{\text{ref}} = \pi D^2 /
 4`$ is the reference area, which in this case is the area of the sphere
 projected to the plane perpendicular to flow in the $`x`$-direction.

Performing a simple calculation,

``` math
\begin{align}
  \mathbf{U}_\infty &= \mathbf{U}_{\text{max}}^{x} = 2 \pi h \\
                   &= 2 \times \pi \times 0.125 \\
                   &\approx 0.7854 \, \frac{\mathrm{m}}{\mathrm{s}},
\end{align}
```

``` math
\begin{align}
  p_{\text{dyn}} &= \frac{1}{2} \rho_{\infty} U_\infty^2 \\
               &= 0.5 \times 1.0 \times (0.6169)^2 \\
               &\approx 0.3084 \, \mathrm{Pa},
\end{align}
```

```math
\begin{align}
  A_{\text{ref}} &= \frac{\pi D^2}{4}
               &= \frac{\pi \times 1^2}{4}
               &\approx 0.7854 \, \mathrm{m}^2,
\end{align}
```

```math
\begin{align}
  \texttt{forceScaling} &= \frac{1}{A_{\text{ref}} p_{\text{dyn}}}
                        &= \frac{1}{0.7854 \times 0.3084}
                        &\approx 4.1285 \, \mathrm{N}^{-1},
\end{align}
```

thus, the total fluid forces (in the $`x`$-direction) acting on the sphere are
 scaled by a factor of $`4.1285 \, \mathrm{N}^{-1}`$, to ultimately yield the
 drag coefficient.

```math
\begin{align}
  C_{\mathrm{d}} = \texttt{forceScaling}\, F_x,
\end{align}
```

The `forceCoeffs` function object is employed to, among other coefficients,
  computes the drag coefficient.

```foam
functions
{
    // Generate drag coefficient, Cd, over the sphere patch
    forceCoeffs
    {
        type            forceCoeffs;
        libs            (forces);

        // Fluid force and moment on the following patches will be calculated
        patches         (sphere);

        // Only report the drag coefficient, Cd, and omit others
        coefficients    (Cd);

        // Diameter of the sphere
        D               1.0;

        // Amplitude of the displacement function
        h               ${{ 0.125*$D }};

        // Max velocity is the x direction
        magUxMax        ${{ 2.0 * pi() * $h }};

        // Freestream velocity magnitude
        magUInf         $magUxMax;

        // Reference length
        lRef            $D;

        // Reference area
        Aref            ${{ 0.25 * pi() * $D * $D }};

        // Center of rotation
        CofR            (0 0 0);

        // Direction of the drag force
        // C_d = F_x/(0.5*ρ*U^2*A)
        dragDir         (1 0 0);

        // Direction of the lift force
        // C_l = F_y/(0.5*ρ*U^2*A)
        liftDir         (0 1 0);

        // Field names
        p               p;
        U               U;
        rho             rhoInf; // Indicates incompressible
        rhoInf          1;      // Redundant for incompressible-flow cases

        writeControl    timeStep;
        writeInterval   1;
        log             off;
    }
}
```

> [!NOTE]
> Since the motion of the sphere in the cubic cavity creates a symmetric flow,
> other components of the force coefficients are zero.

> [!TIP]
> The chosen solver computes these quantities on the fly; however, if the
> outputs are missing, they can be re-extracted after the simulation using:
> `pimpleFoam -postProcess`.

> [!TIP]
> The `forceCoeffs` function stores the resulting
> computations in `${FOAM_CASE}/postProcessing/forceCoeffs/0/coefficient.dat`.

---

## References

[^gilmanov2005]: [A. Gilmanov, F. Sotiropoulos, A hybrid Cartesian/immersed
 boundary method for simulating flows with 3D, geometrically complex, moving
 bodies,J. Comput. Phys. 207 (2005)
 45–492.](https://doi.org/10.1016/j.jcp.2005.01.020)

[^erzincanli2013]: [Erzincanli, B. and Sahin, M. (2013). An arbitrary Lagrangian
 –Eulerian formulation for solving moving boundary problems with large
 displacement and rotations. Journal of Computational Physics 255,
 660–679.](https://doi:10.1016/j.jcp.2013.08.038)
