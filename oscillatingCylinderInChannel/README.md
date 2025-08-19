# The Flow Past an Oscillating Circular Cylinder in a Channel

---

## Tutorial Aims

- Benchmarks the accuracy of the `newtonIcoFluid` coupled Newton-Raphson fluid
  solver against the numerical solution of the flow past an oscillating cylinder
  case.

---

## Case Overview

Wan \& Turek[^wan2006] provided a numerical solution for the 2-D,
 incompressible, laminar, viscous flow past a circular cylinder that oscillates
 sinusoidally in an enclosed channel. The case is well-suited for benchmarking
 the accuracy of a discretization.

### Solution Domain

The rectangular channel has a domain $[0, 2.2] \times [0, 0.41] \mathrm{m}^2$,
 and the cylinder is initially located at $(1.1, 0.2) \mathrm{m}$ with respect
 to the lower left corner of the domain.

![Domain and dimensions](./images/domain_dimensions.png)

### Governing Equations

When the cylinder patch (`cylinder`) begins oscillating, the surrounding
 mesh deforms, and the volume of individual cells change throughout the
 simulation. Consequently, we can no longer assume $\mathrm{d}V / \mathrm{d}t =
 0$ in the governing equations. The continuity and momentum equations therefore
 take the following forms:

$$
\int_{\Omega(t)} \nabla \cdot \mathbf{U} \mathrm{d}V = 0,
$$

$$
  \frac{\partial}{\partial t} \int_{\Omega(t)} \mathbf{U} \mathrm{d}V
\+ \int_{\Omega(t)} \nabla \cdot
  \big[\mathbf{U} \otimes (\mathbf{U} - \mathbf{U}_{\mathrm g})\big] \mathrm{d}V
\- \int_{\Omega(t)} \nu \nabla^2 \mathbf{U} \mathrm{d}V
\- \int_{\Omega(t)} \nabla p \mathrm{d}V
\= 0,
$$

where $\mathbf{U}$ denotes the fluid velocity, $p$ is the kinematic pressure,
 $\nu$ is the kinematic viscosity of the fluid, $\rho$ represents the fluid
 density, and $\mathbf{U}_{\mathrm g}$ is the grid (mesh) velocity. The grid
 velocity $\mathbf{U}_g$ is obtained by solving the **Space Conservation Law
 (SCL)**, also known as the **Grid Conservation Law (GCL)**:

$$
\frac{\partial}{\partial t} \int_{\Omega(t)} \mathrm{d}V
\- \int_{\Omega(t)} \nabla \cdot \mathbf{U}_{\mathrm g} \mathrm{d}V = 0.
$$

```note
The above formulation is known as the Arbitrary Lagrangian–Eulerian (ALE)
 formulation. If $\mathbf{U}_{\mathrm{g}} = 0$, we recover the Eulerian
 description (static mesh). If $\mathbf{U}_{\mathrm{g}} = \mathbf{U}$, we obtain
 the Lagrangian description (mesh moves with fluid particles). Any other case
 corresponds to the ALE description.
```

### Cylinder (Boundary) Motion

The displacement of points on the `cylinder` patch, with an initial cylinder
 center at $(x_0, y_0) = (1.1, 0.2)$ (relative to the lower-left corner of the
 domain), is given by:

$$
\mathbf{x}_{\mathrm cyl}
\= \left(x_{\mathrm{cyl}_0} + A \sin(\omega t)\right) \mathbf{i}
\+ y_{{\mathrm cyl}_0} \mathbf{j},
$$

where $t$ is time, $A = 0.25 \mathrm m$ is the oscillation amplitude, $\omega =
 2\pi f$ is the angular velocity, $f = 0.25 \mathrm Hz$ is the oscillation
 frequency, and $\mathbf{x}_{{\mathrm cyl}_0}$ is the initial location points on
 the `cylinder` patch. This motion is specified in the `0/pointsDisplacement`
 file.

### Boundary condition

A no-slip boundary condition, $\mathbf{U} = \mathbf{0} \mathrm{m/s}$, and a
 zero-pressure-gradient condition, ${\partial p}/{\partial n} = 0$, are applied
 on all solid surfaces. The fluid is initially at rest $\mathbf{U}_0 =
 \mathbf{0} \mathrm{m/s}$. As all pressure boundaries are of gradient type, the
 pressure field is defined up to a constant. To deal with this in the current
 case, the pressure in one cell of the domain, close to the point specified by
 the dictionary entry `system/fvSolution.PIMPLE.pRefPoint`, is fixed to $0
 \mathrm{Pa}$.

The `left` and `right` patches are assigned a fixed-value zero point
 displacement, i.e., _no-slip_, $\mathbf{x}_{\mathrm{g}} = \mathbf{0}
 \mathrm{m}$ while the `up` and `down` patches use a _slip_ boundary condition,
 which constrains motion in the normal direction, $\mathbf{x}_{\mathrm{g}} \cdot
 \mathbf{n} = 0 \mathrm{m}$, while allowing tangential displacement with the
 same displacement as the nearest interior cell, satisfying ${\partial
 (\mathbf{x}_g \cdot \mathbf{t})}/{\partial n} = 0$.

![Initial and boundary conditions](./images/initial_and_boundary_condition.png)

The slip boundary condition preserves the initial orthogonality and minimizes
 mesh skewness during the cylinder's horizontal oscillations: the `up` and
 `down` boundaries allow tangential sliding, while the mesh to the left and
 right of the cylinder undergoes compression and expansion like an accordion,
 preventing excessive distortion and maintaining quality throughout the domain.
 For example, the following figure shows a structured hexagonal mesh
 (`MESH_LEVEL=2`) at time $t = 3.00 \mathrm{s}$, which corresponds to the
 farthest cylinder can move toward left, as $\sin(2 \pi f t) = -1$ .

$$
  x_{\text{cell}}\big|_{t=3.0\mathrm s}
= A \sin(2 \pi f t)
= 0.25 \times \sin(2 \times 0.25 \times \pi \times 3 )
= -0.25 \mathrm{m}
$$

At this instance of time, and similarly at $t = 9.00 \mathrm{s}$, the mesh has
 the highest deformation.

![Slip boundary condition](./images/mesh_t_3.00s.png)

Also, note that since the `cylinder` patch is moving, we need to use the motion
 boundary condition, `movingWallVelocity`, in the `0/U` file.

```{.cpp .numberLines}
boundaryField
{
    cylinder
    {
        type            movingWallVelocity;
        value           uniform (0 0 0);
    }
...
```

### Constant Properties

The fluid is Newtonian and its kinematic viscosity, $\nu = \mu / \rho = 1 \times
 10^{-3} \mathrm{m}^2/\mathrm{s}$, and density, $\rho = 1
 \mathrm{kg}/\mathrm{m}^3$, (needed by `newtonIcoFluid` fluid model) is
 specified in `constant/transportProperties`

```foam
transportModel  Newtonian;

rho             rho [1 -3 0 0 0 0 0] 1; // Required by newtonIcoFluid
nu              nu  [0 2 -1 0 0 0 0] 1e-3;
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
fluidModel ${FLUID_MODEL:-newtonIcoFluid}; // or pimpleFluid

pimpleFluidCoeffs
{}

newtonIcoFluidCoeffs
{
    // Trottenberg gives 1/16 for Stokes flow
    omega           [0 -2 1 0 0 0 0] 0.625;
    localReRef      1;
    omegaExponent   1.0;
    alphaU          1.0;

    // Set the pressure to zero at the centre cell
    pRefPoint       (0.5 0.2 0);
    pRefValue       0.0;

    // PETSc options file used by PETSc SNES
    optionsFile     petscOptions.lu;
}
```

```note
Coefficients for `pimpleFluid` can also be specified in the
 `system/fvSolution.PIMPLE` dictionary.
```

```tip
`${FLUID_MODEL:-newtonIcoFluid}` uses shell parameter expansion, meaning it will
 insert the value of the environment variable `FLUID_MODEL` if it's set and
 non-empty, or fall back to `newtonIcoFluid` if its unset or empty. This allows
 us to set `fluidModel` dynamically without editing the file and by simply
 setting the environment variable before running the case, e.g.,
 `FLUID_MODEL=pimpleFluid solids4Foam` or `FLUID_MODEL=pimpleFluid ./Allrun`
```

---

## Initial Mesh

Four types of meshes can be generated:

- structured quadrilateral uniform mesh (`MESH=QUAD`),
- structured quadrilateral graded mesh (`MESH=QUAD_GRADED`),
- unstructured triangular mesh (`MESH=TRI`),
- unstructured polyhedral mesh (`MESH=POLY`).

As an example, uniform and graded structured hexagonal grids at three different
 mesh density levels are shown in the following figures:

![Various mesh density levels of the uniform structured mesh: (a) `MESH_LEVEL=1`,
 (b) `MESH_LEVEL=2`, (c) `MESH_LEVEL=3`.](./images/structured_uniform_hex_meshes.png)

An example of unstructured polyhedral mesh generated using Gmsh, and other
 OpenFOAM mesh manipulation utilities, is shown in the following figure:

![An example of unstructured graded polyhedral mesh.](./images/unstructured_graded_poly_mesh.png)

---

## Mesh Motion

Dynamics of mesh is specified in `constant/dynamicMeshDict`, where the type of
 dynamic mesh is set to `dynamicMotionSolverFvMesh`, which uses a mesh motion
 solver to compute the mesh motion based on a prescribed motion, which in this
 case is the inline sinusoidal oscillation of the `cylinder` boundary patch.

```foam
dynamicFvMesh       dynamicMotionSolverFvMesh;
motionSolverLibs    (fvMotionSolvers);
```

In this case, among other mesh motion solvers implemented in OpenFOAM, we use
 `displacementLaplacian`, a diffusion-based solver that smoothly propagates
 motion throughout the domain.

```foam
motionSolver        displacementLaplacian;
```

It solves the following Laplace equation for cell displacement
 (`cellDisplacement`) field (stored at cell centers),

$$
\nabla \cdot \left( \Gamma \nabla (\Delta \mathbf{X}_{\text{cell}}) \right) = 0,
$$

where $\Delta \mathbf{X}_{\text{cell}}$ is the cell displacement field
 (`cellDisplacement`) and $\Gamma$ is the motion diffusivity (scalar) field.

`displacementLaplacian` interpolates the cell displacement at cell centers to
 mesh points to obtain the grid point displacement field, $\Delta
 \mathbf{X}_{\mathrm g}$, and computes new positions based on the initial
 points, $\mathbf{X}_{t_0}$. It then moves the grid points to these new
 locations. This solver smoothly propagates prescribed motion at boundaries or
 within mesh zones to the rest of the domain, preventing skewed or tangled
 cells.

$\Gamma$ determines how the motion of the boundary patch is distributed
 throughout the domain. The value for diffusivity is determined by specifying a
 diffusivity model via the `displacementLaplacianCoeffs.diffusivity` dictionary
 entry.

```foam
displacementLaplacianCoeffs
{
    diffusivity     inverseDistance (cylinder);
}
```

In this case we use the `inverseDistance` model which calculates diffusivity
 based on inverse distance to given patches, which in this case is the
 `cylinder` patch. This means that mesh points farther away from the cylinder
 patch are deformed less.

Finally, as we have selected the `displacementLaplacian` motion solver, the
 prescribed motion of the boundary patch `cylinder` must be set up in the
 `0/pointDisplacement` file.

Since the oscillation of cylinder is a rigid body motion,
 `solidBodyMotionDisplacement` is used as the motion boundary type. This
 specifies all points move together without deformation. This condition allows
 the use of a predefined _solid body motion function_ to set fixed
 displacements.

In this case, the `oscillatingLinearMotion` function is applied, causing the
 patch to oscillate back and forth along a straight line with a given
 _amplitude_ (`amplitude`) and _angular velocity_ (`omega`), as defined in the
 `oscillatingLinearMotionCoeffs` dictionary:

$$
\mathbf{X}_{t_n} = \mathbf{X}_{t_0} + \Delta \mathbf{X}_{\mathrm g}(t)
$$

$$
  \Delta \mathbf{X}_{\mathrm g}(t)
= \texttt{amplitude} \times \sin(\texttt{omega} \times t)
$$

```foam
boundaryField
{
    cylinder
    {
        type                    solidBodyMotionDisplacement;
        solidBodyMotionFunction oscillatingLinearMotion;

        // \Delta X = amplitude*sin(omega*t)
        oscillatingLinearMotionCoeffs
        {
            // Amplitude (x-direction) of displacement oscillations [m]
            Ax              0.25;
            // Frequency of displacement oscillations [Hz]
            f               0.25;

            // Amplitude vector [m]
            amplitude       ($Ax 0 0);

            // Angular velocity of oscillations [rad/s]
            omega           ${{ 2*pi()*$f }};
        }
    }
    ...
}
```

```note
`Ax` and `f` are user-defined variables here; they are not OpenFOAM keywords but
 are used to calculate `amplitude` and `omega`.
```

Alternatively, the function can be implemented manually using a
 `codedFixedValue` boundary condition:

```foam
boundaryField
{
    cylinder
    {
        name    pointDisplacement_cylinder;
        type    codedFixedValue;
        code
        #{
            using constant::mathematical::twoPi;

            constexpr scalar Ax = 0.25;     // Amplitude [m]
            constexpr scalar f  = 0.25;     // Frequency [Hz]
            constexpr scalar omega = twoPi*f;

            const scalar t = this->db().time().value();
            operator==(Ax*sin(omega*t)*vector(1, 0, 0));
        #};
        value   $internalField;
    }
    ...
```

For the `left` and `right` patches, a no-slip boundary condition is applied,
 while the `up` and `down` patches use a slip boundary condition.

```foam
boundaryField
{
    ...
    slipWalls
    {
        type    slip;
    }

    noSlipWalls
    {
        type    fixedValue;
        value   $internalField;
    }
...
```

---

## Numerics

A linear solver must be specified for `cellDisplacement` in `system/fvSolution`,

```foam
solvers
{
    ...
    "cellDisplacement|cellDisplacementFinal"
    {
        solver      GAMG;
        tolerance   1e-6;
        relTol      0;
        smoother    GaussSeidel;
    }
}
```

Also, the Laplacian term in the mesh motion Laplace equation,
 `laplacian(diffusivity, cellDisplacement)`, requires a discretization scheme:

```foam
// system/fvSchemes
...
laplacianSchemes
{
    ...
    laplacian(diffusivity,cellDisplacement) Gauss linear corrected;
}
```

The implicit first order `backward` scheme is used for time derivatives,

```foam
// system/fvSchemes
...
ddtSchemes
{
    default backward;
}
```

---

## Running the Case

The case can be run using the included `Allrun` script located in the template
 case (`oscillatingCylinderInChannel/base`) directory:

```bash
# Serial run with default values of environment variables
# Please refer to the top of `./Allrun` to inspect or change the default values
./Allrun

# Specify the mesh density level
# options: $START_MESH <= MESH_LEVEL <= $END_MESH
MESH_LEVEL=4 ./Allrun

# Specify the physics model to use
# options: newtonIcoFluid, pimpleFluid
FLUID_MODEL=pimpleFluid ./Allrun

# Parallel run (currently using `scotch` method)
NUMBER_OF_SUBDOMAINS=8 ./Allrun -p

# Example:
MESH_LEVEL=8 NUMBER_OF_SUBDOMAINS=32 FLUID_MODEL=pimpleFluid ./Allrun -p
```

or by executing the `Allrun` script in the parent `oscillatingCylinderInChannel`
 directory.

```bash
# Serial run
./Allrun

# Parallel run (currently using `scotch` method)
PARALLEL=8 ./Allrun
PARALLEL=8 NUMBER_OF_SUBDOMAINS=32 ./Allrun
```

The `oscillatingCylinderInChannel/Allrun` script runs multiple cases based on the
 configurations and a range of mesh density level given at the top of the
 script, e.g.

```bash
# Define configurations as space-separated strings
configs=(
    "BASE=base NAME=hex.uniform.lu.pimpleFluid GRADED=0 PIMPLEFLUID=1"
)
# Example
# configs=(
#     "BASE=base NAME=hex.uniform.lu.newtonIcoFluid GRADED=0 PIMPLEFLUID=0"
#     "BASE=base NAME=hex.uniform.lu.pimpleFluid    GRADED=0 PIMPLEFLUID=1"
#     "BASE=base NAME=hex.graded.lu.newtonIcoFluid  GRADED=1 PIMPLEFLUID=0"
#     "BASE=base NAME=hex.graded.lu.pimpleFluid     GRADED=1 PIMPLEFLUID=1"
# )

# Mesh level (density) range (integer).
declare -i START_MESH=1
declare -i END_MESH=1
```

The `Allrun` script automates the setup, execution, and post-processing of the
 cases. It creates a unique run directory, iterates over predefined case
 configurations and mesh refinement levels, generates the mesh using the
 appropriate `blockMeshDict` (uniform or graded), and selects the fluid model
 (`newtonIcoFluid` or `pimpleFluid`). Simulations are then run in serial
 (`PARALLEL=0`) or parallel (`PARALLEL=1`), while performance metrics are logged
 using some version of `time` utility. Finally, plots of the force coefficients
 are generated using Gnuplot.

## Expected Results

The contour plots at four different time instances, together with the reference
 calculations of Wan and Turek[^wan2006] and Erzincanli et al.[^erzincanli2013],
 show that the $x$-component of the velocity reaches higher magnitudes when the
 cylinder passes through the channel center ($x = 1.1 \mathrm m$). Furthermore,
 vortices are formed near the channel corners as well as along the upper and
 lower walls.

<!-- ![TODO(abzrg): ADD FIGURE: CONTOUR PLOT AND STREAM TRACES](placeholder) -->

To compre quantitatively the drag coefficient, $C_{\mathrm d}$, and the lift
 coefficient, $C_{\mathrm l}$, are evaluated against the results reported by Wan
 and Turek[^wan2006]. The comparison shows that the present numerical results
 are virtually identical to those in the reference study.

<!-- ![TODO(abzrg): ADD FIGURE: DRAG COEFFICIENT PLOT](placeholder) -->

<!-- ![TODO(abzrg): ADD FIGURE: LIFT COEFFICIENT PLOT](placeholder) -->

### Function Objects for Calculating Forces and Force Coefficients

In this 2-D problem, the $x$-axis is defined as the drag direction, with the
 drag coefficient given by

$$
  C_{\mathrm d}
= \frac{2 F_x}{\rho_\infty \mathbf{U}_\infty^2 L_{\text{ref}} L_z},
$$

and the $y$-axis as the lift direction, with the lift coefficient defined as

$$
  C_{\mathrm l}
= \frac{2 F_y}{\rho_\infty \mathbf{U}_\infty^2 L_{\text{ref}} L_z}.
$$

Here, $F_x$ and $F_y$ are the components of the fluid force acting on the
 cylinder in the $x$- and $y$-directions, respectively. $\rho$ is the freestream
 fluid density ($1 \mathrm{kg/m^3}$), $\left|\mathbf{U}_\infty\right| = 2 \pi A
 f$ is the magnitude of the freestream flow velocity relative to the cylinder,
 which in this case is considered to be the maximum velocity of the cylinder
 (while passing through the channel center). $L_{\text{ref}}$ is the reference
 length, which in this case is the cylinder diameter, and $L_z$ is the span of
 the computational mesh in the $z$-direction. Also, $L_{\text{ref}} L_z$ is the
 reference area $A_{\text{ref}}$.

Performing a simple calculation,

$$
|\mathbf{U}_\infty| = 2 \pi f A_{\text{ref}}
                     = 2 \times \pi \times 0.25 \times 0.25
               \approx 0.3927\ \frac{\mathrm{m}}{\mathrm{s}},
$$

$$
p_{\text{dyn}} = \frac{1}{2} \rho_{\infty} \mathbf{U}_\infty^2
                = 0.5 \times 1.0 \times (0.3927)^2
          \approx 0.0771\ \mathrm{Pa},
$$

$$
A_{\text{ref}} = L_{\text{ref}} L_z = 0.1 \times 0.1 = 0.01 \mathrm{m}^2,
$$

$$
\texttt{forceScaling} = \frac{1}{A_{\text{ref}} p_{\text{dyn}}}
                      = \frac{1}{0.01 \times 0.0771}
                \approx 1297\ \mathrm{N}^{-1},
$$

the total fluid forces (in the $x$- and $y$-directions) acting
 on the cylinder are scaled by a factor with units of $\mathrm{N}^{-1}$, in this
 case, $1927\ \mathrm{N}^{-1}$, to ultimately yield the force coefficients.

$$
C_{\mathrm{d}} = \texttt{forceScaling}\ F_x,
$$

$$
C_{\mathrm{l}} = \texttt{forceScaling}\ F_y.
$$

The following function objects are employed:

- `forces`: computes the three components of the fluid force (and moments) on
  specified patches (here, the `cylinder` patch).
- `forceCoeffs`: computes drag and lift coefficients, as well as other force
  coefficients.

```foam
functions
{
    forces
    {
        type            forces;
        libs            (forces);

        // Fluid force and moment on the following patches will be calculated
        patches         (cylinder);

        // Field names
        p               p;
        U               U;
        rho             rhoInf;
        // Freestream density
        rhoInf          1;

        // Center of rotation (required for moment calculation)
        CofR            (0 0 0);

        writeControl    $/writeControl;
        writeInterval   $/writeInterval;
        log             on;
    }

    // Generate drag (Cd) and lift (Cl) coefficients over the cylinder patch
    forceCoeffs
    {
        type            forceCoeffs;
        libs            (forces);

        // Fluid force and moment on the following patches will be calculated
        patches         (cylinder);

        // Only report drag (Cd) and lift (Cl) coefficients and omit others
        coefficients    (Cd Cl);

        // Frequency and amplitude of cylinder displacement oscillation
        f               0.25;
        A               0.25;
        // NOTE(abzrg): Here I considered the maximum x-velocity of the cylinder
        //              patch: u_max = max(2*\pi*f*A*cos(2*\pi*f*t))
        magUxMax        ${{ 2*pi()*$f*$A }}; // Maximum x velocity magnitude

        // Freestream velocity magnitude
        magUInf         $magUxMax; // 0.39269908169872414

        // Span of mesh in the z (depth) direction
        // See blockMeshDict.vertices
        meshSpanZ       0.1;
        // Diameter of the cylinder
        cylDiameter     0.1;
        // Reference length
        lRef            $cylDiameter;
        // Reference area (projected area per unit depth in 2D)
        Aref            ${{ $meshSpanZ * $cylDiameter }};

        // Center of rotation
        CofR            (0 0 0);

        // Direction of the drag force
        // C_d = F_x/(0.5*\rho*U^2*A)
        dragDir         (1 0 0);

        // Direction of the lift force
        // C_l = F_y/(0.5*\rho*U^2*A)
        liftDir         (0 1 0);

        // Field names
        p               p;
        U               U;
        rho             rhoInf; // Indicates incompressible
        rhoInf          1;      // Redundant for incompressible-flow cases

        writeControl    $/writeControl;
        writeInterval   $/writeInterval;
        log             on;
    }
}
```

```tip
The chosen solver computes these quantities on the fly; however, if the outputs
 are missing, they can be re-extracted after the simulation using: `pimpleFoam
 -postProcess`.
```

The outputs of these function objects are written to the directory
 `${FOAM_CASE}/postProcessing`. The results from the `forces` function object
  are stored in `${FOAM_CASE}/postProcessing/forces/0/force.dat`, and output of
  the `forceCoeffs` function is stored in `${FOAM_CASE}/postProcessing/forceCoeffs/0/coefficient.dat`.

---

## References

[^erzincanli2013]: Erzincanli, B. and Sahin, M. (2013). An arbitrary Lagrangian
 –Eulerian formulation for solving moving boundary problems with large
 displacement and rotations. Journal of Computational Physics 255, 660–679.
 [doi:10.1016/j.jcp.2013.08.038][erzdoi]

[^wan2006]: Wan, D. and Turek S. (2006). Fictitious Boundary and Moving Mesh
 Methods for the Numerical Simulation of Rigid Particulate Flows. Journal of
 Computational Physics 222, 28–56. [doi:10.1016/j.jcp.2006.06.002][wanddoi]

[erzdoi]: https://doi.org/10.1016/j.jcp.2013.08.038
[wanddoi]: https://doi.org/10.1016/j.jcp.2006.06.002
