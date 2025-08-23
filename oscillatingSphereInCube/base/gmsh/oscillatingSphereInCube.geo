// Sphere oscillating in a cubic cavity (geometry + mesh)
// L = 2D cube; centred sphere of diameter D is removed (fluid)
// A global scale factor is computed as:
//    sf = 2 / (2^meshLevel)
// where `meshLevel` can be set on the CLI and defaults to 1.
// Example:
//   gmsh -3 -format msh2 -setnumber meshLevel 3 sphere_in_cube.geo
// ------------------------------------------------------------
SetFactory("OpenCASCADE");

// ---------- controls (override with -setnumber) ----------
If (!Exists(meshLevel)) meshLevel = 1; EndIf
Printf("mesh level = %g", meshLevel);

// geometry sizes (nondimensional base: D = 1)
D  = 1.0;
R  = 0.5*D;          // sphere radius
L  = 2.0*D;          // cube edge length
h  = 0.125*D;        // oscillation amplitude used in the paper (for reference)

// global mesh scale factor (make sure this is floating-point)
sf = 2.0 / (2.0^meshLevel);
Printf("scale factor = %g", sf);

// base mesh sizes (scaled by sf)
lc_far    = 0.20*D*sf;     // away from the sphere
lc_near   = 0.05*D*sf;     // near the sphere

// small geometric tolerance for selecting faces
eps = 1e-6*L;

// ---------- geometry ----------
xmin = -0.5*L;  xmax =  0.5*L;
ymin = -0.5*L;  ymax =  0.5*L;
zmin = -0.5*L;  zmax =  0.5*L;

// Make the solid cube, then subtract the sphere to get the fluid region
Box(1)    = {xmin, ymin, zmin, L, L, L};
Sphere(2) = {0.0, 0.0, 0.0, R};

// Boolean difference: fluid = cube  sphere
out[] = BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; };

// ---------- identify boundary pieces ----------
allBdry[] = Boundary{ Volume{out[0]}; };

// The spherical inner wall (unique surface inside a small BB around the origin)
sphSurf[] = Surface In BoundingBox{-R-eps, -R-eps, -R-eps, R+eps, R+eps, R+eps};

// Outer cube faces via thin bounding boxes around each plane
xMinSurf[] = Surface In BoundingBox{xmin-eps, ymin-eps, zmin-eps, xmin+eps, ymax+eps, zmax+eps};
xMaxSurf[] = Surface In BoundingBox{xmax-eps, ymin-eps, zmin-eps, xmax+eps, ymax+eps, zmax+eps};
yMinSurf[] = Surface In BoundingBox{xmin-eps, ymin-eps, zmin-eps, xmax+eps, ymin+eps, zmax+eps};
yMaxSurf[] = Surface In BoundingBox{xmin-eps, ymax-eps, zmin-eps, xmax+eps, ymax+eps, zmax+eps};
zMinSurf[] = Surface In BoundingBox{xmin-eps, ymin-eps, zmin-eps, xmax+eps, ymax+eps, zmin+eps};
zMaxSurf[] = Surface In BoundingBox{xmin-eps, ymin-eps, zmax-eps, xmax+eps, ymax+eps, zmax+eps};

// ---------- sizing field: refine near the sphere ----------
Field[1] = Distance;
Field[1].SurfacesList = {sphSurf[]};  // distance to the inner sphere surface
Field[1].NumPointsPerCurve = 100;

Field[2] = Threshold;
Field[2].InField = 1;
// Closer than DistMin -> SizeMin; farther than DistMax -> SizeMax; smooth in between
Field[2].DistMin = 0.20*D;
Field[2].DistMax = 0.60*D;
Field[2].SizeMin = lc_near;
Field[2].SizeMax = lc_far;

Background Field = 2;

// (Optional) a little optimisation
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

// ---------- physical groups ----------
Physical Volume("Fluid") = {out[]};

Physical Surface("sphere")  = {sphSurf[]};
Physical Surface("xmin")    = {xMinSurf[]};
Physical Surface("xmax")    = {xMaxSurf[]};
Physical Surface("ymin")    = {yMinSurf[]};
Physical Surface("ymax")    = {yMaxSurf[]};
Physical Surface("zmin")    = {zMinSurf[]};
Physical Surface("zmax")    = {zMaxSurf[]};}