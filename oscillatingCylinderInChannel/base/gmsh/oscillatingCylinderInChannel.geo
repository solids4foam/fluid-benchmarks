//--------------------------------------------------
// Global scaling factor
//--------------------------------------------------
// The "meshLevel" can be set as a command line option,
// > gmsh -3 -format msh2 -setnumber meshLevel 2.0 oscillatingCylinderInChannel.geo
// A larger value corresponds to a finer mesh
// Increasing the meshLevel by 1.0 corresponds to halving the average mesh
// spacing, i.e., this will increase the number of cells by a factor of 4 (as
// this is a 2-D case)
If (!Exists(meshLevel))
  meshLevel = 1.0; // default value
EndIf

// Set the global scale factor
sf = 2.0 / (2.0 ^ meshLevel);

//--------------------------------------------------
// Geometry parameters
//--------------------------------------------------

// Cylinder centre and radius
x0 = 1.1;
y0 = 0.2;
r  = 0.05;

// Domain bounds
x1 = 0.0;
x2 = 2.2;
y1 = 0.0;
y2 = 0.41;
z  = 0.1;

//--------------------------------------------------
// Mesh size parameters
//--------------------------------------------------
cylinderDeltaX   = 0.01 * sf;
outerWallsDeltaX = 0.01 * sf;

Mesh.CharacteristicLengthMin = 0.01 * sf;
Mesh.CharacteristicLengthMax = 0.05 * sf;

//--------------------------------------------------
// Outer rectangle
//--------------------------------------------------
Point(1) = {x1, y1, 0, outerWallsDeltaX};
Point(2) = {x2, y1, 0, outerWallsDeltaX};
Point(3) = {x2, y2, 0, outerWallsDeltaX};
Point(4) = {x1, y2, 0, outerWallsDeltaX};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//--------------------------------------------------
// Cylinder
//--------------------------------------------------
Point(5) = {x0, y0, 0, cylinderDeltaX};
Point(6) = {x0 + r, y0, 0, cylinderDeltaX};
Point(7) = {x0, y0 + r, 0, cylinderDeltaX};
Point(8) = {x0 - r, y0, 0, cylinderDeltaX};
Point(9) = {x0, y0 - r, 0, cylinderDeltaX};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

//--------------------------------------------------
// Surface with cylinder hole
//--------------------------------------------------
Curve Loop(10) = {1, 2, 3, 4}; // outer rectangle
Curve Loop(11) = {5, 6, 7, 8}; // cylinder
Plane Surface(12) = {10, 11};

//--------------------------------------------------
// Mesh control: distance-based refinement
//--------------------------------------------------
Mesh.CharacteristicLengthFromPoints     = 0;
Mesh.CharacteristicLengthFromCurvature  = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

// 1) Walls
Field[1] = Distance;
Field[1].EdgesList = {1, 2, 3, 4};
Field[1].Sampling  = 800 / sf;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].LcMin   = 0.010 * sf;
Field[2].LcMax   = 0.050 * sf;
Field[2].DistMin = 0; //0.20 * sf;
Field[2].DistMax = 0.1; //1.0 * sf;

// 2) Cylinder
Field[3] = Distance;
Field[3].EdgesList = {5, 6, 7, 8};

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].LcMin   = 0.006 * sf;
Field[4].LcMax   = 0.050 * sf;
Field[4].DistMin = 0; //0.15 * sf;
Field[4].DistMax = 0.1; //0.5 * sf;

// 3) Combine
Field[5] = Min;
Field[5].FieldsList = {2, 4};
Background Field = 5;

//--------------------------------------------------
// Extrude to 3D (triangles -> prisms)
//--------------------------------------------------
out[] = Extrude {0, 0, z} {
    Surface{12};
    Layers{1};
    Recombine; // still produces prisms from tris
};

//--------------------------------------------------
// Physical groups
//--------------------------------------------------
Physical Volume("fluidVolume") = {out[1]};

Physical Surface("cylinder") = {out[6], out[7], out[8], out[9]};
Physical Surface("left")     = {37};
Physical Surface("right")    = {29};
Physical Surface("up")       = {33};
Physical Surface("down")     = {25};
Physical Surface("back")     = {12};
Physical Surface("front")    = {54};
