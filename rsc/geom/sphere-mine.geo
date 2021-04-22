// SetFactory("OpenCASCADE");

// Set the geometry order (1, 2, ..., 9)
order = 3;

// Set the element type (3 - triangles, 4 - quadrilaterals)
type = 3;

// Set radius
rad = 1;

// Set point mesh density
lc = 0.1;

// Point
Point(1) = {0, 0, 0, lc};
Point(2) = {rad, 0, -0, lc};
Point(3) = {-rad, 0, 0, lc};
Point(4) = {0, rad, 0, lc};
Point(5) = {0, -rad, 0, lc};
Point(6) = {0, 0, rad, lc};
Point(7) = {0, 0, -rad, lc};

// Frame
Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {2, 1, 6};
Circle(6) = {6, 1, 3};
Circle(7) = {3, 1, 7};
Circle(8) = {7, 1, 2};
Circle(9) = {4, 1, 6};
Circle(10) = {6, 1, 5};
Circle(11) = {5, 1, 7};
Circle(12) = {7, 1, 4};

// Surfaces
Curve Loop(1) = {1, 9, -5};
Surface(1) = {1};
Curve Loop(3) = {6, -2, 9};
Surface(2) = {3};
Curve Loop(5) = {2, 7, 12};
Surface(3) = {5};
Curve Loop(7) = {12, -1, -8};
Surface(4) = {7};
Curve Loop(9) = {8, -4, 11};
Surface(5) = {9};
Curve Loop(11) = {7, -11, -3};
Surface(6) = {11};
Curve Loop(13) = {6, 3, -10};
Surface(7) = {13};
Curve Loop(15) = {10, 4, 5};
Surface(8) = {15};

// Quadrilaterals mesh
If (type == 4)
    Recombine Surface {1};
    Recombine Surface {2};
    Recombine Surface {3};
    Recombine Surface {4};
    Recombine Surface {5};
    Recombine Surface {6};
    Recombine Surface {7};
    Recombine Surface {8};
EndIf

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};
Physical Surface(8) = {8};

// Generate 2D mesh
Mesh 2;
SetOrder order;
Mesh.MshFileVersion = 2.2;

