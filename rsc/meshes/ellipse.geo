Mesh.Algorithm = 6;

// mesh density
lc = 0.1;
// x axis
a = 1;
// y axis
b = 1.5;
// z axis
c = 2;

// Points
Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {a,0.0,0.0,lc};
Point(3) = {0,b,0.0,lc};
Point(4) = {-a,0,0.0,lc};
Point(5) = {0,-b,0.0,lc};
Point(6) = {0,0,-c,lc};
Point(7) = {00,0,c,lc};

// Lines
Ellipse(5) = {7, 1, 4, 4};
Ellipse(6) = {7, 1, 2, 2};
Ellipse(7) = {6, 1, 2, 2};
Ellipse(8) = {6, 1, 4, 4};
Ellipse(9) = {3, 1, 4, 4};
Ellipse(10) = {3, 1, 2, 2};
Ellipse(11) = {5, 1, 2, 2};
Ellipse(12) = {5, 1, 4, 4};
Ellipse(13) = {7, 1, 3, 3};
Ellipse(14) = {6, 1, 3, 3};
Ellipse(15) = {6, 1, 5, 5};
Ellipse(16) = {7, 1, 5, 5};

// Surfaces and volumes
Curve Loop(1) = {10, -6, 13};
Surface(1) = {1};
Curve Loop(2) = {14, 10, -7};
Surface(2) = {2};
Curve Loop(3) = {14, 9, -8};
Surface(3) = {3};
Curve Loop(4) = {9, -5, 13};
Surface(4) = {4};
Curve Loop(5) = {7, -11, -15};
Surface(5) = {5};
Curve Loop(6) = {6, -11, -16};
Surface(6) = {6};
Curve Loop(7) = {16, 12, -5};
Surface(7) = {7};
Curve Loop(8) = {8, -12, -15};
Surface(8) = {8};

Surface Loop(9) = {1, 2, 3, 4, 7, 6, 5, 8};
Volume(10) = {9};
Physical Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Volume(2) = {10};

// Generate 2D mesh
Mesh 2;
Mesh.MshFileVersion = 2.2;