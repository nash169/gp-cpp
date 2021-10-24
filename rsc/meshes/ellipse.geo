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
Point(7) = {0,0,c,lc};

// Lines
Ellipse(1) = {7, 1, 4, 4};
Ellipse(2) = {7, 1, 2, 2};
Ellipse(3) = {6, 1, 2, 2};
Ellipse(4) = {6, 1, 4, 4};
Ellipse(5) = {3, 1, 4, 4};
Ellipse(6) = {3, 1, 2, 2};
Ellipse(7) = {5, 1, 2, 2};
Ellipse(8) = {5, 1, 4, 4};
Ellipse(9) = {7, 1, 3, 3};
Ellipse(10) = {6, 1, 3, 3};
Ellipse(11) = {6, 1, 5, 5};
Ellipse(12) = {7, 1, 5, 5};

Curve Loop(1) = {-9, -6, 2};
Surface(1) = {1};
Curve Loop(2) = {6, -3, 10};
Surface(2) = {2};
Curve Loop(3) = {-10, -5, 4};
Surface(3) = {3};
Curve Loop(4) = {5, -1, 9};
Surface(4) = {4};
Curve Loop(5) = {-4, 8, 11};
Surface(5) = {5};
Curve Loop(6) = {-8, 1, -12};
Surface(6) = {6};
Curve Loop(7) = {12, 7, -2};
Surface(7) = {7};
Curve Loop(8) = {-7, 3, -11};
Surface(8) = {8};

// Generate 2D mesh
Mesh 2;
Mesh.MshFileVersion = 2.2;
