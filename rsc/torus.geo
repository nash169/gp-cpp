// Parameters
lc = DefineNumber[ 0.2, Name "Parameters/lc" ];
rad_ext = DefineNumber[ 3, Name "Parameters/rad_ext" ];
rad_in = DefineNumber[ 1, Name "Parameters/rad_in" ];

// Set the geometry order (1, 2, ..., 9)
order = DefineNumber[ 3, Name "Parameters/order" ];

// Set the element type (3 - triangles, 4 - quadrilaterals)
type = DefineNumber[ 3, Name "Parameters/type" ];

// Points
Point(1) = {0, 0, 0, lc};
Point(22) = {0, 0, rad_in, lc};
Point(23) = {0, 0, -rad_in, lc};

Point(2) = {rad_ext, 0, 0, lc};
Point(3) = {0, rad_ext, 0, lc};
Point(4) = {-rad_ext, 0, 0, lc};
Point(5) = {0, -rad_ext, 0, lc};

Point(6) = {rad_ext + rad_in, 0, 0, lc};
Point(7) = {rad_ext - rad_in, 0, 0, lc};
Point(8) = {rad_ext, 0, rad_in, lc};
Point(9) = {rad_ext, 0, -rad_in, lc};

Point(10) = {-rad_ext + rad_in, 0, 0, lc};
Point(11) = {-rad_ext - rad_in, 0, 0, lc};
Point(12) = {-rad_ext, 0, rad_in, lc};
Point(13) = {-rad_ext, 0, -rad_in, lc};

Point(14) = {0, rad_ext + rad_in, 0, lc};
Point(15) = {0, rad_ext - rad_in, 0, lc};
Point(16) = {0, rad_ext, rad_in, lc};
Point(17) = {0, rad_ext, -rad_in, lc};

Point(18) = {0, -rad_ext + rad_in, 0, lc};
Point(19) = {0, -rad_ext - rad_in, 0, lc};
Point(20) = {0, -rad_ext, rad_in, lc};
Point(21) = {0, -rad_ext, -rad_in, lc};

// Curves
Circle(1) = {20, 22, 12};
Circle(2) = {12, 22, 16};
Circle(3) = {16, 22, 8};
Circle(4) = {8, 22, 20};
Circle(5) = {13, 23, 17};
Circle(6) = {13, 23, 21};
Circle(7) = {21, 23, 9};
Circle(8) = {9, 23, 17};
Circle(9) = {19, 1, 11};
Circle(10) = {11, 1, 14};
Circle(11) = {14, 1, 6};
Circle(12) = {6, 1, 19};
Circle(13) = {15, 1, 7};
Circle(14) = {7, 1, 18};
Circle(15) = {18, 1, 10};
Circle(16) = {10, 1, 15};
Circle(17) = {20, 5, 19};
Circle(18) = {19, 5, 21};
Circle(19) = {21, 5, 18};
Circle(20) = {18, 5, 20};
Circle(21) = {17, 3, 14};
Circle(22) = {14, 3, 16};
Circle(23) = {16, 3, 15};
Circle(24) = {15, 3, 17};
Circle(25) = {12, 4, 11};
Circle(26) = {11, 4, 13};
Circle(27) = {10, 4, 13};
Circle(28) = {12, 4, 10};
Circle(29) = {8, 2, 7};
Circle(30) = {7, 2, 9};
Circle(31) = {9, 2, 6};
Circle(32) = {6, 2, 8};

// Surfaces and volumes
Curve Loop(1) = {2, -22, -10, -25};
Surface(1) = {1};
Curve Loop(2) = {1, 25, -9, -17};
Surface(2) = {2};
Curve Loop(3) = {17, -12, 32, 4};
Surface(3) = {3};
Curve Loop(4) = {-32, 3, 22, -11};
Surface(4) = {4};
Curve Loop(5) = {29, 14, 20, -4};
Surface(5) = {5};
Curve Loop(6) = {15, -28, -1, -20};
Surface(6) = {6};
Curve Loop(7) = {16, -23, -2, 28};
Surface(7) = {7};
Curve Loop(8) = {23, 13, -29, -3};
Surface(8) = {8};
Curve Loop(9) = {-5, -21, 10, -26};
Surface(9) = {9};
Curve Loop(10) = {26, 6, -18, 9};
Surface(10) = {10};
Curve Loop(11) = {18, 7, 31, 12};
Surface(11) = {11};
Curve Loop(12) = {-31, 11, 21, 8};
Surface(12) = {12};
Curve Loop(13) = {5, -24, -16, 27};
Surface(13) = {13};
Curve Loop(14) = {-27, -6, -19, -15};
Surface(14) = {14};
Curve Loop(15) = {19, -14, 30, -7};
Surface(15) = {15};
Curve Loop(16) = {-30, -8, 24, -13};
Surface(16) = {16};

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
    Recombine Surface {9};
    Recombine Surface {10};
    Recombine Surface {11};
    Recombine Surface {12};
    Recombine Surface {13};
    Recombine Surface {14};
    Recombine Surface {15};
    Recombine Surface {16};
EndIf

Surface Loop(1) = {16, 15, 14, 13, 9, 12, 11, 10, 2, 6, 7, 8, 5, 3, 4, 1};
Volume(1) = {1};
Physical Surface(33) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
Physical Volume(34) = {1};

// Generate 2D mesh
Mesh 2;
SetOrder order;
Mesh.MshFileVersion = 2.2;
