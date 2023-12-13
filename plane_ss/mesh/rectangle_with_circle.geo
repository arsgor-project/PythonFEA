// Gmsh project created on Wed Dec 06 14:13:27 2023
SetFactory("OpenCASCADE");
// Gmsh project created on Sun Dec 03 14:45:51 2023
//+
L = 10;
b = 5;
rad = 0.5;

Point(1) = {-L/2, -b/2, 0, 1.0};
Point(2) = {-L/2, b/2, 0, 1.0};
Point(3) = {L/2, -b/2, 0, 1.0};
Point(4) = {L/2, b/2, 0, 1.0};

Point(5) = {0, -b/2, 0, 1.0};
Point(6) = {0, b/2, 0, 1.0};
Point(7) = {0, -rad, 0, 1.0};
Point(8) = {0, rad, 0, 1.0};
Point(9) = {0, 0, 0, 1.0};

Point(10) = {-L/2, 0, 0, 1.0};
Point(11) = {L/2, 0, 0, 1.0};

Point(12) = {-rad, 0, 0, 1.0};
Point(13) = {rad, 0, 0, 1.0};
Point(14) = {rad, 0, 0, 1.0};
Point(15) = {rad, 0, 0, 1.0};

//+
Line(1) = {2, 10};
//+
Line(2) = {1, 10};
//+
Line(3) = {1, 5};
//+
Line(4) = {5, 3};
//+
Line(5) = {3, 11};
//+
Line(6) = {11, 4};
//+
Line(7) = {4, 6};
//+
Line(8) = {6, 2};
//+
Line(9) = {6, 8};
//+
Line(10) = {11, 13};
//+
Line(11) = {5, 7};
//+
Line(12) = {10, 12};
//+
Circle(13) = {12, 9, 8};
//+
Circle(14) = {13, 9, 8};
//+
Circle(15) = {7, 9, 13};
//+
Circle(16) = {7, 9, 12};
//+
Curve Loop(1) = {8, 1, 12, 13, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 9, -14, -10, 6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -15, -11, 4, 5};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, -16, -11, -3, 2};
//+
Plane Surface(4) = {4};
//+
Coherence;
//+
Physical Curve("fix", 1) = {1, 2};
//+
Physical Curve("load", 5) = {14, 15};
//+
Physical Point("exclude", 1) = {9};
