// Gmsh project created on Wed Dec 13 19:01:55 2023
SetFactory("OpenCASCADE");
a = 3;

//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {a, 0, 0, 1.0};
Point(3) = {2*a, 0, 0, 1.0};
Point(4) = {0, a, 0, 1.0};
Point(5) = {a, a, 0, 1.0};
Point(6) = {2*a, a, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 6};
//+
Line(4) = {6, 5};
//+
Line(5) = {5, 4};
//+
Line(6) = {4, 1};
//+
Line(7) = {2, 5};
//+
Line(8) = {2, 6};
//+
Line(9) = {3, 5};
//+
Line(10) = {2, 4};
//+
Line(11) = {1, 5};
//+
Physical Point("Fix", 1) = {1};
//+
Physical Point("Fix2", 2) = {3};
//+
Physical Point("Load", 3) = {4, 5, 6};
