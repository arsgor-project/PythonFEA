// Gmsh project created on Wed Dec 13 18:13:23 2023
SetFactory("OpenCASCADE");
//+
a = 1;
Point(1) = {0.0, 0.0, 0, 1.0};
Point(2) = {a, 0.0, 0, 1.0};
Point(3) = {2*a, 0.0, 0, 1.0};
Point(4) = {0.5*a, a, 0, 1};
Point(5) = {1.5*a, a, 0, 1}; //+
Line(1) = {1, 4};
//+
Line(2) = {4, 2};
//+
Line(3) = {2, 5};
//+
Line(4) = {5, 3};
//+
Line(5) = {3, 2};
//+
Line(6) = {2, 1};
//+
Line(7) = {4, 5};
//+
Physical Point("Fix", 1) = {1};
//+
Physical Point("Fix2", 2) = {3};
//+
Physical Point("Load", 3) = {4};
