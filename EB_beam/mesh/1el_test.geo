// Gmsh project created on Thu Dec 14 18:10:45 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0,  0, 0, 1.0};
Point(2) = {1,  1, 1, 1.0};
//+
Line(1) = {2, 1};
//+
Physical Point("fix", 2) = {2};
//+
Physical Point("load", 1) = {1};
