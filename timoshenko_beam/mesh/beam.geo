// Gmsh project created on Tue Dec 12 19:43:59 2023
SetFactory("OpenCASCADE");
//+
L = 10;
Point(1) = {0.0, 0, 0, 1.0};
Point(2) = {L,0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Physical Point("Left", 1) = {1};
//+
Physical Point("Right", 2) = {2};
//+
Physical Curve("UDL", 3) = {1};
