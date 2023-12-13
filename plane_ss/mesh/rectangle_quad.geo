// Gmsh project created on Wed Nov 29 16:44:53 2023
SetFactory("OpenCASCADE");
Lx = 10;
Ly = 5;
//+
Point(1) = {0, 0,   0, 1.0};
Point(2) = {Lx, 0,  0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0, Ly,  0, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1};
//+
Physical Curve("fix", 1) = {4};
//+
Transfinite Curve {4, 1, 2, 3} = 21 Using Progression 1;
//+
Physical Curve("load", 5) = {2};
