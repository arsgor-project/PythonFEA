// Gmsh project created on Wed Dec 06 15:23:30 2023
SetFactory("OpenCASCADE");
// Gmsh project created on Sun Dec 03 14:45:51 2023
//+
L = 20;
b = 20;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {L,0, 0, 1.0};
Point(3) = {L, b, 0, 1.0};
Point(4) = {L+b, b, 0, 1.0};

//+
Line(1) = {3, 4};
//+
Line(2) = {4, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 3};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
//Transfinite Surface {1};
//+
//Transfinite Curve {4, 3, 2, 1} = 3 Using Progression 1;
//+
Physical Curve("fix", 1) = {3};
//+
Physical Curve("load", 5) = {1};
//+
Point(5) = {30, 0, 0, 1.0};
