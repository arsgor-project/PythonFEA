// Gmsh project created on Tue Dec 12 14:31:49 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0.5, 0.0, 0, 1.0};
Point(2) = {1.0, 0.0, 0, 1.0};
Point(3) = {1.0, 1.0, 0, 1.0};
Point(4) = {0.5, 1, 0, 1.0};


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
Physical Curve("load", 7) = {4};
//+
Physical Curve("uz_fix", 5) = {1, 3};
//+
Transfinite Surface {1};
//+
Transfinite Curve {1, 2, 3, 4} = 2 Using Progression 1;
