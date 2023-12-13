// Gmsh project created on Thu Nov 30 16:17:33 2023
SetFactory("OpenCASCADE");
// Gmsh project created on Wed Nov 29 16:44:53 2023
SetFactory("OpenCASCADE");

r_offset = 0;
z_offset = 0;
Rin = 5;
th = 1;
Nr = 3;
Nt = 20;
//+
Point(1) = {r_offset, z_offset, 0, 1.0};
Point(2) = {r_offset + Rin, z_offset, 0, 1.0};
Point(3) = {r_offset, z_offset + Rin, 0, 1.0};
//+
Point(4) = {r_offset + Rin + th, z_offset, 0, 1.0};
Point(5) = {r_offset, z_offset + Rin + th, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {4, 1, 5};
//+
Line(3) = {5, 3};
//+
Line(4) = {4, 2};
//+
Curve Loop(1) = {1, -3, -2, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1};
//+
Transfinite Curve {4} = Nr+1 Using Progression 1;
//+
Transfinite Curve {3} = Nr+1 Using Progression 1;
//+
Transfinite Curve {2} = Nt+1 Using Progression 1;
//+
Transfinite Curve {1} = Nt+1 Using Progression 1;
//+
Physical Curve("left", 5) = {3};
//+
Physical Curve("right", 6) = {4};
//+
Physical Curve("inner", 7) = {1};
