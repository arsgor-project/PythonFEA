// Gmsh project created on Thu Dec 14 15:01:58 2023
//+
//+
SetFactory("OpenCASCADE");
x1 = 0; y1 = -5; z1 = 0;
x2 = 15; y2 = -2; z2 = 25;


Point(1) = {x1, y1, z1, 1.0};
Point(2) = {x2, y2, z2, 1.0};

//+
Line(1) = {2, 1};
//+

For i In {1:19}
Rotate {{0, 0, 1}, {0, 0, 0}, i*2*Pi/20} {
  Duplicata { Curve{1}; }
}
EndFor

zn = {0, 5 , 10 , 15, 20, 25};
For i In {0:#zn[]-1}
	kn = zn[i]/(z2-z1);
	xn1 = x1 + kn*(x2-x1);
	yn1 = y1 + kn*(y2-y1);
	rn1 = Sqrt(xn1^2 + yn1^2);
	Circle(100+i) = {0, 0, zn[i], rn1, 0, 2*Pi};
EndFor
//+
BooleanFragments{ Curve{4}; Curve{5}; Curve{3}; Curve{6}; Curve{2}; Curve{7}; Curve{1}; Curve{20}; Curve{8}; Curve{19}; Curve{9}; Curve{18}; Curve{10}; Curve{17}; Curve{16}; Curve{11}; Curve{15}; Curve{12}; Curve{14}; Curve{13}; Delete; }{ Curve{100}; Curve{101}; Curve{102}; Curve{103}; Curve{104}; Curve{105}; Delete; }

//+
Physical Curve("top", 1) = {106, 105, 104, 103, 102, 101, 120, 119, 118, 117, 116, 115, 114, 113, 112, 111, 110, 109, 108, 107};
//+
Physical Curve("bottom_fix", 2) = {220, 221, 222, 223, 224, 205, 206, 225, 219, 218, 217, 216, 215, 214, 213, 212, 211, 210, 209, 208, 207};
//+
Physical Point("top_nodes", 1) = {78, 96, 66, 108, 54, 120, 36, 114, 24, 102, 12, 90, 6, 84, 18, 72, 30, 60, 42, 48};
//+
Physical Point("botot_nodes", 2) = {19, 31, 7, 49, 1, 61, 13, 73, 25, 91, 125, 37, 103, 43, 115, 55, 109, 67, 97, 79, 85};
