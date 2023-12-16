// Gmsh project created on Thu Dec 14 15:01:58 2023
//+
//+
SetFactory("OpenCASCADE");
x1 = 0; y1 = -5; z1 = 0;
x2 = 5; y2 = 3; z2 = 25;


Point(1) = {x1, y1, z1, 1.0};
Point(2) = {x2, y2, z2, 1.0};

//+
Line(1) = {2, 1};
//+

For i In {1:29}
Rotate {{0, 0, 1}, {0, 0, 0}, i*2*Pi/30} {
  Duplicata { Curve{1}; }
}
EndFor

zn = {0, 2, 5 , 10 , 12, 15, 17, 20, 22, 25};
For i In {0:#zn[]-1}
	kn = zn[i]/(z2-z1);
	xn1 = x1 + kn*(x2-x1);
	yn1 = y1 + kn*(y2-y1);
	rn1 = Sqrt(xn1^2 + yn1^2);
	Circle(100+i) = {0, 0, zn[i], rn1, 0, 2*Pi};
EndFor
