SetFactory("OpenCASCADE");

// There are 4 parameters controlling a lens Rl, Rr, t, h, meshsize
// Rl: left surface radius of curvature
// Rr: right Surface radius of curvature
// t : thickness of the lens
// h : height of the lens
// meshsize: the size of the mesh 
// CALL: gmsh lens.geo -setnumber Rr 40.0 -setnumber Rl 40.0 -setnumber t 20.0 -setnumber h 25.0 -setnumber meshsize 4.0 -setnumber bl1  1.0 -setnumber bl2 3.0

// making gmsh show line numbers 
// Geometry.LineNumbers = 1;

// Mesh.CharacteristicLengthMax =  meshsize;

// sphere on the right defining the left surface of the lens
Yr = Sqrt(Rr*Rr - h*h);
Xr = 0.0;
Zr = 0.0;

// sphere on the left hand side corresponding to the right surface of the lens
Xl = 0.0;
Yl = - t + (Yr-Rr) - Rl;
Zl = 0.0;

dis = Rl - Sqrt(Rl*Rl - h*h);


Point(1) = {h,0,0};
Point(2) = {0,-Rr+Yr,0};
Point(3) = {0,Yr,0};
Circle(10) = {1,3,2};

Point(4) = {0,-t+(Yr-Rr),0};
Point(5) = {h,-t+(Yr-Rr)-dis,0};
Point(6) = {0,Yl,0};
Circle(12) = {4,6,5};

Point(7) = {h+bl2,-t+(-Rr+Yr)-dis-bl1,0};
Point(8) = {h+bl2,0+bl1,0};

Line(11) = {2,4};
Line(13) = {5,7};
Line(14) = {7,8};
Line(15) = {8,1};

Curve Loop(1) = {10, 11, 12, 13, 14, 15};
//+
Plane Surface(1) = {1};
//+
Physical Surface(1) = {1};
//+
Physical Curve(2) = {11};
//+
Physical Curve(3) = {12,13};
//+
Physical Curve(4) = {14};
//+
Physical Curve(5) = {15,10};
//+ 
Save "mesh.step";

