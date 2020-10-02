SetFactory("OpenCASCADE");

// There are 4 parameters controlling a lens Rl, Rr, t, h, meshsize
// Rl: left surface radius of curvature
// Rr: right Surface radius of curvature
// t : thickness of the lens
// h : height of the lens
// meshsize: the size of the mesh 
// CALL: gmsh lens.geo -setnumber Rr 40.0 -setnumber Rl 40.0 -setnumber t 20.0 -setnumber h 25.0 -setnumber meshsize 4.0

// making gmsh show line numbers 
Geometry.LineNumbers = 1;


// sphere on the right defining the left surface of the lens
Yr = Sqrt(Rr*Rr - h*h);
Xr = 0.0;
Zr = 0.0;

// sphere on the left hand side corresponding to the right surface of the lens
Xl = 0.0;
Yl = t - (Rr-Yr) - Rl;
Zl = 0.0;

dis = Rl - Sqrt(Rl*Rl - h*h);


Point(1) = {h,0,0};
Point(2) = {0,Rr-Yr,0};
Point(3) = {0,-Yr,0};
Circle(10) = {1,3,2};

Point(4) = {0,-t+(Rr-Yr),0};
Point(5) = {h,-t+(Rr-Yr)+dis,0};
Point(6) = {0,-Yl,0};
Circle(12) = {4,6,5};

Line(11) = {2,4};
Line(13) = {5,1};

Curve Loop(1) = {10, 11, 12, 13};
//+
Plane Surface(1) = {1};
//+
Physical Surface(1) = {1};
//+
Physical Curve(2) = {11};
//+
Physical Curve(3) = {12};
//+
Physical Curve(4) = {13};
//+
Physical Curve(5) = {10};
//+ 
Save "mesh.step";

