// Gmsh project created on Thu Apr 13 15:11:37 2023
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Line(1) = {2, 4};
//+
Line(2) = {4, 1};
//+
Line(3) = {2, 1};
//+
Line(4) = {2, 3};
//+
Line(5) = {4, 3};
//+
Curve Loop(1) = {1, 5, -4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -2, -1};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {1} = 10 Using Progression 1;
//+
Transfinite Curve {3, 2, 5, 4} = 10 Using Progression 3;
