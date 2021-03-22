// Gmsh project created on Fri Mar 12 09:54:02 2021
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 100, 0, 0, 10, 2*Pi};
//+
Physical Volume(1) = {1};
//+
Physical Surface(2) = {2};
//+
Physical Surface(1) = {3};
//+
Physical Surface(3) = {1};
