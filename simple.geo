Merge "MeshesMSH/sphere.msh";
Surface Loop(43) = {0};


tol = .05;
avg_len = 7;
lc = 10;


xy_res = 141.7/512;
z_res = 120/151;
xy_bound = 141.7;
z_bound = 120;

Point(1) = {0, 0, 0, lc};
Point(2) = {0, xy_bound, 0, lc};
Point(3) = {xy_bound, xy_bound, 0, lc};
Point(4) = {xy_bound, 0, 0, lc};
Point(5) = {0, 0, z_bound, lc};
Point(6) = {0, xy_bound, z_bound, lc};
Point(7) = {xy_bound, xy_bound, z_bound, lc};
Point(8) = {xy_bound, 0, z_bound, lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {1,5};
Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {7,8};
Line(9) = {8,5};
Line(10) = {2, 6};
Line(11) = {3,7};
Line(12) = {4,8};
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {1, 10, -6, -5};
Curve Loop(3) = {2, 11, -7, -10};
Curve Loop(4) = {3, 12, -8, -11};
Curve Loop(5) = {6, 7, 8, 9};
Curve Loop(6) = {5, -9, -12, 4};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Surface Loop(42) = {1:6};

Volume(0) = {42, 43};

Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

Mesh.CharacteristicLengthMax = avg_len + tol;
Mesh.CharacteristicLengthMin = avg_len - tol;
Mesh.CharacteristicLengthFactor = 1;

