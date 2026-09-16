// inputs
ox = 0.0;
oy = 0.0;
oz = 0.0;
r = 0.5;
numPoints = 23;
elSize = 1. / numPoints;

Point(1) = {ox    , oy    , oz    , elSize};
Point(2) = {ox + r, oy    , oz    , elSize};
Point(3) = {ox    , oy + r, oz    , elSize};
Point(4) = {ox    , oy    , oz + r, elSize};
Point(5) = {ox - r, oy    , oz    , elSize};
Point(6) = {ox - r, oy    , oz    , elSize};
Point(7) = {ox    , oy - r, oz    , elSize};
Point(8) = {ox    , oy    , oz - r, elSize};

Circle(1) = {4, 1, 2};
Circle(2) = {4, 1, 3};
Circle(3) = {2, 1, 3};
Circle(6) = {5, 1, 3};
Circle(7) = {5, 1, 7};
Circle(8) = {7, 1, 2};
Circle(9) = {7, 1, 4};
Circle(10) = {5, 1, 4};
Circle(17) = {8, 1, 7};
Circle(18) = {8, 1, 2};
Circle(19) = {8, 1, 3};
Circle(20) = {8, 1, 5};

Line Loop(5) = {3, -2, 1};
Ruled Surface(5) = {5};
Line Loop(12) = {-9, -1, 8};
Ruled Surface(12) = {12};
Line Loop(14) = {9, -10, 7};
Ruled Surface(14) = {14};
Line Loop(16) = {10, 2, -6};
Ruled Surface(16) = {16};
Line Loop(22) = {-20,-7, 17};
Ruled Surface(22) = {22};
Line Loop(24) = {20, 6, -19};
Ruled Surface(24) = {24};
Line Loop(26) = {19, -3, -18};
Ruled Surface(26) = {26};
Line Loop(28) = {18, -8, -17};
Ruled Surface(28) = {28};
