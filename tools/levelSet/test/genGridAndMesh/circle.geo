// -*- C++ -*-

// inputs
ox = 0.0;
oy = 0.0;
oz = 0.0;
r  = 0.5;
numPoints = 25;
elSize = 1. / numPoints;

Point(1) = { ox,     oy,     oz, elSize };
Point(2) = { ox - r, oy,     oz, elSize };
Point(3) = { ox,     oy - r, oz, elSize };
Point(4) = { ox + r, oy,     oz, elSize };
Point(5) = { ox	,    oy + r, oz, elSize };

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1,2,3,4};
