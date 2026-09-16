// -*- C++ -*-

// inputs
x1 =  -0.65;
y1 =  -1.2;
z1 =  0.0;
x2 =  1.1;
y2 =  0.55;
z2 =  0.0;
numPoints  = 10;

elSize = 1. / numPoints;

Point(0) = { x1, y1, z1, elSize };
Point(1) = { x2, y2, z2, elSize };

Line(10) = { 0, 1 };
