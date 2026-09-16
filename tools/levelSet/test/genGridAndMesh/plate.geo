// -*- C++ -*-

// inputs
ox = 0.0;
oy = 0.0;
oz = 0.5;
lx = 0.4;
ly = 0.4;
numPoints = 10;
nx = numPoints;
ny = numPoints;

// corner points of a rectangle
Point(1) = {ox,      oy,      oz};
Point(2) = {ox + lx, oy,      oz};
Point(3) = {ox + lx, oy + ly, oz};
Point(4) = {ox,      oy + ly, oz};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line "*" = numPoints;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
