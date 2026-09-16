// -*- C++ -*- 

// element size
elSize = 0.01;

xmax = 2.;
ymax = 0.5;

xc = 0.8;
yc = 0.2;
r = 0.1;

// corner points of the rectangle 
Point(0) = {  0.,   0., 0., elSize};
Point(1) = {xmax,   0., 0., elSize};
Point(2) = {xmax, ymax, 0., elSize};
Point(3) = {  0., ymax, 0., elSize};

// outer box
Line(10) = {0, 1};
Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 0};

// circle points
Point(4) = { xc, yc, 0., elSize };
Point(5) = { xc - r, yc, 0., elSize };
Point(6) = { xc, yc + r, 0., elSize };
Point(7) = { xc + r, yc, 0., elSize };
Point(8) = { xc, yc - r, 0., elSize };

// embedded circle
Circle(20) = {5,4,6};
Circle(21) = {6,4,7};
Circle(22) = {7,4,8};
Circle(23) = {8,4,5};


Line Loop(20) = {10, 11, 12, 13};
Line Loop(21) = {20, 21, 22, 23};

Plane Surface(30) = {20,21};

Physical Line("Flow") = {11,13};
Physical Line("Wall") = {10,12,20,21,22,23};
Physical Surface("Domain") = {30};



