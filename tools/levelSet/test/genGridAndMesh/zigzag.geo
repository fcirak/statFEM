// -*- C++ -*- 

// inputs
ox = -0.5;
oy = -0.3;
a  = 0.8;
b  = 0.8;
c  = 0.4;
N  = 10;						
elSize = 1. / N;				

// point coordates
p1x = ox;
p1y = oy;
Point(1) = {p1x, p1y, 0., elSize};
p2x = p1x + b;
p2y = p1y - a;
Point(2) = {p2x, p2y, 0., elSize};
p3x = p1x + 2. * b;
p3y = p1y;
Point(3) = {p3x, p3y, 0., elSize};
p4x = p3x - c;
p4y = p3y + c;
Point(4) = {p4x, p4y, 0., elSize};
p5x = p3x - 2. * c;
p5y = p3y; 
Point(5) = {p5x, p5y, 0., elSize};

// lines
Line(1) = {1, 2};				
Line(2) = {2, 3};				
Line(3) = {3, 4};
Line(4) = {4, 5};
