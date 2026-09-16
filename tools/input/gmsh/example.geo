// -*- C++ -*- 

// element size
elSize = 0.4;

// corner points of the rectangle 
Point(0) = {0., 0., 0., elSize};
Point(1) = {2., 0., 0., elSize};
Point(2) = {2., 1., 0., elSize};
Point(3) = {0., 1., 0., elSize};

Line(10) = {0, 1};
Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 0};

Line Loop(20) = {10, 11, 12, 13};

Plane Surface(30) = {20};

Physical Line("DirichletBoundary") = {10};

Physical Line("NeumannBoundary") = {12};

Physical Surface("Domain") = {30};



