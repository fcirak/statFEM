// Gmsh project created on Tue May 11 10:02:04 2010
length = 100;   // length of each side of plate
numDiv = 5;    // number of divisions along each side of plate

// corners
Point(1) = {0, 0, 0};
Point(2) = {length, 0, 0};
Point(3) = {0, length, 0};
Point(4) = {length, length, 0};

// edges
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

// surface
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// element subdivisions on edges
Transfinite Line {1, 2, 3, 4} = numDiv;
// surface is transfinite
Transfinite Surface {6};
// make quadrilaterals
//Recombine Surface{6};

