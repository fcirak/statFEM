// constants
length = 10.;
height = 10.;
depth = 10.;
meshLength = 3;
meshHeight = 3;
meshDepth = 3;

// corners
Point(11) = {-length, -.5*depth, 0.};
Point(12) = {0., -.5*depth, 0.};
Point(13) = {0., .5*depth, 0.};
Point(14) = {-length, .5*depth, 0.};
rotP11[] = Rotate { {0,1,0}, {0,0,0}, -0.25*Pi } { Duplicata { Point {11}; } };
rotP14[] = Rotate { {0,1,0}, {0,0,0}, -0.25*Pi } { Duplicata { Point {14}; } };

// edges
Line(21) = {11, 12};
Line(22) = {12, 13};
Line(23) = {13, 14};
Line(24) = {14, 11};
Line(25) = {11, rotP11};
Line(26) = {rotP11, rotP14};
Line(27) = {rotP14, 14};

// surface
Line Loop(31) = { 21, 22, 23, 24 };
Plane Surface(41) = { 31 };
Line Loop(32) = { 24, 25, 26, 27 };
Plane Surface(42) = { 32 };
surf[] = {
  41,
  Translate {length,0,0} { Duplicata { Surface {41}; } },
  Rotate { {0,1,0}, {0,0,0}, -.25*Pi } { Duplicata { Surface {41}; } },
  Rotate { {0,1,0}, {0,0,0},  -.5*Pi } { Duplicata { Surface {41}; } },
  Rotate { {0,1,0}, {0,0,0}, -.75*Pi } { Duplicata { Surface {41}; } },
  42,
  Rotate { {0,1,0}, {0,0,0}, -.25*Pi } { Duplicata { Surface {42}; } },
  Rotate { {0,1,0}, {0,0,0},  -.5*Pi } { Duplicata { Surface {42}; } },
  Rotate { {0,1,0}, {0,0,0}, -.75*Pi } { Duplicata { Surface {42}; } }
};

// mesh is transfinite
numPoints = 10; # [50, 100, ]
Mesh.CharacteristicLengthFactor = 100 / numPoints; # 0.5
For i In {0:8}
  Transfinite Surface { surf[i] };
  // make quadrilaterals
  //Recombine Surface { surf[i] };
EndFor

