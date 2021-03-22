lc = .1;

Point(1) = {0, -2, 0, lc};
Point(2) = {15, -2,  0, lc} ;
Point(3) = {15, 2, 0, lc} ;
Point(4) = {0,  2, 0, lc} ;

//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Physical Line("1") = {4};
Physical Line("2") = {2};
Physical Line("3") = {1};
Physical Line("4") = {3};


Physical Surface("1") = {1} ;


