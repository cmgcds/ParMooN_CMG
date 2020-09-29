// ------------------- defining parameters -------------------
d = 1;
w = 4;
a = 4;
Lu = 12*d;
Ld = 42*d;

// ------------------ defining duct boundaries -------------------
zstart = 0;
zend = a;
xstart = 0;
xend = Lu+Ld;
ystart = 0;
yend = w;

// ----------------- defining obstacle boundaries ---------
xostart = Lu-d/2;
xoend = Lu+d/2;
yostart = w/2-d/2;
yoend = w/2+d/2;

// ---------------- defining obstacle geometry -------------
lc1 = d/5;
Point(1) = {xostart,yostart,zstart,lc1};
Point(2) = {xostart,yoend,zstart,lc1};
Point(3) = {xoend,yoend,zstart,lc1};
Point(4) = {xoend,yostart,zstart,lc1};
Point(5) = {xostart,yostart,zend,lc1};
Point(6) = {xostart,yoend,zend,lc1};
Point(7) = {xoend,yoend,zend,lc1};
Point(8) = {xoend,yostart,zend,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line(21) = {1,5};
Line(22) = {2,6};
Line(23) = {3,7};
Line(24) = {4,8};

Line Loop(1) = {1,2,3,4};     // Bottom
Line Loop(2) = {5,6,7,8};     // Top
Line Loop(3) = {1,22,-5,-21}; // Left
Line Loop(4) = {-3,23,7,-24}; // Right 
Line Loop(5) = {4,21,-8,-24}; // Front
Line Loop(6) = {-2,22,6,-23}; // Back

// ------------------ defining duct geometry -------------- 
lc2 = w/5;
Point(9) = {xstart,ystart,zstart,lc2};
Point(10) = {xstart,yend,zstart,lc2};
Point(11) = {xend,yend,zstart,lc2};
Point(12) = {xend,ystart,zstart,lc2};

Point(13) = {xstart,ystart,zend,lc2};
Point(14) = {xstart,yend,zend,lc2};
Point(15) = {xend,yend,zend,lc2};
Point(16) = {xend,ystart,zend,lc2};

Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,9};

Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,13};

Line(17) = {9,13};
Line(18) = {10,14};
Line(19) = {11,15};
Line(20) = {12,16};

Line Loop(7) = {9,10,11,12};        // Bottom
Line Loop(8) = {13,14,15,16};       // Top
Line Loop(9) = {9,18,-13,-17};      // Left
Line Loop(10) = {-11,19,15,-20};    // Right 
Line Loop(11) = {12,17,-16,-20};    // Front
Line Loop(12) = {-10,18,14,-19};    // Back

// ------------------ defining surfaces -----------------

// Duct surface
Plane Surface(1) = {7,1}; // Bottom
Plane Surface(2) = {8,2}; // Top
Plane Surface(3) = {9}; // Left
Plane Surface(4) = {12}; // Back
Plane Surface(5) = {10}; // Right
Plane Surface(6) = {11}; // Front

// Obstacle surface
Plane Surface(7) = {3}; // Left
Plane Surface(8) = {6}; // Back
Plane Surface(9) = {4}; // Right
Plane Surface(10) = {5}; // Front

Surface Loop(1) = {1,2,3,4,5,6,7,8,9,10}; // ????

Physical Surface(1000) = {1};
Physical Surface(1001) = {2};
Physical Surface(1002) = {3};
Physical Surface(1003) = {4};
Physical Surface(1004) = {5};
Physical Surface(1005) = {6};
Physical Surface(1006) = {7};
Physical Surface(1007) = {8};
Physical Surface(1008) = {9};
Physical Surface(1009) = {10};

// ----------------- defining volume ---------------------
Volume(1) = {1};
Physical Volume(1) = {1};
Coherence;
