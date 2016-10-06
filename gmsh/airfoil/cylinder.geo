
Point(1) = {0,0,0,0.0,0.01};
Point(2) = {1,0,0,0.0,0.01};
Point(3) = {0,1,0,0.0,0.01};
Point(4) = {-1,0,0,0.0,0.01};
Point(5) = {0,-1,0,0.0,0.01};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Point(6) = {-5, 5, 0, 1.0,0.1};
Point(7) = {10, 5, 0, 1.0,0.1};
Point(8) = {10, -5, 0, 1.0,0.1};
Point(9) = {-5, -5, 0, 1.0,0.1};
Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 9};
Line(11) = {9, 6};
Line Loop(12) = {2, 3, 4, 1};
Line Loop(13) = {8, 9, 10, 11};
Plane Surface(14) = {12, 13};

Characteristic Length {4, 3, 2, 5, 8, 7, 6, 9} = 0.1;
Characteristic Length {6, 9} = 0.5;
Characteristic Length {7, 8} = 0.5;


Physical Line("inlet",1) = {8};
Physical Line("outlet",2) = {10};
Physical Line("wall",4) = {9,11};
Physical Line("cylinder",3) = {1,2,3,4};
Physical Surface("body",999) = {14};
