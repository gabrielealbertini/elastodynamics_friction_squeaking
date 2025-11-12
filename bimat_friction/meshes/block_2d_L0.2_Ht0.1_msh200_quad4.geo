// use multiples of 10
Nb_elems = 200;

x_size = 0.2;
y_size = 0.1;

h_def = x_size/Nb_elems;

Point(1) = {0,           0, 0, h_def};
Point(2) = {x_size,      0, 0, h_def};
Point(3) = {x_size, y_size, 0, h_def};
Point(4) = {0,      y_size, 0, h_def};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Line("slider_top") = {3};
Physical Line("slider_bottom") = {1};
Physical Line("slider_left") = {4};
Physical Line("slider_right") = {2};

Physical Surface("slider") = {1};

//uniform mesh
Transfinite Surface "*";
Recombine Surface "*"; // quadrangle

Mesh.SecondOrderIncomplete = 1;