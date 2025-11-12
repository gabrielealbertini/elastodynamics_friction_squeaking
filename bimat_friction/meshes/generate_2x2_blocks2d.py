#!/usr/bin/env python
from __future__ import print_function

mesh_file="""
// use multiples of 10
Nb_elems = {Nx};

x_size = {Lx};
y_size = {Ly}; 


"""
mesh_file_block="""
h_def = x_size/Nb_elems;

b_x_size = x_size;
b_y_size = y_size;

// top left

Point(1) = {0,           0, 0, h_def};
Point(2) = {x_size,      0, 0, h_def};
Point(3) = {x_size, y_size, 0, h_def};
Point(4) = {0,      y_size, 0, h_def};

// top right

Point(5) = {x_size,                   0, 0, h_def};
Point(6) = {x_size + b_x_size,        0, 0, h_def};
Point(7) = {x_size + b_x_size,   y_size, 0, h_def};
Point(8) = {x_size,              y_size, 0, h_def};

// bot left

Point(9)  = {0,        -b_y_size, 0, h_def};
Point(10) = {b_x_size, -b_y_size, 0, h_def};
Point(11) = {b_x_size,         0, 0, h_def};
Point(12) = {0,                0, 0, h_def};

// bot right
Point(13) = {x_size,            -b_y_size, 0, h_def};
Point(14) = {x_size + b_x_size, -b_y_size, 0, h_def};
Point(15) = {x_size + b_x_size,         0, 0, h_def};
Point(16) = {x_size,                    0, 0, h_def};


// top left
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// top right
Line(5)  = {5,   6};
Line(6)  = {6,   7};
Line(7)  = {7,   8};
Line(8)  = {8,   5};

// bot left
Line(9)   = {9,    10};
Line(10)  = {10,   11};
Line(11)  = {11,   12};
Line(12)  = {12,   9};

// bot right
Line(13)  = {13,   14};
Line(14)  = {14,   15};
Line(15)  = {15,   16};
Line(16)  = {16,   13};


Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Line Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3};

Line Loop(4) = {13, 14, 15, 16};
Plane Surface(4) = {4};

Physical Point("topleftblock_bottom_left")  = {1};
Physical Point("topleftblock_bottom_right") = {2};

Physical Line("topleftblock_top") = {3};
Physical Line("topleftblock_bottom") = {1};
Physical Line("topleftblock_left") = {4};
Physical Line("topleftblock_right") = {2};

Physical Line("toprightblock_bottom") = {5};
Physical Line("toprightblock_top") = {7};
Physical Line("toprightblock_left") = {8};
Physical Line("toprightblock_right") = {6};

Physical Line("botleftblock_bottom") = {9};
Physical Line("botleftblock_top") = {11};
Physical Line("botleftblock_left") = {12};
Physical Line("botleftblock_right") = {10};

Physical Line("botrightblock_bottom") = {13};
Physical Line("botrightblock_top") = {15};
Physical Line("botrightblock_left") = {16};
Physical Line("botrightblock_right") = {14};

Physical Surface("top_L")  = {1};
Physical Surface("top_R") = {2};
Physical Surface("bot_L")  = {3};
Physical Surface("bot_R") = {4};

Transfinite Surface "*";
Recombine Surface "*";

Mesh.SecondOrderIncomplete = 1;

"""


def main(argv):


    if len(argv)!=3:
        raise RuntimeError("Usage: ./generate_2x2_block3d.py Nx Lx,Ly")
    else:
        Nx = int(argv[1])
        Lx,Ly = [float(i) for i in argv[2].split(',')]

    info={'Nx':Nx,
          'Lx':Lx,
          'Ly':Ly}

    fname = "2x2block_2d_Ll{Lx}r{Lx}_Ht{Ly}b{Ly}_Theta90_msh{Nx}_quad4.geo".format(**info)
    
    with open(fname, "w") as f:
        print(fname)
        print(mesh_file.format(**info)+mesh_file_block,file=f)
        
if __name__ == "__main__":
    import sys
    main(sys.argv)
