#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"

using namespace alglib;


int main(int argc, char **argv)
{
    //
    // We use trilinear spline to interpolate f(x,y,z)=x+xy+z sampled 
    // at (x,y,z) from [0.0, 1.0] X [0.0, 1.0] X [0.0, 1.0].
    //
    // We store x, y and z-values at local arrays with same names.
    // Function values are stored in the array F as follows:
    //     f[0]     (x,y,z) = (0,0,0)
    //     f[1]     (x,y,z) = (1,0,0)
    //     f[2]     (x,y,z) = (0,1,0)
    //     f[3]     (x,y,z) = (1,1,0)
    //     f[4]     (x,y,z) = (0,0,1)
    //     f[5]     (x,y,z) = (1,0,1)
    //     f[6]     (x,y,z) = (0,1,1)
    //     f[7]     (x,y,z) = (1,1,1)
    //
    real_1d_array x = "[0.0, 1.0]";
    real_1d_array y = "[0.0, 1.0]";
    real_1d_array z = "[0.0, 1.0]";
    real_1d_array f = "[0,1,0,2,1,2,1,3]";
    double vx = 0.50;
    double vy = 0.50;
    double vz = 0.50;
    double v;
    spline3dinterpolant s;

    // build spline
    spline3dbuildtrilinearv(x, 2, y, 2, z, 2, f, 1, s);

    // calculate S(0.5,0.5,0.5)
    v = spline3dcalc(s, vx, vy, vz);
    printf("%.4f\n", double(v)); // EXPECTED: 1.2500
    return 0;
}
