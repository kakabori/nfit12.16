Make and run test_ModelCalculator_class. It will create four files named,

test_CCD_struct_factor.dat
test_mosaic_struct_factor.dat
test_rotated_struct_factor.dat
test_struct_factor.dat

Check that calculated values in each test_*.dat and its corresponding 
keep_*.dat are not different. If they significantly differ, then a newly 
implemented change does something strange.

Note that keep_*.dat files were created on 8/20/2013.

test_CCD_struct_factor.dat
This file contains structure factor in CCD space. It is calculated after
mosaic convolution, rotation integration, and beam convolution.

test_mosaic_struct_factor.dat
This file contains structure factor after mosaic convoluted, but before
rotation integration (i.e. qy integration).

test_rotated_struct_factor.dat
This file contains structure factor after rotation integration, which is 
before beam convolution but after mosaic convolution.

test_struct_factor.dat
This file contains pure structure factor without any extra convolution or
integration. This is the most fundamental and must be unchanged unless
one modifies the height difference correlation function.
