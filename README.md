This version implements 2D mosaic integration.
For the theory behind this implementation, see
Kiyotaka Akabori, Ph.D. thesis, Appendix.

The goal of this version is:
1) implement 2D mosaic integration using 2D integrator
2) implement multithread capability within ModelCalculator class

Item (1) will be accomplished using GLS library.
Item (2) will be a bit complicated b/c the program
needs to calculate a big map twice, once for a pure
structure factor and another for a mosaic integration.

What you need to compile:
- tcl8.4-dev
- tk8.4-dev
- blt2.4z-dev

blt2.4z-dev can be downloaded from: http://sourceforge.net/projects/blt/

As of 2014, the blt package in the Ubuntu repository requires tcl/tk8.6, 
which is incompatible
with this version of NFIT. To avoid a blt version conflict, it is probably 
best to uninstall the blt package from the Ubuntu repository (use synaptic
or apt-get), and manually install the blt package downloaded from the above link.

To compile a downloaded blt package, you might need to add paths to 
tclConfig.sh and tkConfig.sh. Do that by:

./configure --with-tcl=/usr/lib/tcl8.4 --with-tk=/usr/lib/tk8.4

The exact paths might differ. You need to find out for yourself. 

To have both versions of blt (one with 8.4 and another with 8.6), 
you might want to install a dowloaded blt to a user-specified directory. 
Then, add an appropriate search path in the makefile.
