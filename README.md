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
