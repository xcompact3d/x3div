# x3div

The tri-diagonal solvers used for computing the derivatives in Xcompact3d do not currently vectorise
well, particularly in the x-direction, this project extracts just these subroutines in order to
investigate possible improvements.

Note this code is a very stripped down version of Xcompact3d, it is intended for profiling only.
