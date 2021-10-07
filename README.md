# x3div

The tri-diagonal solvers used for computing the derivatives in `Xcompact3d` do not currently vectorise
well, particularly in the x-direction.
To test improvements to the `Xcompact3d` code, `x3div` extracts a simplified timestep based on the
Euler equations, including solution of the Poisson equation using 2decomp&fft.

Building the code requires a working `mpi` installation, depending on your compilers build using
```
make CMP=${compiler}
```
where `${compiler}` is one of `intel`, `gcc`, `cray` or `nagfor`.
Running the produced binary will run the default benchmark for 5 seconds on a `16x16x16` grid, the
run can be customised on the commandline as:
```
mpirun -n ${N} xcompact3d $nx $ny $nz $prow $pcol $rt
```
where `$nx`, `$ny`, `$nz` specify the grid size, `$prow <= $pcol`, `$prow * $pcol = $N` specifies
parallel decomposition and `$rt` specifies how many seconds to run the benchmark for - note the
total time will be greater than this.

Note this code is a very stripped down version of Xcompact3d, it is intended for profiling only.
