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
mpirun -n ${N} xcompact3d $nx $ny $nz $prow $pcol $rt $test
```
where `$nx`, `$ny`, `$nz` specify the grid size, `$prow <= $pcol`, `$prow * $pcol = $N` specifies
parallel decomposition, `$rt` specifies how many seconds to run the benchmark for - note the
total time will be greater than this - and `$test` enables or disables computing the error in the
derivatives.
Each variable has a default value:
```
nx = ny = nz = 16
prow = pcol = 0
rt = 5
test = 0
```
setting `$prow=$pcol=0` the code will attempt to determine a "good" parallel decomposition, `$rt`
must be an integer number of seconds and `$test` is an integer with zero interpreted as logical
`false` and any other value as `true`.

For benchmarking `$test` mode should be disabled, its intention is to validate custom
implementations of the compact finite difference scheme solvers.

Note this code is a very stripped down version of Xcompact3d, it is intended for profiling only.
