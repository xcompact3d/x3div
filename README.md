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

## 2decomp&fft

Like Xcompact3d, x3div builds upon the 2decomp&fft library - rather than copy the code, we now
have a git submodule tracking our fork [2decomp&fft](https://github.com/xcompact3d/2decomp_fft)
which tracks the upstream repo [2decomp&fft-upstream](https://github.com/numericalalgorithmsgroup/2decomp_fft).
This means we can share code more easily (and benefit from others' contributions).

For detailed instructions on using git submodules see [git docs](https://git-scm.com/book/en/v2/Git-Tools-Submodules)
however if you are not planning to work on 2decomp&fft itself the following should suffice:

1) Initial clone of x3div ``git clone https://github.com/xcompact3d/x3div``
2) Initialise and update the 2decomp&fft submodule ``cd x3div && git submodule init && git submodule update``

After which you can continue to build as normal (running ``make`` will first call ``make`` on the 2decomp&fft submodule
and then link the resulting library into x3div).
To ensure you receive the latest changes to 2decomp&fft run ``git submodule update --remote`` periodically, note that you
can also work within the ``decomp2d/`` directory as though it were a standalone git project.

Note that variables are passed down by ``make``, therefore if you have installed ``ffte`` at
``${FFTE_DIR}`` then you can build against this by running

``
make FFT=ffte FFTE_PATH=${FFTE_DIR}
``

where ``FFTE_PATH`` is used by 2decomp&fft to link the appropriate library - see
``decomp2d/src/Makefile.inc`` for different FFT library options.
By default ``FFT=generic`` and no external libraries are required.
