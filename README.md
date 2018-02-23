# phasefield
Code for the Allen-Cahn benchmark problem here:

https://pages.nist.gov/pfhub/benchmarks/benchmark7.ipynb/

I first tried the Forward-Euler finite differences scheme (with second order in space), but measuring convergence in time was tricky because the spatial error always dominates anywhere the scheme is stable. So, I switched to a Forward-Backward scheme with implicit linear part and explicit linear part. The scheme is stable enough to allow large time steps to meaninfully measure time convergence.

## Compilation
I wrote my code in C++14 so this requires a relative new compiler. Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page) is used for the sparse matrix solver. You can download Eigen from the site and put the Eigen folder on your /usr/local/include. If you wish to put Eigen in a custom folder, add the -I{your custom folder} to the CFLAGS option in the Makefile. Then, you can just call
```bash
make
```

## Part a) convergence test
Here is my plot for the spatial error, as specified in the benchmark problem.

![alt text](plots/FBEspace.png)

The spatial order of accuracy matches the theoretical one (second order).

![alt text](plots/FBEtime.png)

The temporal order of accuracy matches the theoretical one (first order). However, when the time step size is small enough, the spatial error starts dominating.

## Part b)
I'm planning to make additional changes to my code for measuring performance. I'm going to try an adaptive scheme with higher order of accuracy. Will update this part.