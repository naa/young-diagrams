# young-diagrams
Some code to compute and draw generalized Young diagrams for Lie algebras so(2n+1).
This code was used to produce the illustrations in the paper *A. Nazarov, P. Nikitin, O. Postnova, Limit shape for infinite rank limit of non simply-laced Lie algebras of series so(2n+1), [arXiv:2010.16383](http://arxiv.org/abs/2010.16383)*

The folder `bn-young-diagrams` contains Python 3 code. The module `bnyoungdiagrams.py` can be used to compute the most probable diagram for the given rank `n` of Lie algebra `so(2n+1)` and the tensor power `N` of the spinor representation. It also plots the limit shape for `c=(N+2n-1)/n`. Use it as a standalone program to output the coordinates `{a}` (see the paper) of the rotated diagram and to plot the diagram and the limit shape:
```
python3 bnyoungdiagrams.py -n 20 -N 40 -file plot.pdf
```
The functions from this module can also be imported and used. The most useful are `findmax(n,N)`, which returns the coordinates `{a}` of the most probable diagram, `rho(x,c)` which computes the limit particle density at the point `x` for a given `c`, and two plotting functions: `plot_rotated_diagram(a)` plots the diagram for the list of coordinates `{a}` and `plot_limit_shape(n,N)` which uses `rho` to compute the limit density, then integrates it and plots the limit shape. The package uses `numpy` for computations and `matplotlib` for plots, code is very simple. Jupyter notebook demo is available in the file `Demo.ipynb`. 

Roughly equivalent code in Mathematica is in the file `bn-young-diagrams.nb`. The function to find the most probable diagram is called `findmax`, the function to plot the most probable diagram and the limit shape for the given values of `n` and `N` is `showDiagramAndLimitShape[n,N]`.

The file `bn-steapest-descent.c` is a very basic C program, that computes the most probable diagram. It is much faster then Python and Mathematica implementations and can be used if rank and tensor power are very large. It has the same parameters `n, N` and outputs the coordinates `{a}` to the standard output. 

