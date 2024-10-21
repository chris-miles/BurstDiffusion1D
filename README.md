# BurstDiffusion1D

 Code for " Incorporating spatial diffusion into bursty stochastic models of transcription."

Requires MATLAB PDE Toolbox. Tested on MATLAB 2024a. 

Any questions or concerns can be directed to chris.miles@uci.edu.

## Usage

1. `fig_sweep_xsource.m` generates the histograms Fig. 2a with $\kappa=\infty$. 
2. `fig_sweep_kappainf.m` plots various quantities for different parameters and $\kappa=\infty$, shown in Figs. 2b, 3, and 4b.
3. `fig_sweep_xsource.m` generates histograms for different sources, Fig 4a.
4. `fig_sweep_finitekappa.m` plots various quantities for finite kappa,  Fig. 5.
5. `fig_comparefano.m` plots various Fano factors for different parameters, Fig. 6 
6. `fig_diffhetgeom.m` generates a random domain, solves the PDE (using MATLAB's toolbox), and generates Fig. 7.
7. `fig_infer_robin.m` performs the profile likelihood analysis and generates Fig. 8.

## Credits

`cmap` colors from https://github.com/tsipkens/cmap

`Poissbeta` Gauss-Jacob integration of Poisson-Beta PDF from https://github.com/cellfate/4DNucleomeEquation

`kldiv` from https://www.mathworks.com/matlabcentral/fileexchange/13089-kldiv

