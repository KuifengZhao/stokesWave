# stokesWave
List of files 

1. StokesTheoryCombined.mlx 
This file is used to derive solutions for Stokes waves following the manuscript 'On Stokes Wave Solutions'. To be used with Matlab Live Editor for best view of equations. 

2. toUseStokesTheory.m 
This file demonstrate how to use the functions prepared for Stokes wave.

3. StokesDispSolver.m 
This function is a solver for nonlinear dispersion relationship. Different cases were provided to find wave length, for example (1) given a, h, T, find H and k; (2) given H, h, T, find a and k; (3) given aw, h, T, find Hw and k, etc.

4. StokesEta.m 
With the wave amplitude and wave number solved using Function 3, calculate the free surface using Stokes theory.

5. StokesU.m
With the wave amplitude and wave number solved using Function 3, calculate the velocity using Stokes theory.

6. normalizeFreeSurface.mlx
This file is to be opened using Matlab Live Editor. This file shows the process of converting the free surface equations in S and H, and Fenton to the format in the present manuscript for comparison. 

7. normalizePotentialFunction.mlx
This file shows the process of converting the potential functions in S and H, and Fenton to the format in the present manuscript for comparison. 

8. FentonU.m
This is a function to find velocity using Fenton's equations.

9. FentonEta.m
This is a function to find free surface profiles using Fenton's equations.

10. FentonDispSolver.m
This function is the solver to find wave number using the nonlinear dispersion relationship in Fenton's paper. 

Others
hyperIden.m
This is a function to change format of hyperbolic function to \cosh () using the hyperbolic identities.

trigIden.m
This is a function to change format of trigonometric function to \sin () or \cos () using the trigonometric identities.

