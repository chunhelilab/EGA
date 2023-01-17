Code for "An improved approach for calculating energy landscape of gene networks from moment equations"
——————————————————————————————————————————————————

Code_multistability includes all the codes of the 2-dimensional genetic circuit momel
Code_limitcycle includes all the codes of the synthetic oscillatory network

In the first folder (Code_multistability):

`WSGA_and_EGA' includes the codes for WSGA and EGA:
	The main programs are `Land.m' and `Landnew.m'. They describe the landscape of WSGA and EGA, respectively.
	`ode_cov.m' and `ode_covnew.m' are corresponding functions of moments.
	`intnew.m' is a kind of high performance numerical integration approach in 2-D.
	`fUn.m' is the characteristic function of EGA.

`LE_simulation' includes the codes for LE simulation:
	The main programs is `pareuler.m', which describes the landscape of LE simulation
	`force.m' is the corresponding function of Explict Euler method.

In the second folder (Code_limitcycle):

`WSGA_and_EGA' includes the codes for WSGA and EGA:
	The main programs are `Land_limit3.m' and `Landnew_limit3.m'. They describe the landscape of WSGA and EGA, respectively.
	`ode_cov_limit3.m' and `ode_covnew_limit3.m' are corresponding functions of moments.
	`momentdata.mat' is a separate data from LE to be the initial value of moment equations in WSGA and EGA.
	`int3new.m' is a kind of high performance numerical integration approach in 3-D.
	`fUn3.m' is the characteristic function of EGA.

`LE_simulation' includes the codes for LE simulation:
	The main programs is `pareulernew.m', which describes the landscape of LE simulation
	`forcenew.m' is the corresponding function of Explict Euler method.
