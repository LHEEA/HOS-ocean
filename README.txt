HOS-ocean

=========

High-Order Spectral method for oceanic simulations

This readme file describes the different cases that might be computed with the periodic HOS model
and gives instructions to set the numerical parameters n1 to p2

Setting the value of integers n1, n2, M, p1 and p2 in variables_3D.f90

For a 2D simulation, 
	Compile with n2=1 AND p2=1 to adjust the memory allocation to minimum
	
	If partial dealiasing is used, compile with p1 set to maximal required value (total dealiasing is obtained
		with p1=M but it can be reduced if p1 is further set to a value below M

For a 3D simulation, 
	Compile with n2\=1 AND p2 set to required value
	
	If partial dealiasing is used in x-direction, compile with p1 set to maximal required value (total dealiasing is obtained
		with p1=M but it can be reduced if p1 is further set to a value below M
	
	If partial dealiasing is used in y-direction, compile with p2 set to maximal required value (total dealiasing is obtained
		with p2=M but it can be reduced if p2 is further set to a value below M

Setting the value of integer i_case.

i_case = 1
	starts from rest
i_case = 2 and 21
	starts with a natural mode, either progressive (2) or stationary (21). 
	The number of the mode, its amplitude and phase have to be chosen in the module initial_condition.f90.
	The x-component may be positive or negative for the progressive case.
	Should be positive for the stationary case.
	The y-component may be positive or negative for both cases.
	In the progressive case, the potential on the free surface is set accordingly to linear theory.
	In the stationary case, it is set to zero (no velocity at t=0).
i_case = 3 and 31 and 32
	Irregular directional sea-state (linear) from spectrum specified with 
		input file parameters (3)
		WAVEWATCH III® spectrum file (31)
		previous HOS-ocean simulation (32)
i_case = 8...
 	Rienecker and Fenton test
	this means that xlen and ylen are interpreted as integer and represent the number of wavelength 
	in x- and y- direction respectively
	i_case = 81
		steepness 0.1
	i_case = 82
		steepness 0.2
	i_case = 83
		steepness 0.3
	i_case = 84
		steepness 0.4
	i_case = 89
		steepness 0.09

Further details about input file, output of the code... may be find at the Wiki page of HOS-ocean project: https://github.com/LHEEA/HOS-ocean/wiki
