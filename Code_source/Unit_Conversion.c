/*----------Conversions of some units--------

References:
Massman, W.J. (1998) A review of the molecular diffusivities of H2O, CO2, CH4, CO,
O3, SO2, NH3, N2O, NO, and NO2 in air, O2 and N2 near STP. Atmospheric Environment,
32, 1111-1127.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "variables_parametres.h"

extern double calcul_pph2o(double);

void Unit_conversion(void)
{
	LAI_yellow = LAI_total - LAI_green;
	
	exa = calcul_pph2o(Ta);
	ea = RHa * calcul_pph2o(Ta) / 100.0;
	PAR = max(0.001 , Rg * 2.0) ;
	VPDa = exa - ea;

	if (uref < 0.6)
		uref = 0.6;

	// Diffusivities are calculated as a function of air temperature as in Massman (1998) 
	diffusivity_H2O = diffusivity_H2O_0C * 1.0 * pow((Ta + T0C) / T0C, 1.81);
	diffusivity_O3 = diffusivity_O3_0C * 1.0 * pow((Ta + T0C) / T0C, 1.81);
	diffusivity_NH3 = diffusivity_NH3_0C * 1.0 * pow((Ta + T0C) / T0C, 1.81);

	Sc = Viscosity_cin / diffusivity_H2O;	
	Sc_NH3 = Viscosity_cin / diffusivity_NH3;
	Sc_O3 = Viscosity_cin / diffusivity_O3;
	Sc_Pest = Viscosity_cin / diffusivity_Pest;
	
}