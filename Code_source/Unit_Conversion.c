/*
This file is part of SurfAtm software
Copyright(c) 2023, – UMR ECOSYS, AgroParisTech INRAe, France

== GNU General Public License Usage ==

SurfAtm is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SurfAtm is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SurfAtm. If not, see <http://www.gnu.org/licenses/>.

== Other Usage ==
Other Usage means a use of SurfAtm that is inconsistent with the GPL
license, and requires a written agreement between You, AgroParisTech and INRAe.
Licensees for Other Usage of SurfAtm may use this file in accordance
with the terms contained in the written agreement between You, AgroParisTech and INRAe.
*/
/*
@author Erwan Personne <erwan.personne@agroparistech.fr>
@technical support : support.surfatm@agroparistech.fr
*/

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
