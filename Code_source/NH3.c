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

The IDDN number assigned to your deposit is: IDDN.FR.001.180010.000.S.P.2022.000.30100.

The electronic certificate of your deposit is available at the following address:
https://secure2.iddn.org/app.server/certificate/?sn=2022180010000&key=87195b4227e4c9142f5e7ff3d81c52388a1492c19ea63f08929474f4a7d4ee1d&lang=fr

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

/*----------Program for NH3 flux calculations--------

References:
Personne, E., Loubet, B., Herrmann, B., Mattson, M., Schjoerring, J.K., Nemitz, E.,
Sutton, M.A., Cellier, P. (2009) SURFATM-NH3: a model combining the surface energy
balance and bi-directional exchanges of ammonia applied at the field scale.
Biogeosciences,	6, 1371-1388.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "variables_parametres.h"
#include "fonctions.h"

void NH3_resistance_calculations(void);
void NH3_Flux_calculations(void);
void NH3_compensation_caluculation(void);

void Flux_NH3()
{
	NH3_compensation_caluculation();
	NH3_resistance_calculations();
	NH3_Flux_calculations();
}

void NH3_resistance_calculations()
{
	// Stomatal conductance of green leaves
	Rstom_NH3 = diffusivity_H2O / diffusivity_NH3 * Rstom_H2O;
	gstom_NH3 = 1.0 / Rstom_NH3;

	// Stomatal conductance of yellow leaves
	Rstom_inactive_NH3 = diffusivity_H2O / diffusivity_NH3 * Rstom_inactive_H2O;
	gstom_inactive_NH3 = 1.0 / Rstom_inactive_NH3;

	// Cuticular conductance
	Rcut_NH3 = Rcutmin_NH3 * exp((100.0 - RH_canopy) / kcut_NH3);
	gcut_NH3 = 1.0 / Rcut_NH3;

	// Soil conductance
	Rsoil_NH3 = R_litter_NH3 + r_soil_surface_H2O * diffusivity_H2O / diffusivity_NH3;
	gsoil_NH3 = 1.0 / Rsoil_NH3;


	// Soil and canopy boundary layer conductances
	Rb_leaf_NH3 = diffusivity_H2O / diffusivity_NH3 * Rb_leaf;
	if (Rb_leaf_NH3 > 5000.0)
		Rb_leaf_NH3 = 5000.0;
	Rb_soil_NH3 = 2.0 / (karman * u_star_ground) * pow((Sc_NH3 / Prandt), 2 / 3);
	if (Rb_soil_NH3 > 5000.0)
		Rb_soil_NH3 = 5000.0;
	
	

	// Controle of the values of the resistances/conductances for very low canopies
	if ((zh - z0_soil) <= 0.03)
	{
		gstom_NH3 = 0.000001;
		Rstom_NH3 = 1000000.0;
		gstom_inactive_NH3 = 0.000001;
		Rstom_inactive_NH3 = 1000000.0;
		Rcut_NH3 = 1000000.0;
		gcut_NH3 = 0.000001;
		Rb_leaf_NH3 = 1000000.0;
	}
}

void NH3_compensation_caluculation()
{
	double KhilNH3, KhiSNH3; // local variable for the calculation of the NH3 concentration in the soil or "leaf" layer


	KhilNH3 = KHenry25C_NH3 * KAcidB25C_NH3 * exp((enthalpAcidB_NH3 + enthalpVap_NH3) / R*(1.0 / 298.15 - (1.0 / (T_leaf+T0C)))) * GammaLeaf_NH3; //compensation point in the substomatal cavity [mol /m3]
	KhilNH3 = KhilNH3 * 17.e-3 *1.e9;																												// Conversion in µg/m3
	
	KhiSNH3 = KHenry25C_NH3 * KAcidB25C_NH3 * exp((enthalpAcidB_NH3 + enthalpVap_NH3) / R*(1.0 / 298.15 - (1.0 / (T_soilwetdry +T0C)))) * GammaSoil_NH3; // compensation point in the wet soil layer [mol /m3]
	KhiSNH3 = KhiSNH3 * 17.e-3 *1.e9;                                                                                                 // Conversion in µg/m3
	
	NH3_i = KhilNH3;
	NH3_soil = KhiSNH3;
}

void NH3_Flux_calculations()
{


	NH3_canopy = (((1./Rb_leaf_NH3) + gstom_NH3 + gcut_NH3)*((1.0/Ra) * concentration_NH3 + (1. / (Rac + Rb_soil_NH3 + Rsoil_NH3)) * NH3_soil) + (1. / Rb_leaf_NH3)*(gstom_NH3) * NH3_i) /
																	(((1. / Rb_leaf_NH3) + gstom_NH3 + gcut_NH3)*((1./Ra) + (1. / (Rac + Rb_soil_NH3 + Rsoil_NH3)) + (1. / Rb_leaf_NH3)) - (1. / Rb_leaf_NH3)*(1. / Rb_leaf_NH3));

	NH3_leaf = ((1./Rb_leaf_NH3) * NH3_canopy + gstom_NH3 * NH3_i) / ((1. / Rb_leaf_NH3) + gstom_NH3 + gcut_NH3); // [NH3] just above the leaf surface (below Rb)

	FNH3_soil = -(1. / (Rac + Rb_soil_NH3 + Rsoil_NH3)) * (NH3_canopy - NH3_soil);            
	FNH3_canopy = -(1./Rb_leaf_NH3) * (NH3_canopy - NH3_leaf);                    
	FNH3_cut = -gcut_NH3*NH3_leaf;                             
	FNH3_stom = -gstom_NH3 * (NH3_leaf - NH3_i);

	FNH3_tot = -(1./Ra) * (concentration_NH3 - NH3_canopy);   

	// Controle of the values of the fluxes/concentrations for very low canopies
	if ((zh - z0_soil) <= 0.03 )
	{
		NH3_canopy = 0.000001;
		FNH3_soil = -(1. / (Ra + Rb_soil_NH3 + Rsoil_NH3)) * (concentration_NH3 - NH3_soil);
		FNH3_tot = FNH3_soil;
		FNH3_canopy = -0.000001;
		NH3_leaf = 0.000001;
		FNH3_cut = -0.000001;
		FNH3_stom = -0.000001;
	}

}
