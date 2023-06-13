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

The IDDN number assigned to your deposit is: IDDN.FR.001.180010.000.S.P.2022.000.30100.

The electronic certificate of your deposit is available at the following address:
https://secure2.iddn.org/app.server/certificate/?sn=2022180010000&key=87195b4227e4c9142f5e7ff3d81c52388a1492c19ea63f08929474f4a7d4ee1d&lang=fr

*/
/*
@author Erwan Personne <erwan.personne@agroparistech.fr>
@technical support : support.surfatm@agroparistech.fr
*/

/*----------Program for energy balance--------

Here are the calculations of each terms of the energy balance (H, LE, and G) and their
partitionning between soil and vegation. Temperature and humidity at the vegetation and
soil surfaces are also calculated


References:
Personne, E., Loubet, B., Herrmann, B., Mattson, M., Schjoerring, J.K., Nemitz, E.,
Sutton, M.A., Cellier, P. (2009) SURFATM-NH3: a model combining the surface energy
balance and bi-directional exchanges of ammonia applied at the field scale.
Biogeosciences,	6, 1371-1388.

Choudhury, B.J., Monteith, J.L., 1988. A four-layer model for the heat budget of homogeneous land surfaces. 
Quarterly Journal of the Royal Meteorological Society 114, 373–398. https://doi.org/10.1002/qj.49711448006

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "variables_parametres.h"

extern double
calcul_PTr(double, double)
, calcul_pph2o(double)
, calcul_Pprime(double)
, calcul_latente(double)
;

extern void Calcul_deltaprime123(void);


void Coeff_energy_resolution(void);
void Calcul_temperature(void);


void Energy_balance(void)
{

	Coeff_energy_resolution();

		LEv = (Pprime1 * Rn_veg + (rho_air*Cp_air * VPD_canopy / Rb_leaf)) / (Pprime1 + Gamma1);
		Hv = (Rn_veg   * Gamma1 - (rho_air*Cp_air * VPD_canopy / Rb_leaf)) / (Pprime1 + Gamma1);

		LEs = (Pprime2*etha*Rn_soil + (Pprime2*rho_air*Cp_air*(T_soilref - T_canopy) / r_thermal_wetlayer) + rho_air*Cp_air*coeffnu* VPD_canopy / (Rac + Rb_soil)) / (Pprime2 + Gamma2);
		Hs = Rn_soil
			+ ((rho_air*Cp_air*(T_soilref - T_canopy) - ((Rac + Rb_soil)*Rn_soil)) / (coeffnu*r_thermal_wetlayer / etha))
			- ((Pprime2*Rn_soil*etha + (Pprime2*rho_air*Cp_air*(T_soilref - T_canopy) / r_thermal_wetlayer) + (rho_air*Cp_air*coeffnu*VPD_canopy / (Rac + Rb_soil))) / ((coeffnu / etha) * (Pprime2 + Gamma2)));

		// Taking into account the energy balance, G is calculated
		G = Rn_soil - Hs - LEs;
		LE = LEv + LEs;
		H = Hv + Hs;

		Calcul_temperature();

}

void Coeff_energy_resolution()
{
		double Q1, Q2, Q3, Q4, Q5;
		double qu1, qu2, qu3, qu4, qu5;
		double GrandQ;

		// Choudhury-Monteith, annexe
		etha = 1.0 / (1.0 + r_thermal_drylayer / (Rac + Rb_soil));
		coeffnu = (1 + (r_thermal_drylayer + (Rac + Rb_soil)) / r_thermal_wetlayer) * etha;
		Gamma1 = PsychroCst * (1 + Rstom_tot_H2O / Rb_leaf);
		Gamma2 = PsychroCst * (1 + r_soil_surface_H2O / (Rac + Rb_soil)) * coeffnu;

		Q5 = Pprime2*rho_air*Cp_air / (r_thermal_wetlayer*(Pprime2 + Gamma2));
		Q4 = (Pprime1 * Rn_veg / (Pprime1 + Gamma1)) + ((Pprime2*Rn_soil*etha) / (Pprime2 + Gamma2));
		Q3 = rho_air*Cp_air*coeffnu / ((Rac + Rb_soil)*(Pprime2 + Gamma2)) + (rho_air*Cp_air / (Rb_leaf*(Pprime1 + Gamma1)));
		Q2 = Q3 * Pprime3 - Q5;
		Q1 = Q4 + (Q5 * (T_soilref - Ta)) + (Q3 * VPDa);

		qu5 = etha * rho_air*Cp_air / (coeffnu *r_thermal_wetlayer) - (etha * Pprime2 * rho_air*Cp_air / (coeffnu* r_thermal_wetlayer*(Pprime2 + Gamma2)));
		qu4 = Rn_soil + Rn_veg * Gamma1 / (Pprime1 + Gamma1) - etha * (Rac + Rb_soil) * Rn_soil / (coeffnu * r_thermal_wetlayer) - (etha*etha*Pprime2*Rn_soil / (coeffnu*(Pprime2 + Gamma2)));
		qu3 = -rho_air*Cp_air / (Rb_leaf * (Pprime1 + Gamma1)) - (rho_air*Cp_air*etha / ((Rac + Rb_soil) * (Pprime2 + Gamma2)));
		qu2 = qu3 * Pprime3 - qu5;
		qu1 = qu4 + qu5*(T_soilref - Ta) + qu3*VPDa;

		GrandQ = (rho_air*Cp_air / (PsychroCst*Ra) + Q3) * (rho_air*Cp_air / Ra - qu2) + Q2 * qu3;

		e_canopy = ea + ((Q1  * (rho_air*Cp_air / Ra - qu2)) + (Q2 * qu1)) / GrandQ;
		T_canopy = Ta + ((qu1 * (rho_air*Cp_air / (PsychroCst*Ra) + Q3)) - (Q1 * qu3)) / GrandQ;


		ex_canopy = calcul_pph2o(T_canopy);
		VPD_canopy = ex_canopy - e_canopy;
		RH_canopy = min(e_canopy / ex_canopy, 1.0)*100.0;

}


void Calcul_temperature()
{
	// Knowing the fluxes, temperature and humidity at the leaf and soil surfaces can
	// be calculated

	T_leaf = Rb_leaf * Hv / (rho_air*Cp_air) + T_canopy;
	T_soil = (Rac + Rb_soil) * Hs / (rho_air*Cp_air) + T_canopy;
	T_soilwetdry = r_thermal_wetlayer * G / (rho_air*Cp_air) + T_soilref;

	if ((zh - z0_soil) <= 0.03 && LAI_total < 0.5)
	{
		T_leaf = Ta;
	}

	ex_leaf = calcul_pph2o(T_leaf);
	ex_soil = calcul_pph2o(T_soil);
	ex_soilwetdry = calcul_pph2o(T_soilwetdry);


	e_soil = (((e_canopy*18.0) / (R*(T_canopy + T0C)*1000.0)) + ((LEs*(Rac + Rb_soil)) / calcul_latente(T_canopy)))*((R*(T_soil + T0C)*1000.0) / 18.0);
	e_leaf = (((e_canopy*18.0) / (R*(T_canopy + T0C)*1000.0)) + ((LEv*Rb_leaf) / calcul_latente(T_canopy)))*((R*(T_canopy + T0C)*1000.0) / 18.0);

	RH_soil = min(e_soil / ex_soil,1.0) * 100.0;
	RH_leaf = min(e_leaf / ex_leaf,1.0) * 100.0;

	if ((zh - z0_soil) <= 0.03 && LAI_total < 0.5)
	{
		e_leaf = ea;
		ex_leaf = exa;
		RH_leaf = RHa;

		e_canopy = ea;
		ex_canopy = exa;
		RH_canopy = RHa;
	}

}

