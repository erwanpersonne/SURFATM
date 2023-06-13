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

/*----------Initialization of variable and parameters--------

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



void initialization(void)
{


	SWC_SoilTotal = Soil_Depth *rho_soil * teta_ini;
	SWC_mulch = 0.0;
	beta_soil = 7.0;
	SWC_WetLayer = SWC_SoilTotal;
	teta_WetLayer = teta_ini;
	depth_drysoil = 0.0;												// in meters
	WetLayer_Depth = Soil_Depth;										// in meters
	depth_threshold_wetsoil = Soil_Depth;								// in meters
	SWC_Wiltingpoint = Soil_Depth * rho_soil * teta_wiltingpoint;
	SWC_Dry = Soil_Depth * teta_minthreshold;
	SWC_Fieldcapacity = Soil_Depth * rho_soil *teta_fieldcapacity;	

	Rain_t3 = 0.0;														// in mm
	Rain_t2 = 0.0;														// in mm
	Rain_t1 = 0.0;														// in mm

	evap_t3 = 0.0;														// in mm
	evap_t2 = 0.0;														// in mm
	evap_t1 = 0.0;														// in mm
	Cumul_evap_t3 = 0.0;


	// parameters
	critmaxStab = 0.01;													// threshold for iteration convergence of the energy balance focused on T° - [in °C]  												
	incrmaxStab = 15;													// maximal number of itertion for stability temperature  
	critmaxRn = 0.05;													// threshold for iteration convergence of the radiation budget focused on T° - [in °C]  												
	incrmaxRn = 100;													// maximal number of itertion for radiation budget   

	critmaxPprime = 0.01;												// threshold for iteration convergence of PPrime calculation focused on T° - [in °C]
	incrmaxPprime = 10;													// maximal number of itertion for Pprime calculation 

	Mulch_threshold = 0.5;												// ratio of the "mulch filling" before destruction during rain event - ranging between 0 and 1 (i.e., 0 and 100%)		 


	// initiatilization of the pesticid quantity on and in the plant and soil
	Pest_quantity_leaf = 0.0;
	Pest_quantity_leaf_adsorb = 0.0;
	Pest_quantity_leaf_non_adsorb = 0.0;
	Pest_quantity_leaf_Tissue = 0.0;
	Pest_quantity_soil = 0.0;
}
