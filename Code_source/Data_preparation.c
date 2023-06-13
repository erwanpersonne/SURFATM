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

/*----------Program for preparation of the input data--------

Some checks are done to avoid errors in calculation
Calculation of vegetation and canopy albedos are also performed

References:
Montes, C., Lhomme, J.P., Demarty, J., Prévot, L., Jacob, F. (2014) A three-source
	SVAT modeling of evaporation: Application to the seasonal dynamics of a grassed 
	vineyard. Agricultural and Forest Meteorology, 191, 64-80.

	ESSAI 1
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "variables_parametres.h"


void Data_preparation(void)
{
	if (LAI_total <= 0.0)
		LAI_total = 0.001;
	
	if (LAI_green <= 0.0)
		LAI_green = 0.001;

	if (LAI_yellow <= 0.0)
		LAI_yellow = 0.001;

	if (displacement_height >= zh)
		displacement_height = 0.75 * zh;

	if (z0_canopy <= 0.0 || z0_canopy <= z0_soil)
		z0_canopy = 0.1 * zh;

	if (displacement_height <= 0.0)
		displacement_height = 0.75 * zh;
		
	if (LAI_total < 0.5 && zh < 0.03)
	{
		displacement_height = 0.0;
		z0_canopy = z0_soil;
	}
	
	if (zh <= 0.0 || zh <= (1.5* z0_soil) || displacement_height >= zh || z0_canopy >= displacement_height || z0_canopy >= zh)
	{
		z0_canopy = z0_soil + 0.001;
		zh = z0_canopy + 0.03;
		displacement_height = 0.0;
	}
	
	canopy_cover = 1.0 - exp(-attenuationray*LAI_total); // Coefficient "attenuationray" is identical here and for Rn partitionning. Cf Montes et al. (2014)

	albedo_veg = (LAI_green / LAI_total) * albedo_veg_green + (LAI_yellow / LAI_total) * albedo_veg_yellow ;
	albedo_canopy = canopy_cover*albedo_veg + (1 - canopy_cover)*albedo_soil;

}
