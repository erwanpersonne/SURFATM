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