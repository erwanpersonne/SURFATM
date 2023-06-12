/*----------Partitionning of total Rn in vegetation and soil components--------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "variables_parametres.h"

void Rn_partitionning()
{
	{
		Rn_soil = Rn * exp(-attenuationray*LAI_total);
		Rn_veg = Rn - Rn_soil;
	}

}
