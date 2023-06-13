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

/*----------Some useful functions for the whole program--------

*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "variables_parametres.h"
#include <time.h>
#include <string.h>


void filecopy(void)
{
	/*----------Copy of the results file----------------*/
	FILE* source_file = NULL;
	source_file = fopen("../output/resultats.csv", "rb");

	char filename[999];
	sprintf(filename, "../output/%s_resultats.csv", buffer_time);
	FILE* destination_file = NULL;
	destination_file = fopen(filename, "wb");

	char buffer_copy[55512];
	int Nbread;
	while ((Nbread = fread(buffer_copy, 1, 55512, source_file)) != 0)
		fwrite(buffer_copy, 1, Nbread, destination_file);

	fclose(source_file);
	fclose(destination_file);

	/*----------Copy of the data file----------------*/
	source_file = fopen("../data/data.txt", "rb");

	sprintf(filename, "../output/%s_data.txt", buffer_time);
	destination_file = fopen(filename, "wb");

	while ((Nbread = fread(buffer_copy, 1, 55512, source_file)) != 0)
		fwrite(buffer_copy, 1, Nbread, destination_file);

	fclose(source_file);
	fclose(destination_file);

	/*----------Copy of the global parameters file----------------*/
	source_file = fopen("../data/global_parametres.txt", "rb");

	sprintf(filename, "../output/%s_global_parametres.txt", buffer_time);
	destination_file = fopen(filename, "wb");

	while ((Nbread = fread(buffer_copy, 1, 55512, source_file)) != 0)
		fwrite(buffer_copy, 1, Nbread, destination_file);

	fclose(source_file);
	fclose(destination_file);

	/*----------Copy of the NH3 parametres file----------------*/
	source_file = fopen("../data/NH3_parametres.txt", "rb");

	sprintf(filename, "../output/%s_NH3_parametres.txt", buffer_time);
	destination_file = fopen(filename, "wb");

	while ((Nbread = fread(buffer_copy, 1, 55512, source_file)) != 0)
		fwrite(buffer_copy, 1, Nbread, destination_file);

	fclose(source_file);
	fclose(destination_file);

	/*----------Copy of the O3 parametres file----------------*/
	source_file = fopen("../data/O3_parametres.txt", "rb");

	sprintf(filename, "../output/%s_O3_parametres.txt", buffer_time);
	destination_file = fopen(filename, "wb");

	while ((Nbread = fread(buffer_copy, 1, 55512, source_file)) != 0)
		fwrite(buffer_copy, 1, Nbread, destination_file);

	fclose(source_file);
	fclose(destination_file);

	/*----------Copy of the Pesticides parametres file----------------*/
	source_file = fopen("../data/Pesticides_parametres.txt", "rb");

	sprintf(filename, "../output/%s_Pesticides_parametres.txt", buffer_time);
	destination_file = fopen(filename, "wb");

	while ((Nbread = fread(buffer_copy, 1, 55512, source_file)) != 0)
		fwrite(buffer_copy, 1, Nbread, destination_file);

	fclose(source_file);
	fclose(destination_file);

	/*----------Copy of the simulation options file----------------*/
	source_file = fopen("../data/simulation_options.txt", "rb");

	sprintf(filename, "../output/%s_simulation_options.txt", buffer_time);
	destination_file = fopen(filename, "wb");

	while ((Nbread = fread(buffer_copy, 1, 55512, source_file)) != 0)
		fwrite(buffer_copy, 1, Nbread, destination_file);

	fclose(source_file);
	fclose(destination_file);

}


double calcul_pph2o(double Thar)
{
double pph2o ;

/* Formula of 'TETENS' */

 pph2o = 610.78 * exp((17.269 * Thar)/(Thar + 237.3)) ;

 return(pph2o) ;
} 



/* calculation of PTr if necessary from Th and Ta */
double calcul_PTr(double Th, double Ta)
{
double PTr ;

 PTr = calcul_pph2o(Th) - PsychroCst * (Ta - Th) ;

 return(PTr) ;
}

/* calculation of P' */
double calcul_Pprime(double Ta)
{
float Pprime ;

 Pprime = calcul_pph2o(Ta + 0.5) - calcul_pph2o(Ta - 0.5) ;

 return(Pprime) ;
}


/* calculation of the latente heat of vaporization as a function of Ta (in °C) */
double calcul_latente(double Ta)
{
double latente ;

 latente = 2500841.0 - 2358.6 * Ta ;

 return(latente) ;
}

/* --------------------------------------------------------------------------------- */


void Calcul_Rn(double tested_temperature)
{
	Rn_calculated = (1 - albedo_canopy)*Rg + emissivity_canopy * (Ratm - SIGMA * pow(tested_temperature + T0C, 4.0));

if (Flag_Rn == 0)	
	Rn = Rn_calculated;


if (Flag_Rn == 1)
	Rn = Rn_measured;

}

void Calcul_deltaprime123(void)
{
	// Calculation of deltaPprimes accounting for calculated temperatures

	Pprime1 = (ex_leaf - ex_canopy) / (T_leaf - T_canopy);
	Pprime2 = (ex_soilwetdry - ex_canopy) / (T_soilwetdry - T_canopy);
	Pprime3 = (ex_canopy - exa) / (T_canopy - Ta);

	if ((zh - z0_soil) <= 0.03 && LAI_total < 0.5)
		Pprime1 = Pprime2 = Pprime3 = calcul_pph2o(Ta + 0.5) - calcul_pph2o(Ta - 0.5);
}


