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

/*----------Main program of SURFATM--------

It includes:
	- reading of input data and parameter files
	- the algorithmic succession for energy and pollutant exchanges
	- writting the output files

Main references of Surfam (chronological order):
Personne, E., Loubet, B., Herrmann, B., Mattson, M., Schjoerring, J.K., Nemitz, E., 
	Sutton, M.A., Cellier, P. (2009) SURFATM-NH3: a model combining the surface energy
	balance and bi-directional exchanges of ammonia applied at the field scale. 
	Biogeosciences,	6, 1371-1388.
Stella, P., Personne, E., Loubet, B., Lamaud, E., Ceschia, E., Béziat, P., Bonnefond, J.M.,
	Irvine, M., Keravec, P., Mascher, N., Cellier, P. (2011) Predicting and partioning
	ozone fluxes to maize crops from sowing to harvest: the Surfatm-O3 model. Biogeosciences,
	8, 2869-2886.
Stella, P., Personne, E., Lamaud, E., Loubet, B., Trebs, I., Cellier, P. (2013). Assessment
	of the total, stomatal, cuticular, and soil 2 year ozone budgets of an agricultural field 
	with winter wheat and maize crops. Journal of Geophysical Research - Biogeosciences, 118,
	1-13, doi:10.1002/jgrg.20094.
Lichiheb, N., Personne, E., Bedos, C., Barriuso, E., 2014. Adaptation of a resistive model to
	pesticide volatilization from plants at the field scale: Comparison with a dataset. 
	Atmospheric Environment 83, 260–268. https://doi.org/10.1016/j.atmosenv.2013.11.004
Lichiheb, N., Bedos, C., Personne, E., Benoit, P., Bergheaud, V., Fanucci, O., Bouhlel, J., 
	Barriuso, E., 2015. Measuring Leaf Penetration and Volatilization of Chlorothalonil and 
	Epoxiconazole Applied on Wheat Leaves in a Laboratory-Scale Experiment. 
	Journal of Environment Quality 44, 1782. https://doi.org/10.2134/jeq2015.03.0165

	
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "variables_parametres.h"
#include "fonctions.h"
#include <string.h>


extern double
calcul_PTr(double, double)
, calcul_pph2o(double)
, calcul_Pprime(double)
, calcul_latente(double)
;
extern void Calcul_Rn(double);


void main()
{

	int a, b; // declaration and initialization for header readings of parameter and data files. Usefulness for the cutting of the 1 or 2 first lines of the input files. 
	a = 0;
	b = 0;

	int c; // Variable necessary to detect the end of file during the lecture. Contains the first character of each line of the data file


	//Open input files
	FILE* global_parameter_file = NULL;
	global_parameter_file = fopen("../data/global_parametres.txt", "r");
	if (global_parameter_file != NULL)
		printf("Global parameter file opened.txt\n");
	else
	{
		printf("Impossible to open global parameter file. Check file directory and file name\n");
		getchar();
	}

	FILE* simulation_options_file = NULL;
	simulation_options_file = fopen("../data/simulation_options.txt", "r");
	if (simulation_options_file != NULL)
		printf("Simulation options file opened.txt\n");
	else
	{
		printf("Impossible to open simulation options file. Check file directory and file name\n");
		getchar();
	}
	
	FILE* O3_parameter_file = NULL;
	O3_parameter_file = fopen("../data/O3_parametres.txt", "r");
	if (O3_parameter_file != NULL)
		printf("O3 parameter file opened.txt\n");
	else
	{
		printf("Impossible to open O3 parameter file. Check file directory and file name\n");
		getchar();
	}

	FILE* Pesticides_parameter_file = NULL;
	Pesticides_parameter_file = fopen("../data/Pesticides_parametres.txt", "r");
	if (Pesticides_parameter_file != NULL)
		printf("Pesticides parameter file opened.txt\n");
	else
	{
		printf("Impossible to open Pesticides parameter file. Check file directory and file name\n");
		getchar();
	}

	FILE* NH3_parameter_file = NULL;
	NH3_parameter_file = fopen("../data/NH3_parametres.txt", "r");
	if (NH3_parameter_file != NULL)
		printf("NH3 parameter file opened.txt\n");
	else
	{
		printf("Impossible to open NH3 parameter file. Check file directory and file name\n");
		getchar();
	}

	FILE* data_file = NULL;
	data_file = fopen("../data/data.txt", "r");
	if (data_file != NULL)
		printf("Data file opened.txt\n");
	else
	{
		printf("Impossible to open data file. Check file directory and file name\n");
		getchar();
	}

	FILE* output_file = NULL;
	output_file = fopen("../output/resultats.csv", "w");
	if (output_file != NULL)
		printf("Output file created.txt\n");
	else
	{
		printf("Impossible to create output file. Check file directory and file name\n");
		getchar();
	}
	

	for (a = 0 ; a < 3 ; a++)
	{
		sscanf(fgets(headers, 800000, global_parameter_file), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &delta_time, &lf, &Soil_Depth, &depth_drysoil_threshold, &rho_soil, &attenuationray, &wind_attenuation, &albedo_veg_green, &albedo_veg_yellow, &albedo_soil, &emissivity_canopy, &teta_ini, &teta_minthreshold, &teta_fieldcapacity, &teta_wiltingpoint, &Soil_Saturat_teta, &Soil_Residual_teta, &nVG, &alphaVG, &Conduct_thermal_wetsoil, &Conduct_thermal_drysoil, &depth_soil_Tsoilref, &soil_porosity, &tortuosity_coeff, &Coefficient_LAI_efficient, &alphaglight, &gmax, &gmin, &Tmin, &Tmax, &Topt, &VPDmax, &VPDmin, &SWPmax, &SWPmin);
		sscanf(fgets(headers, 800000, simulation_options_file), "%d %d %d", &Flag_Rn, &Flag_hydric_stress, &Flag_copy);
		sscanf(fgets(headers, 800000, O3_parameter_file), "%lf %lf %lf", &kcut_O3, &RH0_O3, &clay_content);
		sscanf(fgets(headers, 800000, Pesticides_parameter_file), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %s %lf %lf", &Vaporization_enthalpy_Pest, &ex_Pest_Tref, &diffusivity_Pest, &KH_Pest_Tref, &Tref_Pest, &Molar_mass_Pest, &Solub_Water_Pest, &Solub_Octanol_Pest, &logKow_Pest, &K_photodegrad_Pest_ref, &K_washoff_Pest, &K_leaf_dissip_Pest, &K_stom_absorption_Pest, &Lambda_Pest, &Pest_Application_Date, &Pest_Application_Time, &Pest_quantity_applied, &Rcut_Pest);
		sscanf(fgets(headers, 800000, NH3_parameter_file), "%lf %lf", &Rcutmin_NH3, &kcut_NH3);
	}


	initialization();

	while ((c = fgetc(data_file)) != EOF)
	{		
		b = b + 1;
		sscanf(fgets(headers, 800000, data_file), "%s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &Date, &Time, &Rg, &Ratm, &Rn_measured, &T_soilref, &Ta, &RHa, &uref, &LAI_total, &LAI_green, &displacement_height, &z0_soil, &z0_canopy, &zh, &zref, &Rain, &concentration_O3, &concentration_Pest, &Pest_soil, &concentration_NH3, &GammaLeaf_NH3, &GammaSoil_NH3, &R_litter_NH3);
		sprintf(&Date_Time, "%c%s %s", c, Date, Time); // "c" is concatenated with "Date" to have the complete date (see comment concerning the description of variable "c" above)


		if (b == 1) // Reading the header line --> no calculations
			fprintf(output_file, "Date_Time,Rg,Rn_measured,Rn_calcultated,Rn_used,T_soilref,Ta,RHa,uref,LAI_total,LAI_green,displacement_height,z0_soil,z0_canopy,zh,zref,Rain,concentration_O3,concentration_Pest,Pest_soil,concentration_NH3,GammaLeaf_NH3,GammaSoil_NH3,R_litter_NH3,"
				"H,LE,G,Hv,Hs,LEv,LEs,"
				"T_leaf,T_canopy,T_soil,u_star,u_star_ground,RH_canopy,RH_leaf,RH_soil,"
				"Ra,Rac,Rb_leaf,Rb_soil,r_soil_surface_H2O,r_thermal_drylayer,"
				"depth_drysoil,"
				"WetLayer_Depth,teta_WetLayer,teta_minthreshold,teta_tot,Qmulch,SWC_WetLayer,SWC_SoilTotal,SWP,Cumul_Rain_t3,Cumul_evap_t3,"
				"gstom_tot_H2O,glight,gtemp,gVPD,gSWP,gSWC,"
				"FO3_soil,FO3_cut,FO3_stom,FO3_stom_inactive,FO3_canopy,FO3_tot,O3_canopy,O3_leaf,"
				"FPest_soil,FPest_leaf_penetration_adsorb,FPest_canopy,FPest_leaf_dissip,FPest_washoff,FPest_tot,Pest_canopy,Pest_quantity_leaf,Pest_quantity_leaf_Tissue,Pest_quantity_leaf_adsorb,Pest_quantity_leaf_non_adsorb,"
				"FNH3_soil,FNH3_cut,FNH3_stom,FNH3_canopy,FNH3_tot,NH3_canopy,NH3_leaf,NH3_soil\n");
			
		else if (b == 2) // Reading the unit line--> no calculations
			fprintf(output_file, "[GMT],[W.m-2],[W.m-2],[W.m-2],[W.m-2],[°C],[°C],[%],[m.s-1],[m2.m-2],[m2.m-2],[m],[m],[m],[m],[m],[mm],[µg.m-3],[µg.m-3],[µg.m-3],[µg.m-3],[-],[-],[s.m-1],"
				"[W.m-2],[W.m-2],[W.m-2],[W.m-2],[W.m-2],[W.m-2],[W.m-2],"
				"[°C],[°C],[°C],[m.s-1],[m.s-1],[%],[%],[%],"
				"[s.m-1],[s.m-1],[s.m-1],[s.m-1],[s.m-1],[s.m-1],"
				"[m],"
				"[m],[kg.kg-1],[kg.kg-1],[kg.kg-1],[mm],[mm],[mm],[MPa],[mm],[mm],"
				"[m.s-1],[-],[-],[-],[-],[-],"
				"[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-3],[µg.m-3],"
				"[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-3],[µg.m-2],[µg.m-2],[µg.m-2],[µg.m-2],"
				"[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-2.s-1],[µg.m-3],[µg.m-3],[µg.m-3]\n");
		

		else // Excecute Surfatm (3rd line of the input file, i.e., at the 1st line with data)
		{
			
			/*------------------------------------------------------------------
			This "if" condition is here to break at a specific date and hour
			and print some variables and parameters on the screen. 
			Adapt and uncomment if you want to use it
			-------------------------------------------------------------------*/
			
			if (strcmp("26/07/2008 00:00", Date_Time) == 0)
			{
				printf(" \n\n ** %f", diffusivity_O3/diffusivity_H2O);
				printf(" ** %f", diffusivity_NH3 / diffusivity_H2O);
				printf(" ** %f", Sc);
				printf(" ** %f", Sc_O3);
				printf(" ** %f", Qmulch);
				printf(" **  %s\n", Date_Time);

				getchar();
				}
			

				// First resolution under neutral conditions

				Data_preparation();
				Unit_conversion();
				
				Calcul_Rn(Ta);	// if Rn is available as input data, there is no calculation. If not available, Rn is calculated, but check parameter simulation file to ensure it

				Rn_partitionning();
				Calcul_resistances_air_neutral();
				Calcul_stomatal_resistance();
				Calcul_resistances_soil();
				Pprime1 = Pprime2 = Pprime3 = calcul_Pprime(Ta);
				Energy_balance();

				/* Initialization and calculation for the Pprime1,2,3 integrating the new Tcanopy instead of Pprime1 = Pprime2 = Pprime3 */
				critPprime = 0.0;
				incrPprime = 0;
				do
				{
					intermediatetemperature = T_canopy;
				
					Calcul_Rn(T_canopy);  // if Rn is available as input data, there is no calculation. If not available, Rn is calculated, but check parameter simulation file to ensure it
					
					Calcul_deltaprime123();
					Energy_balance();
					incrPprime++;
					critPprime = fabs(T_canopy - intermediatetemperature);
					if (incrPprime > incrmaxPprime)
						break;
				} while ((critPprime > critmaxPprime));
				if ((critPprime > critmaxPprime))
					printf("non convergence Pprime :%s abs T_canopy-T_canopyintemediate :%f incr :%d\n", Date_Time, critPprime, incrPprime);
				

				/* Iterative process for stability/unstability of the atmosphere */
				if (neutral == 0) // index to check atmospherci neutrality : 1 = neutral || 0 = stable or unstable [-]
				{
					/* Initialization of the iterative process :*/
					critStab = 0.0;
					incrStab = 0;
					do
					{
						/* Initialization of the temperature close to the "top" of the canopy */
						intermediatetemperature = T_canopy;

						Calcul_Rn(T_canopy);  // if Rn is available as input data, there is no calculation. If not available, Rn is calculated, but check parameter simulation file to ensure it

						Calcul_deltaprime123();
						Calcul_resistances_air_stab();
						Calcul_stomatal_resistance();
						Calcul_resistances_soil();
						Energy_balance();
						incrStab++;
						critStab= fabs(T_canopy - intermediatetemperature);
						if (incrStab > incrmaxStab)
							break;
					} while ((critStab > critmaxStab));
					if ((critStab > critmaxStab))
						printf("non convergence of stability at :%s abs T_canopy-T_canopyintemediate :%f incr :%d\n", Date_Time, critStab, incrStab);
				}

				SWC_evolution();

				Flux_O3();

				Flux_Pesticide();

				Flux_NH3();

				printf(" **  %s\n", Date_Time);


				/* ************* */

				fprintf(output_file,"%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
									Date_Time, Rg, Rn_measured, Rn_calculated, Rn, T_soilref, Ta, RHa, uref, LAI_total, LAI_green, displacement_height, z0_soil, z0_canopy, zh, zref, Rain, concentration_O3, concentration_Pest, Pest_soil, concentration_NH3, GammaLeaf_NH3, GammaSoil_NH3, R_litter_NH3,
									H, LE, G, Hv, Hs, LEv, LEs, 
									T_leaf, T_canopy, T_soil, u_star, u_star_ground, RH_canopy, RH_leaf, RH_soil, 
									Ra, Rac, Rb_leaf, Rb_soil, r_soil_surface_H2O, r_thermal_drylayer, 
									depth_drysoil, 
									WetLayer_Depth, teta_WetLayer, teta_minthreshold, teta_tot, Qmulch, SWC_WetLayer, SWC_SoilTotal, SWP, Cumul_Rain_t3, Cumul_evap_t3, 
									gstom_tot_H2O, glight, gtemp, gVPD, gSWP,gSWC,
									FO3_soil, FO3_cut, FO3_stom, FO3_stom_inactive, FO3_canopy, FO3_tot, O3_canopy, O3_leaf,
									FPest_soil, FPest_leaf_penetration_adsorb, FPest_canopy, FPest_leaf_dissip, FPest_washoff, FPest_tot, Pest_canopy, Pest_quantity_leaf, Pest_quantity_leaf_Tissue, Pest_quantity_leaf_adsorb, Pest_quantity_leaf_non_adsorb,
									FNH3_soil, FNH3_cut, FNH3_stom, FNH3_canopy, FNH3_tot, NH3_canopy, NH3_leaf, NH3_soil
									);
		}
	

	}

	fclose(global_parameter_file);
	fclose(simulation_options_file);
	fclose(O3_parameter_file);
	fclose(Pesticides_parameter_file);
	fclose(NH3_parameter_file);
	fclose(data_file);
	fclose(output_file);


	
	
	time_t timestamp = time(NULL);
	strftime(buffer_time, sizeof(buffer_time), "%Y_%m_%d_%HH%M", localtime(&timestamp));

	printf("\n\n **********   Finished run. Congratualations **************\n\n");
	printf("%s\n", buffer_time);
	getchar();

	if (Flag_copy == 1)
	filecopy();


}


