/*----------Evolution of the Soil Water Content--------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "variables_parametres.h"

extern double relationtetafi(double, double, double);


void SWC_evolution(void)
{

	double
		evapmm	= 0.0		// Evaporation in mm                        [mm ou kg/m²]
		, transpmm	= 0.0	// Transpiration in mm                      [mm ou kg/m²]
		, Mulch_Max_SWC	= 0.0 //											[mm ou kg/m²]
		;

// Conversion from [W/m2] to [kg/m2/timestep]
	evapmm = LEs / 2460000.0 * 3600.0 * delta_time; // Soil 
	transpmm = LEv / 2460000.0 * 3600.0 * delta_time; // Plant
	if (evapmm   < 0.0) evapmm = 0.0; 
	if (transpmm < 0.0) transpmm = 0.0; 
	
// Memory of the last 3 rainfall events, and their cumulated sum on the last 3 timesteps
	Rain_t3 = Rain_t2; 
	Rain_t2 = Rain_t1;
	Rain_t1 = Rain ;	
	Cumul_Rain_t3 = Rain_t3 + Rain_t2 + Rain_t1;

// Memory of the last 3 evaporation amounts, and their cumulated sum on the last 3 timesteps
	evap_t3 = evap_t2;
	evap_t2 = evap_t1;
	evap_t1 = evapmm;
	Cumul_evap_t3 = evap_t3 + evap_t2 + evap_t1;

	
	SWC_SoilTotal = SWC_SoilTotal + Rain - evapmm - transpmm;
	teta_tot = SWC_SoilTotal / (Soil_Depth * rho_soil);

	
	if (teta_tot <= teta_minthreshold || teta_WetLayer <= teta_minthreshold)
	{
		teta_tot = teta_minthreshold;
		teta_WetLayer= teta_minthreshold;
		SWC_SoilTotal = teta_minthreshold * Soil_Depth * rho_soil;
		evapmm = 0.0;
		transpmm = 0.0;
		LEs = 0.0;
		LEv = 0.0;
		r_soil_surface_H2O = 100000.0;
	}
	
	Mulch_Max_SWC = depth_drysoil * rho_soil * (teta_fieldcapacity - teta_minthreshold); // Maximum amount of water that can be contained by the mulch

	
	if ((Cumul_Rain_t3 - Cumul_evap_t3) > (Mulch_threshold * Mulch_Max_SWC))   // Rainfall event (at least on rainfall on the last 3 last timesteps) which fill the muclh at at least XX% (see value of Mulch_threshold in initialization.c)
	{
		depth_drysoil = 0.0;
		WetLayer_Depth = Soil_Depth;
		teta_WetLayer = teta_tot;
		SWC_WetLayer = SWC_SoilTotal;
		Qmulch = 0.0;
	}

	else if (depth_drysoil >= depth_drysoil_threshold)
	{
		depth_drysoil = depth_drysoil_threshold;
		Qmulch = depth_drysoil * teta_minthreshold * rho_soil;
		WetLayer_Depth = Soil_Depth - depth_drysoil;
		SWC_WetLayer -= (evapmm - Rain);
//		SWC_WetLayer = SWC_SoilTotal - Qmulch; 
	}

	else
	{
		DryWet_depth_change = (evapmm - Rain) / ((teta_WetLayer - teta_minthreshold)* rho_soil);
		depth_drysoil += DryWet_depth_change;
		if (depth_drysoil <= 0.0)
			depth_drysoil = 0.0;
		//WetLayer_Depth -= DryWet_depth_change; 
		WetLayer_Depth = Soil_Depth - depth_drysoil;
		Qmulch = depth_drysoil * teta_minthreshold * rho_soil;
		SWC_WetLayer = WetLayer_Depth * teta_WetLayer * rho_soil;
		if (depth_drysoil > Soil_Depth)
			depth_drysoil = Soil_Depth;
	} // End of the evaporation phase

	SWC_WetLayer -= transpmm;
	teta_WetLayer = SWC_WetLayer / (WetLayer_Depth * rho_soil);

	if (SWC_WetLayer > SWC_Fieldcapacity)
	{
		Drain_Runoff = SWC_WetLayer - SWC_Fieldcapacity;
		SWC_WetLayer = Soil_Depth * rho_soil * teta_fieldcapacity;
		SWC_SoilTotal = Soil_Depth * rho_soil * teta_fieldcapacity;
		teta_WetLayer = teta_fieldcapacity;
		teta_tot = teta_fieldcapacity;
		depth_drysoil = 0.0;
		Qmulch = 0.0;
		WetLayer_Depth = Soil_Depth;
	}
	else
	{
		Drain_Runoff = 0.0;
	}

	if (SWC_SoilTotal <= (teta_minthreshold*Soil_Depth*rho_soil))
	{
		printf("warning : check LE coherence --> soil wilting point threshold reached\n");
	}
		
	SWP = relationtetafi(teta_WetLayer, SWP, rho_soil);

}
