/*----------Program for O3 flux calculations--------

References:
Lamaud, E., Loubet, B., Irvine, M., Stella, P., Personne, E., Cellier, P.
	(2009) Partitionning of ozone deposition over a developed maize crop between
	stomatal and non stomatal uptakes, using eddy-covariance flux measurements
	and modelling. Agricultural and Forest Meteorology, 149, 1385-1396.
Stella, P., Lamaud, E., Loubet, B., Laville, P, Cellier, P. (2011) Ozone
	deposition onto bare soil: A new parameterisation. Agricultural and Forest
	Meteorology, 151, 669-681.
Stella, P., Personne, E., Loubet, B., Lamaud, E., Ceschia, E., Béziat, P.,
	Bonnefond, J.M., Irvine, M., Keravec, P., Mascher, N., Cellier, P., 2011.
	Predicting and partitioning ozone fluxes to maize crops from sowing to harverst:
	the SURFATM-O3 model. Biogeosciences, 8, 2869-2886.
Stella, P., Loubet, B., de Berranger, C., Charrier, X., Ceschia, E., Gerosa, G., 
	Finco, A., Lamaud, E., Serça, D., George, C., Ciuraru, R. (2019) Soil ozone 
	deposition: Dependence of soil resistance to soil texture. Atmospheric Environnemnt,
	199, 202-209.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "variables_parametres.h"
#include "fonctions.h"

void O3_resistance_calculations(void);
void O3_Flux_calculations(void);


void Flux_O3()
{
	O3_resistance_calculations();
	O3_Flux_calculations();
}

void O3_resistance_calculations()
{
	// Stomatal conductance of green leaves
	Rstom_O3 = diffusivity_H2O / diffusivity_O3 * Rstom_H2O;
	gstom_O3 = 1.0 / Rstom_O3;

	// Stomatal conductance of yellow leaves - Cf Stella et al. (2013)
	Rstom_inactive_O3 = diffusivity_H2O / diffusivity_O3 * Rstom_inactive_H2O;
	gstom_inactive_O3 = 1.0 / Rstom_inactive_O3;

	// Cuticular conductance - Cf Lamaud et al. (2009)
	Rcutmin_O3 = 5000.00 / LAI_total;
	if (RH_canopy <= RH0_O3)
		gcut_O3 = 1.0 / (Rcutmin_O3);
	else
		gcut_O3 = 1.0 / (Rcutmin_O3*exp(-kcut_O3*(RH_canopy - RH0_O3)));
	if (RH_canopy >= 100)
		gcut_O3 = 1.0 / (Rcutmin_O3*exp(-kcut_O3*(100.0 - RH0_O3)));

	// Vegetation conductance
	g_leaf_O3 = gcut_O3 + gstom_O3 + gstom_inactive_O3;
	R_leaf_O3 = 1.0 / g_leaf_O3;

	// Canopy conductance (vegetation and its boundary layer)
	Rb_leaf_O3 = diffusivity_H2O / diffusivity_O3 * Rb_leaf;
	g_canopy_O3 = 1.0 / (Rb_leaf_O3 + R_leaf_O3);

	// Soil conductance - Cf Stella et al. (2011)
	Rsoilmin_O3 = 702.0 * pow(clay_content, -0.98);
	ksoil_O3 = 0.0118 * exp(0.0266 * clay_content);
	Rsoil_O3 = Rsoilmin_O3 * exp(ksoil_O3 * RH_soil);
	Rb_soil_O3 = 2.0 / (karman * u_star_ground) * pow((Sc_O3 / Prandt), 2 / 3);
	if (Rb_soil_O3 > 3000.0)
		Rb_soil_O3 = 3000.0;
	gsoil_O3 = 1.0 / (Rac + Rb_soil_O3 + Rsoil_O3);

	// Total conductance (soil + vegetation)
	gc_O3 = gsoil_O3 + g_canopy_O3;
	Rc_O3 = 1.0 / gc_O3;

	// Controle of the values of the resistances/conductances for very low canopies
	if ((zh - z0_soil) <= 0.03)
	{
		gstom_O3 = 0.000001;
		Rstom_O3 = 1000000.0;
		gstom_inactive_O3 = 0.000001;
		Rstom_inactive_O3 = 1000000.0;
		gcut_O3 = 0.000001;
		Rb_leaf_O3 = 1000000.0;
		g_leaf_O3 = 0.000001;
		R_leaf_O3 = 1000000.0;
		g_canopy_O3 = 0.000001;
		gc_O3 = gsoil_O3;
		Rc_O3 = 1.0 / gc_O3;
	}
}
void O3_Flux_calculations()
{
	// Deposition velocity
	Vd_O3 = 100.0 / (Ra + Rc_O3);

	// Flux calculations
	// The O3 concentrations at the different levels are calculated from flux-profile relationships (e.g., FO3tot = (O3canopy - O3air)/Ra ) 
	FO3_tot = -concentration_O3 * Vd_O3 / 100.0;
	O3_canopy = FO3_tot * Ra + concentration_O3;
	FO3_soil = -gsoil_O3 * O3_canopy;
	FO3_canopy = FO3_tot - FO3_soil;
	O3_leaf = FO3_canopy * Rb_leaf_O3 + O3_canopy;
	FO3_cut = -gcut_O3 * O3_leaf;
	FO3_stom = -gstom_O3 * O3_leaf;
	FO3_stom_inactive = -gstom_inactive_O3 * O3_leaf;

	// Controle of the values of the fluxes/concentrations for very low canopies
	if ((zh - z0_soil) <= 0.03)
	{
		O3_canopy = 0.000001;
		FO3_soil = FO3_tot;
		FO3_canopy = -0.000001;
		O3_leaf = 0.000001;
		FO3_cut = -0.000001;
		FO3_stom = -0.000001;
		FO3_stom_inactive = -0.000001;
	}
}