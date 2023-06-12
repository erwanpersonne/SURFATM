/*----------Program for calculation of air, stomatal and soil resistances--------

The resistances calculated here only concern the transfers of heat, 
water vapor and momentum





References:
Thom, A.S. (1975) Momentum, mass and heat exchange of plant communities. In: J.L.
	Monteith (Editor), Vegetation and the atmosphere. Academic Press, New York, pp. 57-109.
Choudhury, B.J., Monteith, J.L. (1988) A four-layer model for the heat budget of homogeneous
	land surfaces, Quarterly Journal of the Royal Meteorological Society, 114, 373-398.
Lagouarde, J.-P., et al. (1995) Le bilan d'énergie d'un couvert végétal,In: INRA (Editor),
	Actes de l'Ecole-Chercheurs INRA en bioclimatologie, pp. 383-404.
Emberson, L.D., Simpson, D., Tuovinen, J.P., Ashomore, M.R., Cambridge, H.M. (2000)
	Towards a model of ozone deposition and stomatal uptake over Europe. EMEP, Research
	Note No. 42.
Jarvis, P.G. (1976) The interpretation of the variations in leaf water potential and
	stomatal conductance found in canopies in the field. Philosophical Transactions of
	the Royal Society of London Series B, 273, 593-610.
Jones, H.G. (1992) Plants and Microclimate: A quantitative Approach to Environmental
	Plant Physiology. 2nd Edition. Cambridge University Press, Cambridge.
Tuovinen, J.P., Ashmore, M.R., Emberson, L.D., Simpson, D. (2004) Testing and
	improving the EMEP ozone deposition module. Atmospheric Environment, 38, 2373-2385.
Zhou, M.C., Ishidaira; H., Hapuarachchi, H.P., Magome, J., Kiem, A.S., Takeuchi, K. (2006)
	Estimating potential evapotranspiration using Shuttleworth-Wallace model and
	NOAA-AVHRR NDVI data to feed a distributed hydrological model over the Mekong River
	Basin. Journal of Hydrology, 327, 151-173.
Stella, P., Personne, E., Loubet, B., Lamaud, E., Ceschia, E., Béziat, P.,
	Bonnefond, J.M., Irvine, M., Keravec, P., Mascher, N., Cellier, P. (2011)
	Predicting and partitioning ozone fluxes to maize crops from sowing to harverst:
	the SURFATM-O3 model. Biogeosciences, 8, 2869-2886.
Rochette, P., Pattey, E., Desjardins, R.L., Dwyer, L.M., Stewart, D.W., Dubé, P.A.
	(1991) Estimation of maize (Zea mays L.) canopy conductance by scaling up leaf
	stomatal conductance, Agricultural and Forest Meteorology, 54, 241–261.
J. Farahani, H., R. Ahuja, L., 1996. Evapotranspiration Modeling of Partial Canopy/Residue-covered Fields.
	Transactions of the ASAE 39, 2051–2064. https://doi.org/10.13031/2013.27708
Griend, A.A. van de, Owe, M., 1994. Bare soil surface resistance to evaporation by vapor diffusion under semiarid conditions.
	Water Resources Research 30, 181–188. https://doi.org/10.1029/93WR02747

*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "variables_parametres.h"

void Calcul_resistances_air_neutral()
{
	
	// calculation of windspeed profile and friction velocity 
	u_star = karman*uref / log((zref - displacement_height) / z0_canopy);
	if (u_star < 0.2)
		u_star = 0.2;
	Keh = karman * u_star * (zh - displacement_height);
	uh = (u_star / karman) * (log((zh - displacement_height) / z0_canopy));
	

	// Calculation of Ra, Rb et Rac in neutral condition
	Ra = pow(log((zref - displacement_height) / z0_canopy), 2.0) / (karman * karman * uref);
	if (Ra > 3000.0)
		Ra = 3000.0;

	Rb_leaf = 1.0 / (LAI_total * 2 * 0.01 / wind_attenuation * pow(uh / lf, 0.5) * (1.0 - exp(-wind_attenuation / 2.0)));
	if (Rb_leaf > 3000.0)
		Rb_leaf = 3000.0;

	u_star_ground = u_star * u_star * exp(1.2 * LAI_total * (z0_soil / zh - 1.0));
	u_star_ground = pow(u_star_ground, 0.5);
	if (zh < 0.03)
		u_star_ground = u_star;

	Rb_soil = 2.0 / (karman * u_star_ground) * pow((Sc / Prandt), 2 / 3);
	if (Rb_soil > 3000.0)
		Rb_soil = 3000.0;

	Rac = zh * exp(wind_attenuation) / (wind_attenuation * Keh)*(exp(-wind_attenuation * z0_soil / zh) - exp(-wind_attenuation * (displacement_height + z0_canopy) / zh));
	if (Rac <= 0.01)
		Rac = 0.01;
	if (Rac > 1000.0)
		Rac = 1000.0;

}

void Calcul_resistances_air_stab()
// u_star_ground, Rb_soil, Keh et Rac are not recalculated, we implictly consider the neutrality inside the canopy
{
// Atmospheric stability
	L_MO = -(0.5*(Ta + T_canopy) + T0C) / g * pow(u_star, 3.0) / karman * rho_air * Cp_air / H;
	Zeta = (zref - displacement_height) / L_MO;
	x = pow(1 - 16 * Zeta, 0.25);
	if (Zeta < 0) // Unstable
	{
		PsiH = 2 * log((1 + x*x) / 2);
		PsiM = log((1 + x*x) / 2) + 2 * log((1 + x) / 2) - 2 * atan(x) + M_PI / 2;
		stab = 0;
	}
	else           // Stable
	{
		PsiH = max(-5.0, -5.2*Zeta); 
		PsiM = max(-5.0, -5.2*Zeta);
		stab = 1;
	}

// Ra, Rb_canopy, u_star and uh accounting for atmospheric stability
	u_star = (karman * uref) / (log((zref - displacement_height) / (z0_canopy)) - PsiM);
	if (u_star < 0.2)
		u_star = 0.2;
	uh = (u_star / karman) * ((log((zh - displacement_height) / z0_canopy)) - PsiM);
	if (LAI_total <= 0.5)
		uh = (u_star / karman) * (log(zh / z0_soil) - PsiM);
	if (uh < 0.1)
		uh = 0.1;

	u_star_ground = u_star * u_star * exp(1.2 * LAI_total * (z0_soil / zh - 1.0));
	u_star_ground = pow(u_star_ground, 0.5);
	if (zh < 0.03)
		u_star_ground = u_star;

	
	Ra = (log((zref - displacement_height) / (z0_canopy)) - PsiH) * (log((zref - displacement_height) / (z0_canopy)) - PsiM) / (karman  *karman * uref);	
	if (Ra > 2500.0)
		Ra = 2500.0;

	Rb_leaf = 1.0 / (LAI_total * 2 * 0.01 / wind_attenuation * pow(uh / lf, 0.5) * (1.0 - exp(-wind_attenuation / 2.0)));
	if ((zh - z0_soil) <= 0.03)
		Rb_leaf = 3000.0;
	if (Rb_leaf > 3000.0)
		Rb_leaf = 3000.0;

// Comparaison between Ra under neutral condition et Ra with atmospheric stability correction	
	diff_Ra = (pow(log((zref - displacement_height) / z0_canopy), 2.0) / (karman * karman * uref)) - Ra;
	if (abs(diff_Ra) > 300.0)
		Ra = pow(log((zref - displacement_height) / z0_canopy), 2.0) / (karman * karman * uref);
}

void Calcul_stomatal_resistance()
{
// Calculation of efficient LAI, i.e., fraction of the canopy contributing mainly to stomatal exchanges
	LAI_green_efficient = Coefficient_LAI_efficient * LAI_green;
	LAI_yellow_efficient = Coefficient_LAI_efficient * LAI_yellow;
	LAI_total_efficient = Coefficient_LAI_efficient * LAI_total;


// Leaf stomatal conductance from Jarvis multiplicative model
	glight = 1.0 - exp(-alphaglight * PAR);
	
	if (Ta < Tmin)
		gtemp = gmin;
	else if (Ta > Tmax)
		gtemp = gmin;
	else
		gtemp = 1 - pow((Ta - Topt) / (Topt - Tmin), 2);

	if ((VPDa / 1000.0) < VPDmax) // Warning !! VPDa is in Pa while VPDmax is in kPa
		gVPD = 1;
	else if ((VPDa / 1000.0) > VPDmin) // Warning !! VPDa is in Pa while VPDmin is in kPa
		gVPD = gmin;
	else
		gVPD = ((1 - gmin) / (VPDmax - VPDmin)) * (VPDa / 1000.0) + (gmin - ((1 - gmin) / (VPDmax - VPDmin)) * VPDmin); // Warning !! VPDa is in Pa while VPDmax and VPDmin are in kPa


	if (Flag_hydric_stress == 1)
	{
		gSWC = 1.0;
		if (SWP < SWPmax)
			gSWP = gmin;
		else if (SWP > SWPmin)
			gSWP = 1.0;
		else
			gSWP = ((1 - gmin) / (SWPmin - SWPmax))*SWP + (gmin - ((1 - gmin) / (SWPmin - SWPmax))*SWPmax);
	}

	if (Flag_hydric_stress == 0)
	{
		gSWP = 1.0;
		if (teta_WetLayer < teta_wiltingpoint)
			gSWC = gmin;
		else if (teta_WetLayer > teta_fieldcapacity)
			gSWC = 1.0;
		else
			gSWC = (1.0 - gmin) / (teta_fieldcapacity - teta_wiltingpoint) * teta_WetLayer + (1.0 - (teta_fieldcapacity * (1.0 - gmin) / (teta_fieldcapacity - teta_wiltingpoint)));		
	}

	gstom_H2O_leaf = (gmax * glight * max(gmin, gtemp * gVPD * gSWP * gSWC)) / 41000.0;

// Stomatal closure if soil hydric stress, wathever the condition	
	if (teta_WetLayer <= teta_wiltingpoint)
		gstom_H2O_leaf = 0.000001;

// Canopy stomatal conductances, i.e., upscaling from leaf to canopy
	gstom_H2O = LAI_green_efficient * gstom_H2O_leaf;				// Canopy stomatal conductance of green leaves
	Rstom_H2O = 1.0 / gstom_H2O;

	gstom_inactive_H2O = LAI_yellow_efficient * gstom_H2O_leaf;		// Canopy stomatal conductance of yellow leaves
	Rstom_inactive_H2O = 1.0 / gstom_inactive_H2O;

	gstom_tot_H2O = LAI_total_efficient * gstom_H2O_leaf;			// Total canopy stomatal conductance (green + yellow leaves). Required for LE vegetation calculation
	Rstom_tot_H2O = 1.0 / gstom_tot_H2O;
	
// Control of the values of resistances/conductances for very low and sparse canopies
	if (LAI_green <= 0.001 || gstom_H2O_leaf <= 0.0)
	{
		gstom_H2O = 0.000001;
		Rstom_H2O = 1000000.0;
	}

	if (LAI_yellow <= 0.001 || gstom_H2O_leaf <= 0.0)
	{
		gstom_inactive_H2O = 0.000001;
		Rstom_inactive_H2O = 1000000.0;
	}
	
	if (LAI_total <= 0.001 || gstom_H2O_leaf <= 0.0)
	{
		gstom_tot_H2O = 0.000001;
		Rstom_tot_H2O = 1000000.0;
	}

}

void Calcul_resistances_soil()
{


	// Resistance to heat transfert in the soil

	r_thermal_wetlayer = rho_air*Cp_air * (depth_soil_Tsoilref - depth_drysoil) / Conduct_thermal_wetsoil;
	r_thermal_drylayer = rho_air*Cp_air * depth_drysoil / Conduct_thermal_drysoil;

	depth_threshold_drysoil = 0.00;

	if (depth_drysoil < depth_threshold_drysoil)
		r_soil_surface_H2O = 0.0;

	else
	{
		r_soil_surface_H2O = tortuosity_coeff * (depth_drysoil - depth_threshold_drysoil) / (soil_porosity * diffusivity_H2O);
		if (depth_drysoil > 0.05)
		r_soil_surface_H2O = r_soil_surface_H2O * exp(-beta_soil * (((teta_minthreshold + teta_WetLayer) / 2.0) * rho_soil/1000.0)/teta_fieldcapacity); // Inspired from Adriaan et Van de Griend 1994: Farahani et Ahuja 1996 - importance of if 0.05
		r_soil_surface_H2O = (r_soil_surface_H2O < 0) ? 0 : r_soil_surface_H2O;
	}

	if (depth_drysoil >= depth_drysoil_threshold) //depth_threshold_wetsoil)
	{
		r_soil_surface_H2O = 10000.0;
		printf("probable problem of evaporation\n");
	}

}








