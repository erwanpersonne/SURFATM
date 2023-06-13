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


/*-----------------Declarations of the constants and their values--------------*/

#define M_PI			3.1415926535
#define diffusivity_H2O_0C	2.178e-5		// H2O molecular diffusivity in air at 0°C and 1 atm	[m2/s] cf Massman (1998)
#define diffusivity_O3_0C  1.444e-5			// O3 molecular diffusivity in air at 0°C and 1 atm		[m2/s] cf Massman (1998)
#define diffusivity_NH3_0C  1.978e-5		// NH3 molecular diffusivity in air at 0°C and 1 atm	[m2/s] cf Massman (1998)
#define karman			0.41				// von Karman constant								[-]
#define Prandt			0.72				// Prandt constant for air							[-]
#define T0C				273.15				// Kelvin temperature at 0°C					    [°C]
#define g				9.81				// Gravitational constant		                    [m/s2]
#define rho_air			1.164				// Air density					                    [kg/m3]
#define Cp_air			1010.0				// Air specific heat			                    [kJ/kg/K]
#define R				8.318				// Perfect gaz constant								[Pa.m3.mol-1.K-1]
#define PsychroCst		66.5				// psychrometric constant at 101.325kPa and 0°C		[Pa.k-1]
#define	Viscosity_cin	1.56e-5				// Cinematic viscosity of air						[m2/s]
#define KHenry25C_NH3	7.24e-4				// Henry constante de Henry at 25°C					[-]
#define KAcidB25C_NH3  5.62e-7				// dissociation constante Acide Base at 25°C		[mol/m3]
#define enthalpVap_NH3  34.18e3				// free enthalpy for vaporization of NH3			[J/mol]
#define enthalpAcidB_NH3  52.21e3			// free enthalpy fo acid-base dissociation			[J/mol]
#define SIGMA			5.67e-8				// Stefan-Boltzmann Constant						[W.m-2.K-4]


/*-----------------Declaration of the flags for simulation options --------------*/
int
	Flag_Rn,				// Flag to choose Rn measured (= 1) or calculated ( = 0)
	Flag_hydric_stress,		// Flag to choose how to calculate hydric stress on gstom: 1 = with SWP || 0 = with SWC
	Flag_copy;				// Flag to choose to copy (= 1) or not (= 0) all the simulation files


/*-----------------Declaration of the variables--------------*/
int
	stab,					// index to check atmospheric stability : 0 = unstable || 1 = stable	[-]
	neutral					// index to check atmospherci neutrality : 1 = neutral || 0 = stable or unstable [-]
	;

char headers[800000];
char Date[150];
char Time[150];
char Date_Time[150];
char Date_Time_Pest_application[150];

double
Rg,						// Shortwave global radiation - downward						[W/m2]
Ratm,					// Long Wave Radiation - downward (equivalent to atmosph. radiation) [W/m2]
RHa,					// Air relative humidity at zref								[%]
Ra,						// Aerodynamic resistance above the canopy						[s/m]
Rac,					// Aerodynamic resistance inside the canopy						[s/m]
Rb_leaf,				// Canopy boundary layer resistance								[s/m]
Rb_leaf_O3,				// Canopy boundary layer resistance for ozone					[s/m]
Rb_leaf_NH3,			// Canopy boundary layer resistance for NH3						[s/m]
Rb_soil,				// Soil boundary layer resistance								[s/m]
Rb_soil_O3,				// Soil boundary layer resistance for ozone						[s/m]
Rb_soil_NH3,			// Soil boundary layer resistance for NH3						[s/m]
Rstom_H2O,				// Green leaves stomatal resistance for H2O						[s/m]
gstom_H2O,				// Green leaves stomatal conductance for H2O					[s/m]
gstom_tot_H2O,			// Canopy stomatal conductance for H2O (green + yellow leaves)	[m/s]
Rstom_tot_H2O,			// Canopy stomatal resistance for H2O (green + yellow leaves)	[m/s]
gstom_H2O_leaf,			// Leaf stomatal conductance for H2O							[m/s]
Rstom_O3,				// Green leaves stomatal resistance for ozone					[s/m]
gstom_O3,				// Green leaves stomatal conductance for ozone					[s/m]
Rstom_NH3,				// Green leaves stomatal resistance for NH3						[s/m]
gstom_NH3,				// Green leaves stomatal conductance for NH3					[s/m]
Rstom_inactive_H2O,		// Yellow leaves stomatal resistance for H2O					[s/m]
gstom_inactive_H2O,		// Yellow leaves stomatal conductance for H2O					[s/m]
Rstom_inactive_O3,		// Yellow leaves stomatal resistance for ozone					[s/m]
gstom_inactive_O3,		// Yellow leaves stomatal conductance for ozone					[s/m]
Rstom_inactive_NH3,		// Yellow leaves stomatal resistance for NH3					[s/m]
gstom_inactive_NH3,		// Yellow leaves stomatal conductance for NH3					[s/m]
Rsoil_O3,				// Soil resistance to ozone deposition							[s/m]
gsoil_O3,				// soil conductance to ozone deposition							[m/s]
Rsoil_NH3,				// Soil resistance to NH3 deposition							[s/m]
gsoil_NH3,				// soil conductance to NH3 deposition							[m/s]
R_litter_NH3,			// additional resist. of the litter covering the soil for NH3 exchange	[s/m] 
Rcut_NH3,				// cuticular resistance to NH3 deposition						[s/m]
gcut_NH3,				// cuticular conductance to NH3 deposition						[m/s]
Rcutmin_O3,				// Minimal cuticular resistance to ozone deposition				[s/m]
Rc_O3,					// Canopy resistance to ozone deposition						[s/m]	
gc_O3,					// Canopy conductance to ozone deposition						[m/s]
u_star,					// Friction velocity											[m/s]
u_star_ground,			// Friction velocity near the soil								[m/s]
RH_leaf,				// Relative humidity of the air near the leaf surface			[%]
RH_soil,				// Relative humidity of the air neat the soil surface			[%]
e_soil,					// air vapor pressure at the soil surface						[Pa]
ex_soil,				// air saturated vapor pressure at the soil surface				[Pa]
e_canopy,				// air vapor pressure at the canopy level (above Rb_leaf)		[Pa]
ex_canopy,				// air saturated vapor pressure at the canopy level (above Rb_leaf)		[Pa]
VPD_canopy,				// air vapor pressure deficit at the canopy level (above Rb_leaf)		[Pa]
RH_canopy,				// air relative humidity at the canopy level (above Rb_leaf)			[Pa]
ea,						// air vapor pressure at zref									[Pa]
exa,					// air saturated vapor pressure at zref							[Pa]
VPDa,					// Vapour pressure deficit of the air at zref					[Pa]
T_leaf,					// Leaf surface temperature (below Rb_leaf)						[°C]
T_soil,					// soil surface temperature										[°C]
T_soilref,				// soil temperature at the deeply soil - ie, at "depth_soil_Tsoilref"	[°C]
T_soilwetdry,			// soil temperature at the interface of wet and dry soil layers	[°C]
ex_soilwetdry,			// air saturated vapor pressure at the interface of wet and dry soil layers [Pa]
e_leaf,					// air vapor pressure at the leaf level (below Rb_leaf)			[Pa]
ex_leaf,				// air saturated vapor pressure at the leaf level (below Rb_leaf)		[Pa]
LAI_total,				// Total leaf area index of the canopy							[m2/m2]
LAI_total_efficient,	// Leaf area index of leaves contributing to stomatal exchanges	[m2/m2]
LAI_green,				// Leaf area index of green leaves								[m2/m2]
LAI_green_efficient,	// Leaf area index of green leaves contributing to stomatal exchanges	[m2/m2]
LAI_yellow,				// Leaf area index of yellow leaves								[m2/m2]
LAI_yellow_efficient,	// Leaf area index of yellow leaves contributing to stomatal exchanges	[m2/m2]
zh,						// Canopy height												[m]
canopy_cover,			// Fraction the soil surface covered by canopy (ground cover ratio)	[-]		

concentration_O3,		// Ozone concentration at zref									[µg/m3]
Vd_O3,					// Ozone deposition velocity at zref							[cm/s]
g_leaf_O3,				// Ozone deposition conductance of the leaf (cuticle + stomata) without Rb_leaf	[m/s]
R_leaf_O3,				// Ozone deposition resistance of the leaf (cuticle + stomata) without Rb_leaf	[m/s]
g_canopy_O3,			// Ozone deposition conductance of the leaf (cuticle + stomata) with Rb_leaf	[m/s]
gcut_O3,				// Ozone cuticular deposition conductance						[m/s]
FO3_tot,				// Total ozone flux at zref										[µg/m2/s]
FO3_soil,				// Soil ozone flux												[µg/m2/s]
FO3_stom,				// Ozone flux to green leaves									[µg/m2/s]
FO3_stom_inactive,		// Ozone flux to yellow leaves									[µg/m2/s]
FO3_cut,				// Cuticular ozone flux											[µg/m2/s]
FO3_canopy,				// Ozone flux to canopy (stomata + cuticular)					[µg/m2/s]
O3_canopy,				// Ozone concentration above the canopy (above Rb)				[µg/m3]
O3_leaf,				// Ozone concentration near the leaves surface (below Rb)		[µg/m3]

concentration_Pest,		// Concentration of the studied pesticide at zref				[µg/m3]
Rstom_Pest,				// Green leaves stomatal resistance for studied pesticide		[s/m]
Rstom_inactive_Pest,	// Yellow leaves stomatal resistance for studied pesticide		[s/m]
Rb_soil_Pest,			// Soil boundary layer resistance for studied pesticide			[s/m]
Rb_leaf_Pest,			// Canopy boundary layer resistance for studied pesticide		[s/m]
FPest_soil,				// Soil flux of the studied pesticide							[µg/m2/s]
FPest_canopy,			// Flux to canopy (stomata + cuticular)	 of the studied pesticide	[µg/m2/s]
FPest_tot,				// Total vertical surface-atmosphere flux at zref of the studied pesticide					[µg/m2/s]
FPest_washoff,			// Flux due to washoff process of the pesticid after a rain		[µg/m2/s]
FPest_leaf_photodegrad_adsorb, // Flux due to photodegrad on the adsorbed pesticid on the leaves [µg/m2/s] 
FPest_leaf_photodegrad_non_adsorb, // Flux due to photodegrad on the non-adsorbed pesticid on the leaves [µg/m2/s]
FPest_leaf_penetration_adsorb, // Flux due to penetration of the pesticid located on the adsorbed layer on the leaves [µg/m2/s]
FPest_leaf_penetration_non_adsorb, // Flux due to penetration of the pesticid located on the non-adsorbed layer on the leaves [µg/m2/s]
FPest_leaf_dissip, // Flux due to dissipation of the pesticid contained inside the leaves [µg/m2/s]
Pest_canopy,			// Pesticide concentration above the canopy (above Rb)			[µg/m3]
Pest_leaf,				// Pesticide concentration near the leaves surface (below Rb)	[µg/m3]
Pest_soil,				// Pesticide concentration at the soil surface (below Rb)		[µg/m3]
Sc_Pest,				// Schmidt number for studied pesticide							[-]
Pest_quantity_leaf,		// Quantity of pesticide of the leaf (inside + on the leaves) 	[µg/m2]
Pest_quantity_leaf_adsorb, // Quantity of pesticide adsorbed on the leaf surface and non volatilizable [µg/m2]
Pest_quantity_leaf_non_adsorb, // Quantity of pesticide non adsorbed on the leaf surface and available for volatilization [µg/m2]
Pest_quantity_leaf_Tissue, // Quantity of pesticide penetrated in the plant tissue		[µg/m2]
Pest_quantity_soil,		// Quantity of pesticide on the soil							[µg/m2]
ex_Pest,				// Saturated vapor pressure of the studied pesticide			[Pa]
K_photodegrad_Pest,		// photodegradation rate of the Pesticid						[s-1]
K_penetration_Pest,		// coeeficient for penetration of the pesticide in the plant tissue [s-1]

concentration_NH3,		// Concentration of NH3 at zref									[µg/m3]
NH3_canopy,				// NH3 concenrtation above the canopy (above Rb)				[µg/m3]
NH3_i,					// NH3 concentration inside the substomatal cavity				[µg/m3]
NH3_soil,				// NH3 concentration in the soil (in the wet soil layer)		[µg/m3]
NH3_leaf,				// NH3concentration near the leaves surface (below Rb)			[µg/m3]
FNH3_tot,				// Total NH3 flux at zref										[µg/m2/s]
FNH3_soil,				// Soil NH3 flux												[µg/m2/s]
FNH3_stom,				// NH3 flux to green leaves										[µg/m2/s]
FNH3_stom_inactive,		// NH3 flux to yellow leaves									[µg/m2/s]
FNH3_cut,				// Cuticular NH3 flux											[µg/m2/s]
FNH3_canopy,			// NH3 flux to canopy (stomata + cuticular)						[µg/m2/s]
GammaLeaf_NH3,			// NH3 emission potentiel for the leaf given by [NH4+]/[H+]		[-]
GammaSoil_NH3,			// NH3 emission potentiel for the soil given by [NH4+]/[H+]		[-]

displacement_height,	// Displacement height											[m]
z0_canopy,				// Canopy roughness length										[m]
uref,					// Wind speed at zref											[m/s]
uh,						// Wind speed at the top of the canopy							[m/s]
lf,						// characteristic length of the leaves							[m]
Keh,					// Eddy diffusivity at the topp of the canopy					[/s]
L_MO,					// Monin-Obukhov length											[m]
Ta,						// Air temperature at zref										[°C]
T_canopy,				// Air temperature above the canopy (above Rb_leaf)				[°C]
Zeta,					// (z-d)/L														[-]
x,						// Variable for stability calculations							[-]
PsiM,					// Integrated stability function for momentum					[-]
PsiH,					// Integrated stability function for heat						[-]
diff_Ra,				// Ra in neutral condition - Ra corrected for stability			[s/m]
glight,					// Response function of gstom to light							[-]
PAR,					// Photosynthetic Active Radiation								[µmol/m2/s]
gtemp,					// Response function of gstom to temperature					[-]				
gVPD,					// Response function of gstom to VPD							[-]
SWP,					// Soil Water Potential											[MPa]
gSWP,					// Response function of gstom to SWP							[-]
gSWC,					// Response function of gstom to SWC							[-]

r_thermal_wetlayer,		// resistance to heat transfer though the wet soil layer		[s/m]
r_thermal_drylayer,		// resistance to heat transfer though the dry soil layer (mulch)	[s/m]
r_soil_surface_H2O,		// resistance to the diffusion of water vapour through the dry soil layer (through the mulch)	[s/m]
depth_soil_Tsoilref,	// depth of the measure for the temperature in the soil 		[m]
depth_drysoil,			// depth of the soil dried out by the soil evaporation (depth of the mulch) [m]
depth_drysoil_threshold, // threshold for the increase of the Mulch						[m]

LEs, LEv,				// Latent heat flux for soil (s) and vegetation (v)				[W/m2]
Hs, Hv,					// Sensible heat flux for soil (s) and vegetation (v)			[W/m2]
LE, H,					// Total latent and sensible heat flux							[W/m2]
G,						// Conductive soil heat flux									[W/m2]

Rain,					// Amount of water from rainfall during the timestep			[mm]
Drain_Runoff,			// Amount of water drained or stream (run off) on the soil surface (i.e., not entering in the soil) 

Qmulch,					// quantite d eau du mulch     
Mulch_threshold,		// seuil de remplissage du mulch avant destruction (ratio : entre 0 et 1 - proche de 0 sol sableux a priori - proche de 1 = argile apriori)
Soil_Depth,				// Depth of the soil				[m]
WetLayer_Depth,			// thickness of the soil non-affected by the depth_drysoil (Soil_Depth - depth_drysoil) [m]  
DryWet_depth_change,	// variable intermediaire pour calculer le changement d epaisseur dela couche humide avec la pluie - cas d un remplissage par le bas (faible pluie ne detruisant pas le mulch) [m]
SWC_mulch,				// Soil Water content in the Mulch Layer (depth_drysoil)								[mm]
SWC_WetLayer,			// Soil Water Content in the wet soil layer (WetLayer_Depth)							[mm]
SWC_SoilTotal,			// Soil Water Content in the soil (WetLayer_Depth+DryLayer_Depth = in the Soil_Depth)	[mm]
SWC_Wiltingpoint,		// Soil Water Content under which the wilting point of the vegetation is reached		[mm]	
SWC_Dry,				// Soil Water Content under which the soil can not be more dried						[mm]
SWC_Fieldcapacity,		// Soil Water Content at the Field Capacity												[mm]
	teta_ini,			// soil moisture for the soil	: surface => Soil_Depth		** initial soil moisture = at the beginning of the simulation	[kg/kg]
	teta_tot,			//								: surface => Soil_Depth		** mean soil moisture for each time step						[kg/kg]
	teta_WetLayer,		//	soil moisture for the soil wet layer (calculated for the wet layer thickness) **										[kg/kg]
	teta_wiltingpoint,	// soil moisture under which the wilting point of the vegetation is reached													[kg/kg]	
	teta_minthreshold,	// Soil moisture under which the soil can not be more dried																	[kg/kg]
	teta_fieldcapacity, // Soil moisture at the Field Capacity

Rn_measured,			// Net radiation (model input)										[W.m-2]
Rn_calculated,			// Net radiation calculated											[W.m-2]
Rn,						// Net radiation (model input or calculated)						[W.m-2]
Rn_soil,				// Net radiation at the ground surface								[W.m-2]
Rn_veg,					// Net radiation available at the big leaf surface 					[W.m-2]
albedo_veg,				// result of the fration of albedoveg_green and Yellow				[-]
albedo_canopy,			// albedo of the canopy (soil + veg)								[-]

Rain_t3,				// Rainfall at the -3 timestep										[mm]
Rain_t2,				// Rainfall at the -2 timestep										[mm]
Rain_t1,				// Rainfall at the -1 timestep										[mm]
Cumul_Rain_t3,			// Cumulated rainfall during the last 3 timesteps					[mm]
evap_t3,				// Evaporation at the -3 timestep									[mm]
evap_t2,				// Evaporation at the -2 timestep									[mm]
evap_t1,				// Evaporation at the -1 timestep									[mm]
Cumul_evap_t3			// Cumulated evaporation during the last 3 timesteps				[mm]
;



/*-----------------Déclaration des paramètres--------------*/

	double
		ksoil_O3,					// Empirical coefficient to calculate Rsoil_O3					[-]
		Rsoilmin_O3,				// Empirical coefficient to calculate Rsoil_O3					[s/m]
		clay_content,				// Percentage of soil clay content for Rsoil_O3 calculation		[%]
		kcut_O3,					// Empirical coefficient to calculate Rcut_O3					[-]
		RH0_O3,						// Leaf air relative humidity threshold to calculate Rcut_O3	[%]
		Rcutmin_NH3,				// Minimal cuticular resistance to NH3 deposition				[s/m]
		kcut_NH3,					// Empirical coefficient to calculate Rcut_NH3					[-]
		Vaporization_enthalpy_Pest,	// Enthalpy of vaporization of the studied pesticide			[J.mol-1]
		ex_Pest_Tref,				// Saturated vapor pressure of the studied pesticide at Tref	[Pa]
		diffusivity_Pest,			// Molecular diffusivity in air of the studied pesticide		[m2/s]
		KH_Pest_Tref,				// Henry Constant at Tref										[-]
		Tref_Pest,					// Reference Temperature for KH and ex Values					[°C]
		Molar_mass_Pest,			// Molar mass of the studied Pesticide							[g.mol-1]
		Solub_Water_Pest,			// Pesticid Solubility in water									[g.L-1]
		Solub_Octanol_Pest,			// Pesticid Solubility in octanol								[-]
		logKow_Pest,				// Constant for adsoption of the Pesticid on the cuticle : adsorption properties on octan-water [-]
		K_photodegrad_Pest_ref,		// reference photodegradation constant of the Pesticid			[s-1]
		K_washoff_Pest,				// constant for washoff process of the pesticid after a rain	[mm-1]
		K_leaf_dissip_Pest,			// dissipation/degradation constant of the pesticid after penetration in the cuticle [s-1]
		K_stom_absorption_Pest,		// coefficient to estimate stomatal resistance to pesticide	: 1 = normal stomatal flux || +inf = no stomatal flux	[-]
		Lambda_Pest,				// coefficient related of the covering of the leaf by one layer of molecule of the pesticide. Usually equal to 1 kg/ha [µg.m-2]
		Pest_quantity_applied,		// Quantity of pesticide applied at the time of application		[µg/m2]
		Rcut_Pest,					// Cuticular resistance of the studied pesticide				[s/m]

		zref,						// Reference or measurement height								[m]
		wind_attenuation,			// Attenuation factor of wind speed								[-]
		z0_soil,					// Soil roughness length										[m]
		Coefficient_LAI_efficient,	// Coefficient determining efficient LAI (between 0 and 1)		[-]	
		albedo_veg_green, albedo_veg_yellow, albedo_soil, // albedo of the green vegetation, yellow vegetation and soil [-]										[-]
		emissivity_canopy,			// emissivity of the canopy										[-]

		gmax,						// Maximal vegetation stomatal conductance for H2O				[mmol_H2O/m2/s]
		gmin,						// Minimal stomatal conductance coefficient						[-]
		alphaglight,				// Coefficient of gstom response to light						[-]
		Tmin,						// Minimum temperature for stomata opening						[°C]
		Topt,						// Temperature for maximum stomata opening						[°C]
		Tmax,						// Maximum temperature for stomata closure						[°C]
		VPDmax,						// VPD below which stomata opening is maximum					[kPa]
		VPDmin,						// VPD above which stomata opening is minimum					[kPa]
		SWPmax,						// SWP below which stomata opening is minimum					[MPa]
		SWPmin,						// SWP above which stomata opening is maximum					[MPa]

		attenuationray,				// coefficient of Radiation attenuation through the canopy		[-] (often close to 0.7)
		
		beta_soil,					// coefficient for soil resistance in case of very dry soil


		Sc,							// Schmidt number for water							[-]
		Sc_O3,						// Schmidt number for ozone							[-]
		Sc_NH3,						// Schmidt number for NH3							[-]
		diffusivity_H2O,			// H2O molecular diffusivity in air at Tair and 1 atm	[m2/s] 
		diffusivity_O3,				// O3 molecular diffusivity in air at Tair and 1 atm	[m2/s] 
		diffusivity_NH3,			// NH3 molecular diffusivity in air at Tair and 1 atm	[m2/s] 


		Conduct_thermal_wetsoil,	// thermal conductivity for the wet soil layer									[W.m-1.K-1] 
		Conduct_thermal_drysoil,	// thermal conductivity for the dry soil layer									[W.m-1.K-1]
		depth_threshold_wetsoil,	// depth of the maximum thickness of the soil wet layer							[m]
		depth_threshold_drysoil,	// depth of the minimum thickness of the soil dry layer	to have a significant effect on soil water vapor transfers		[m]
		tortuosity_coeff,			// coefficent of the soil tortuosity											[-] (range 1.5-3.0 - standard value = 2.5)
		soil_porosity,				// porosity of the soil	- directly linked to the vomumetric soil water moisture [-] 
		rho_soil,					// dry bulk density	of the soil													[kg/m3]	
		Soil_Saturat_teta,			// volumetric soil moisture in saturated condition								[m3.m-3]
		Soil_Residual_teta,			// volumetric soil moisture for residual water in dry soil (for van Genuchten equation)		[m3.m-3]
		nVG,						// Van Genuchten n parameter													[-]
		alphaVG,					// Van Genuchten alpha parameter												[m-1]

		delta_time					// Timestep											******* WARNING ***********	[hour]
		;

	char Pest_Application_Date[150];
	char Pest_Application_Time[150];
	char Pest_Application_Date_Time[150];

	char buffer_time[256];

		/* intermediates variables or parameters */ 
	double
		Pprime1, Pprime2, Pprime3,		// 
		etha, coeffnu, Gamma1, Gamma2,
		intermediatetemperature, 
		critStab,					// threshold for iterative convergence of Energy balance closure focused on T° - [°C]
		critmaxStab,				// number of the maximum iteration for the iterative process (energy balance)	 [-] 
		critPprime,					
		critmaxPprime,
		critmaxRn					// acceptability threshold for the calculation of the net radiation in the iterative process (in case of net radiation calculation) - focused on T° [°C]  
		;

	int
		incrStab,
		incrmaxStab,
		incrPprime,
		incrmaxPprime,
		incrmaxRn
		; 


