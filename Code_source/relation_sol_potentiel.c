/* ---------Fonctions to convert teta in fi or fi in teta -----------------
   
   References :
   van Genuchten, M.Th., 1980. Soil Science Society of America Journal 44, 892–898. https://doi.org/10.2136/sssaj1980.03615995004400050002x
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "variables_parametres.h"


double relationtetafi(double humid,double potentiel, double rho_soil)
{
  double tetam ,teta1;
  double m ;
  double DappSoil;

  DappSoil = rho_soil / 1000.0;

  humid = humid * DappSoil ; // To convert weighted humidity to volumic humidity
  m = 1.0- 1.0 / nVG ; 
  
  tetam = pow(((Soil_Saturat_teta - Soil_Residual_teta)*DappSoil / ( humid - Soil_Residual_teta*DappSoil) ), 1.0/ m );
  teta1 = pow((tetam -1.0) , 1.0/nVG ) ;

  potentiel = - 1.0 / alphaVG * teta1 ; // Potential in m or 10E4Pa
  potentiel = potentiel / 100.0 ; // Potential in MPa

  return (potentiel);
} 


double relationfiteta(double potentiel, double humid, double rho_soil)
{
  double tetadiff, an , anpui , anpuim ;
  double m;
  double DappSoil;

  DappSoil = rho_soil / 1000.0;

  potentiel = potentiel * 100.0 ; // Conversion of potential from MPa to m or 10E4Pa

  m = 1.0-1.0 / nVG ;
  tetadiff = Soil_Saturat_teta - Soil_Residual_teta;
  an = fabs ( alphaVG * potentiel );
  anpui = pow ( an , nVG );
  anpuim = pow ( 1.0 + anpui , (-m));
  
  humid = Soil_Residual_teta + tetadiff * anpuim ;
  humid = humid / DappSoil; // Coversion of volumic humidity to weighted humidity: kg water / kg soil

  return(humid);
}