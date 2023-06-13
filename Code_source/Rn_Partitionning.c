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
