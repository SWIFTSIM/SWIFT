/*******************************************************************************
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_GEAR_STAR_FORMATION_IACT_H
#define SWIFT_GEAR_STAR_FORMATION_IACT_H

/**
 * @file GEAR/star_formation_iact.h
 * @brief Density computation
 */

/**
 * @brief do star_formation computation after the runner_iact_density (symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_star_formation(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) 
    {
	/*!the goal here is to estimate the local turbulence*/
	/*!value used to evaluate the SPH kernel*/
    float wi; 
    float wj;
    /*!evaluation of the SPH kernel*/
    kernel_eval(sqrt(r2)/hi,&wi);
    kernel_eval(sqrt(r2)/hj,&wj);
    /*!square of the velocity norm*/
	float norm_v2=pow(pi->v[0]-pj->v[0],2)+pow(pi->v[1]-pj->v[1],2)+pow(pi->v[2]-pj->v[2],2);
	/*!if(pi->id==(long long int)1)
	{
	message("**************"); 
	message("ancien %e",pi->starform_data.sigma);  }*/
	/*!estimation of local turbulence for pi and pj*/
	/*!if(pow(hi,-3)*norm_v2*wi*hydro_get_mass(pj)>pi->starform_data.sigma)
	{
		message("ATTENTION : %e",pi->starform_data.sigma);
		message("nouveau %e",pow(hi,-3)*norm_v2*wi*hydro_get_mass(pj));
	}*/
		/*!estimation of local turbulence for pi and pj*/
	/*!if(pow(hj,-3)*norm_v2*wj*hydro_get_mass(pi)>pj->starform_data.sigma)
	{
		message("ATTENTION : %e",pj->starform_data.sigma);
		message("nouveau %e",pow(hj,-3)*norm_v2*wj*hydro_get_mass(pi));
	}*/
	if(pi->starform_data.voisins==0 && pi->starform_data.sigma!=0.)
	{
		message("Debut %e",pi->starform_data.sigma);
	}
	if(pj->starform_data.voisins==0 && pj->starform_data.sigma!=0.)
	{
		message("Debut %e",pj->starform_data.sigma);
	}
	pi->starform_data.sigma=pi->starform_data.sigma+pow(hi,-3)*norm_v2*wi*hydro_get_mass(pj); //pi->starform_data.sigma+
	pj->starform_data.sigma=pj->starform_data.sigma+pow(hj,-3)*norm_v2*wj*hydro_get_mass(pi); //pj->starform_data.sigma+
	pi->starform_data.voisins+=1;
	pj->starform_data.voisins+=1;
	//message("W %e",wi);
	//message("norm v %e",norm_v2);
	//message("h-3 %e",pow(hi,-3));
	//message("mass %e",hydro_get_mass(pj));
	/*!if(pi->id==(long long int)1) {
	printf("id1 %lld \n",pi->id);
	printf("id2 %lld \n",pj->id);
	message("nouveau %e",pow(hi,-3)*norm_v2*wi*hydro_get_mass(pj));}*/
									  }

/**
 * @brief do star_formation computation after the runner_iact_density (non
 * symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_star_formation(float r2, const float *dx, float hi, float hj,
                                  struct part *restrict pi,
                                  const struct part *restrict pj, float a,
                                  float H) 
    {
	/*!the goal here is to estimate the local turbulence*/
	/*!value used to evaluate the SPH kernel*/
    float wi;
    /*!square of the velocity norm*/
	float norm_v2=pow(pi->v[0]-pj->v[0],2)+pow(pi->v[1]-pj->v[1],2)+pow(pi->v[2]-pj->v[2],2); 
	/*!evaluation of the SPH kernel*/
	kernel_eval(sqrt(r2)/hi,&wi);
	/*!estimation of local turbulence for pi*/
		/*!estimation of local turbulence for pi and pj*/
	if(pow(hi,-3)*norm_v2*wi*hydro_get_mass(pj)>pi->starform_data.sigma)
		if(pi->starform_data.voisins==0 && pi->starform_data.sigma!=0.)
	{
		message("ATTENTION Debut %e",pi->starform_data.sigma);
	}
	/*!{
		message("ATTENTION : %e ",pi->starform_data.sigma);
		message("nouveau %e",pow(hi,-3)*norm_v2*wi*hydro_get_mass(pj));
	}*/
	pi->starform_data.voisins+=1;
	pi->starform_data.sigma=pi->starform_data.sigma+pow(hi,-3)*norm_v2*wi*hydro_get_mass(pj); //pi->starform_data.sigma+
									  }

#endif /* SWIFT_GEAR_STAR_FORMATION_IACT_H */
