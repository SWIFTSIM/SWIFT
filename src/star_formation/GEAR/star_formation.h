/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
 *******************************************************************************/
#ifndef SWIFT_GEAR_STAR_FORMATION_H
#define SWIFT_GEAR_STAR_FORMATION_H

/* Local includes */
#include "cosmology.h"
#include "error.h"
#include "hydro_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"
#include "star_formation_struct.h"
#include "random.h"
#include "cooling.h"
#include "engine.h"
#define PI 3.14159265358979323846


/**
 * @brief Compute the temperature of a #part based on the cooling function.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
INLINE static float get_temperature(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp) {

  /* Physical constants */
  const double m_H = phys_const->const_proton_mass; //d
  const double k_B = phys_const->const_boltzmann_k; //d

  /* Gas properties */
  const double T_transition = hydro_props->hydrogen_ionization_temperature; //d
  const double mu_neutral = hydro_props->mu_neutral; //d
  const double mu_ionised = hydro_props->mu_ionised; //d

  /* Particle temperature */
  const double u = hydro_get_physical_internal_energy(p, xp, cosmo); //d

  /* Temperature over mean molecular weight */
  const double T_over_mu = hydro_gamma_minus_one * u * m_H / k_B; //d

  /* Are we above or below the HII -> HI transition? */
  if (T_over_mu > (T_transition + 1.) / mu_ionised)
    return T_over_mu * mu_ionised;
  else if (T_over_mu < (T_transition - 1.) / mu_neutral)
    return T_over_mu * mu_neutral;
  else
    return T_transition;
}
/**
 * @brief Calculate if the gas has the potential of becoming
 * a star.
 *
 * @param starform the star formation law properties to use.
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 *
 */
INLINE static int star_formation_is_star_forming(  //eneleve le const sur xp REMIS
    const struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling) 
    { //begining
		//starform->cooling=cooling;
		const double G = phys_const->const_newton_G; //d
		const double kb = phys_const->const_boltzmann_k; //d
		const double mH = phys_const->const_proton_mass; //d
		float sigma=0.f; //f
		if(starform->with_sigma>0)
		{
			sigma=p->starform_data.sigma;
		}
		double const T=get_temperature(phys_const,hydro_props,us,cosmo,cooling,p,xp); ///////d
		xp->sf_data.temperature=T;
		//const double gamma=starform->EOS_gamma_effective;//supprime pour le moment, on le trouve d'une autre facon
		const int N=starform->Njeans;
		const double h=p->h; //smoothing length d
		const double mu=hydro_props->mu_neutral; /////juste de prendre mu neutre ???? d
		const double T0=starform->Max_temperature; ///////a changer d
		const double physical_density = hydro_get_physical_density(p, cosmo); //d
		const double rho_sfr=PI*(hydro_gamma*kb*T/mu/mH+sigma*sigma)/h/h/G/pow(N,2./3.); //d
		//message("T %e",hydro_gamma*kb*T/mu/mH);
		//message("sigma %e",sigma*sigma);
		//message("**************************333");
		//message("temperature %e",T);
		//message("t0 %e",T0);
		//message("densite %e",physical_density);
		//message("rhosfr %e",rho_sfr);
		//message("kb %e",kb);
		//message("mh %e",mH);
		//message("h %e",h);
		//message("G %e",G);
		//message("mu AAAA  %e",mu);
		if(T>T0) //the temperature is to high to form a star
		{
			//message("trop chaud");
			return 0;
		}
		else
		{
			if(physical_density>rho_sfr) //we can form a star
			{
				//message("temperature %e",T);
				//message("temperature_0 %e",T0);
				return 1;
			}
			else //the density is to low to form a star
			{
				//message("pas assez dense");
				return 0;
			}
		}
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart.
 *
 * Nothing to do here.
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR(
    const struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const double dt_star) 
    {
		//first check
		if (dt_star == 0.) {
			xp->sf_data.SFR = 0.;
		}
		else
		{
		const double G = phys_const->const_newton_G; //d
		const double c_star=starform->star_formation_rate; //d
		const double physical_density = hydro_get_physical_density(p, cosmo); //d
		if(physical_density!=0)
		{
			//double tff=sqrt(3*PI/32/G/physical_density);
			xp->sf_data.SFR=c_star*physical_density*sqrt(physical_density*32.f*G)/sqrt(3*PI);
			//xp->sf_data.SFR=1.-exp(-c_star*dt_star*sqrt(physical_density*32.f*G)/sqrt(3*PI)); //c_star*hydro_get_mass(p) / dt_star;
		}
		else
		{
			xp->sf_data.SFR=0.;
		}
		//xp->sf_data.proba=1.- exp(xp->sf_data.SFR * dt_star * (-1.)/physical_density);
		//if(xp->sf_data.SFR>1)
		//{
			//message("ATTENTION, SFR > 1 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
			//message("densite %e",physical_density);
		//}
		//xp->sf_data.SFR=xp->sf_data.SFR*hydro_get_mass(p) /dt_star; //pour des questions d'unites

	}
	return;
	}

/**
 * @brief Decides whether a particle should be converted into a
 * star or not.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 * @param e The #engine (for random numbers).
 * @param dt_star The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */ //eneleve le const sur xp REMIS RETIRE a nouveau (7.4.19)
INLINE static int star_formation_should_convert_to_star(
    const struct part* p, struct xpart* xp,
    const struct star_formation* starform, const struct engine* e,
    const double dt_star) {
	/* Calculate the propability of forming a star */
	//a voir si il faut garder ce terme dt/m
	
	const double prob = 1.- exp(xp->sf_data.SFR * dt_star * (-1.) / hydro_get_physical_density(p,e->cosmology));
	xp->sf_data.proba=1.- exp(xp->sf_data.SFR * 0.5*(e->dt_min+e->dt_max) * (-1.) / hydro_get_physical_density(p,e->cosmology));
	/* Get a unique random number between 0 and 1 for star formation */
	const double random_number =
    random_unit_interval(p->id, e->ti_current, random_number_star_formation); 
    if(random_number>prob) //no star
	{
		//message("pas eu de chance");
		return 0;
	}
	else //we form a new star particle
	{
		return 1;
	}
}

/**
 * @brief Update the SF properties of a particle that is not star forming.
 *
 * Nothing to do here.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param e The #engine.
 * @param starform The properties of the star formation model.
 * @param with_cosmology Are we running with cosmology switched on?
 */
INLINE static void star_formation_update_part_not_SFR(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct star_formation* starform, const int with_cosmology) {

  /* Check if it is the first time steps after star formation */
  if (xp->sf_data.SFR > 0.) {

    /* Record the current time as an indicator of when this particle was last
       star-forming. */
     // message("Pas de formation stellaire pour le moment !!!!!!KKKKKKK");
    if (with_cosmology) {
		//message("Cosmo %e", -(e->cosmology->a));
      xp->sf_data.SFR = -(e->cosmology->a);
    } else {
		//message("Pas de cosmo %e",-(e->time));
      xp->sf_data.SFR = -(e->time);
    }
  }
}

/**
 * @brief Copies the properties of the gas particle over to the
 * star particle.
 *
 * Nothing to do here.
 *
 * @param e The #engine
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sp the new created star particle with its properties.
 * @param starform the star formation law properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 */
 //juste pour tester pour le moment je copie l'autre
INLINE static void star_formation_copy_properties(
    const struct part* p, const struct xpart* xp, struct spart* sp,
    const struct engine* e, const struct star_formation* starform,
    const struct cosmology* cosmo, const int with_cosmology) {
	
  /* Store the current mass */
  sp->mass = hydro_get_mass(p);

  /* Store the current mass as the initial mass */
  sp->mass_init = hydro_get_mass(p);

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }

  /* Store the chemistry struct in the star particle */
  sp->chemistry_data = p->chemistry_data;

  /* Store the tracers data */
  sp->tracers_data = xp->tracers_data;

  /* Store the birth density in the star particle */
  sp->birth_density = hydro_get_physical_density(p, cosmo);
  /* Store the birth temperature*/
  sp->birth_temperature = get_temperature(starform->phys_const,starform->hydro_props,starform->us,cosmo,e->cooling_func,p,xp);
  //message("temperature %e",sp->birth_temperature);
  //message("densite %e",sp->birth_density);
}

/**
 * @brief initialization of the star formation law
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param starform the star formation law properties to initialize
 *
 */
INLINE static void starformation_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct hydro_props* hydro_props,
    struct star_formation* starform) {
		//noms des parametres a changer peut-etre plus tard.
		starform->Njeans = parser_get_param_int(
      parameter_file, "GEARStarFormation:jeans_number");
      starform->star_formation_rate = parser_get_param_double(
      parameter_file, "GEARStarFormation:star_formation_efficiency");
      starform->Max_temperature=parser_get_param_double(parameter_file, "GEARStarFormation:Max_temperature")*units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);;
      starform->hydro_props=hydro_props;
      starform->us=us;
      starform->phys_const=phys_const;
      starform->with_sigma = parser_get_param_int(
      parameter_file, "GEARStarFormation:with_turbulence_estimation");
      //starform->cooling=cooling;
      //quelques initialisations sans interet
      //xp->sf_data.star_created=false;
      
		}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {

  message("Star formation law is 'GEAR'");
}

/**
 * @brief Finishes the density calculation.
 *
 * @param p The particle to act upon
 * @param cd The global star_formation information.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void star_formation_end_density(
    struct part* restrict p, const struct star_formation* cd,
    const struct cosmology* cosmo) {
		//message("sigma_provisoire %e",p->starform_data.sigma);
		//message("density %e",hydro_get_physical_density(p,cosmo));
		p->starform_data.sigma=p->starform_data.sigma/hydro_get_physical_density(p,cosmo);
		}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #star_formation containing star_formation informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
star_formation_part_has_no_neighbours(struct part* restrict p,
                                      struct xpart* restrict xp,
                                      const struct star_formation* cd,
                                      const struct cosmology* cosmo) {
		p->starform_data.sigma=0.f;
										 
										  }

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state.
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param cosmo The current cosmological model.
 * @param data The global star_formation information used for this run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_part(const struct phys_const* restrict phys_const,
                               const struct unit_system* restrict us,
                               const struct cosmology* restrict cosmo,
                               const struct star_formation* data,
                               const struct part* restrict p,
                               struct xpart* restrict xp) {
								   }

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state.
 *
 * @param p Pointer to the particle data.
 * @param data The global star_formation information.
 */
__attribute__((always_inline)) INLINE static void star_formation_init_part(
    struct part* restrict p, const struct star_formation* data) {
		 p->starform_data.sigma=0.f;
		 //message("ATTTENTION : SIGMA INIT");
		}

#endif /* SWIFT_GEAR_STAR_FORMATION_H */
