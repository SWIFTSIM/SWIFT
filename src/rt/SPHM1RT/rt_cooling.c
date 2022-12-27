/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
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

/**
 * @file src/rt/SPHM1RT/rt_cooling.c
 * @brief SPHM1RT cooling functions
 */

/* Some standard headers. */
#include <cvode/cvode.h>
#include <cvode/cvode_direct.h> /* access to CVDls interface            */
#include <string.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

/* Local includes. */
#include "rt_cooling.h"
#include "rt_cooling_rates.h"
#include "rt_getters.h"
#include "rt_setters.h"

/**
 * @brief Main function for the thermochemistry step.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param dt The time-step of this particle.
 */
void rt_do_thermochemistry(struct part* restrict p, struct xpart* restrict xp,
                           struct rt_props* rt_props,
                           const struct cosmology* restrict cosmo,
                           const struct hydro_props* hydro_props,
                           const struct phys_const* restrict phys_const,
                           const struct unit_system* restrict us,
                           const double dt) {

  /* Nothing to do here? */
  if (rt_props->skip_thermochemistry == 1) return;
  if (dt == 0.0) return;

  rt_check_unphysical_elem_spec(p, rt_props);

  struct rt_part_data* rpd = &p->rt_data;

  struct RTUserData data; /* data for CVODE */

  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  const double conv_factor_internal_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  const double conv_factor_frad_to_cgs =
      conv_factor_internal_energy_to_cgs * units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);     
  const double conv_factor_opacity_from_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_MASS) /
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) /
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);

  // TK remark: here we only consider the on-the-spot approximation

  /**************************/
  /* INITIALIZATION         */
  /**************************/

  int useparams = rt_props->useparams;

  /* adopt on the spot approximation (OTSA) by default */
  /* TODO: currently there is no non-OTSA implementation; to do in the future */
  int onthespot = rt_props->onthespot;
  data.onthespot = onthespot;

  int coolingon = rt_props->coolingon;
  data.coolingon = coolingon;

  int fixphotondensity = rt_props->fixphotondensity;
  data.fixphotondensity = fixphotondensity;

  int smoothedRT = rt_props->smoothedRT;
  data.smoothedRT = smoothedRT;   

  double metal_mass_fraction[rt_chemistry_element_count];

  for (int elem = 0; elem < rt_chemistry_element_count; elem++) {
    metal_mass_fraction[elem] = (double)(rpd->tchem.metal_mass_fraction[elem]);
    data.metal_mass_fraction[elem] = metal_mass_fraction[elem];
  }

  const double X_H = metal_mass_fraction[rt_chemistry_element_H];

  const double m_H_cgs = phys_const->const_proton_mass *
                         units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double proton_mass_cgs_inv = 1.0 / m_H_cgs;
  data.m_H_cgs = m_H_cgs;

  const double k_B_cgs = phys_const->const_boltzmann_k *
                         units_cgs_conversion_factor(us, UNIT_CONV_ENERGY) /
                         units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  data.k_B_cgs = k_B_cgs;

  const double cred_phys = rt_get_physical_cred(p, cosmo->a);
  const double cred_cgs =
      cred_phys * units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);
  data.cred_cgs = cred_cgs;

  /* Get particle density [ and convert to g * cm^-3] */
  const double rho = hydro_get_physical_density(p, cosmo);
  double rho_cgs = rho * units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  data.rho_cgs = rho_cgs;

  /* Hydrogen number density (X_H * rho / m_p) [cm^-3] */
  const double n_H_cgs = X_H * rho_cgs * proton_mass_cgs_inv;
  data.n_H_cgs = n_H_cgs;

  /* Current energy (in internal units) */
  float urad[RT_NGROUPS];
  rt_get_physical_urad_multifrequency(p, cosmo, urad);

  /* Current flux (in internal units)*/
  float frad[RT_NGROUPS][3];
  rt_get_physical_frad_multifrequency(p, cosmo, frad);


  /* need to convert to cgs */
  double ngamma_cgs[3];
  double fgamma_cgs[3][3];

  /* for now, the 0th bin for urad is 0-HI, so we ignore it */
  for (int g = 0; g < 3; g++) {
    ngamma_cgs[g] =
        (double)(rho_cgs * urad[g + 1] * conv_factor_internal_energy_to_cgs /
                 rt_props->ionizing_photon_energy_cgs[g]);
    data.ngamma_cgs[g] = ngamma_cgs[g];
    fgamma_cgs[g][0] =
        (double)(rho_cgs * frad[g + 1][0] * conv_factor_frad_to_cgs /
        rt_props->ionizing_photon_energy_cgs[g]);
    data.fgamma_cgs[g][0] = fgamma_cgs[g][0];
    fgamma_cgs[g][1] =
        (double)(rho_cgs * frad[g + 1][1] * conv_factor_frad_to_cgs /
        rt_props->ionizing_photon_energy_cgs[g]);
    data.fgamma_cgs[g][1] = fgamma_cgs[g][1];
    fgamma_cgs[g][2] =
        (double)(rho_cgs * frad[g + 1][2] * conv_factor_frad_to_cgs /
        rt_props->ionizing_photon_energy_cgs[g]);
    data.fgamma_cgs[g][2] = fgamma_cgs[g][2];
  }



  /* overwrite the photon density if we choose to fix it */
  for (int i = 0; i < 3; i++) {
    if ((rt_props->Fgamma_fixed_cgs[i] > 0.0) &&
        (rt_props->fixphotondensity == 1)) {
      ngamma_cgs[i] = rt_props->Fgamma_fixed_cgs[i] / cred_cgs;
      data.ngamma_cgs[i] = ngamma_cgs[i];
      urad[i + 1] = (float)(data.ngamma_cgs[i] / rho_cgs /
                            conv_factor_internal_energy_to_cgs *
                            rt_props->ionizing_photon_energy_cgs[i]);
    }
  }


  float uradinj[RT_NGROUPS], uradratepro[RT_NGROUPS];
  /* Current energy injection (not rate; in internal units) */
  rt_get_physical_urad_injection(p, cosmo, uradinj);
  /* Current energy propagation rate (in internal units) */
  rt_get_physical_urad_propagation_rate(p, cosmo, uradratepro); 

  
  float fradinj[RT_NGROUPS][3], fradratepro[RT_NGROUPS][3];
  /* Current flux injection rate (in internal units) */
  rt_get_physical_frad_injection(p, cosmo, fradinj);
  /* Current flux injection rate (in internal units) */
  rt_get_physical_frad_propagation_rate(p, cosmo, fradratepro);

  double ngamma_inject_rate_cgs[3], fgamma_inject_rate_cgs[3][3];
  /* for now, the 0th bin for urad is 0-HI, so we ignore it */
  if (rt_props->smoothedRT == 1) {
    for (int g = 0; g < 3; g++) {
      ngamma_inject_rate_cgs[g] =
          (double)(rho_cgs * (uradinj[g + 1]/dt+uradratepro[g+1]) * conv_factor_internal_energy_to_cgs /
                  rt_props->ionizing_photon_energy_cgs[g] / 
                  units_cgs_conversion_factor(us, UNIT_CONV_TIME));
      data.ngamma_inject_rate_cgs[g] = ngamma_inject_rate_cgs[g];
    }
    for (int g = 0; g < 3; g++) {    
      fgamma_inject_rate_cgs[g][0] =
          (double)(rho_cgs * (fradinj[g + 1][0]/dt+fradratepro[g + 1][0]) * conv_factor_frad_to_cgs /
                  rt_props->ionizing_photon_energy_cgs[g] /
                  units_cgs_conversion_factor(us, UNIT_CONV_TIME));
      data.fgamma_inject_rate_cgs[g][0] = fgamma_inject_rate_cgs[g][0];
      fgamma_inject_rate_cgs[g][1] =
          (double)(rho_cgs * (fradinj[g + 1][1]/dt+fradratepro[g + 1][1]) * conv_factor_frad_to_cgs /
                  rt_props->ionizing_photon_energy_cgs[g] /
                  units_cgs_conversion_factor(us, UNIT_CONV_TIME));
      data.fgamma_inject_rate_cgs[g][1] = fgamma_inject_rate_cgs[g][1];
      fgamma_inject_rate_cgs[g][2] =
          (double)(rho_cgs * (fradinj[g + 1][2]/dt+fradratepro[g + 1][2]) * conv_factor_frad_to_cgs /
                  rt_props->ionizing_photon_energy_cgs[g] /
                  units_cgs_conversion_factor(us, UNIT_CONV_TIME));
      data.fgamma_inject_rate_cgs[g][2] = fgamma_inject_rate_cgs[g][2];
    }
  } else {
    for (int g = 0; g < 3; g++) {
      ngamma_inject_rate_cgs[g] = 0.0;
      data.ngamma_inject_rate_cgs[g] = ngamma_inject_rate_cgs[g];
    }    
    for (int g = 0; g < 3; g++) {    
      fgamma_inject_rate_cgs[g][0] = 0.0;
      data.fgamma_inject_rate_cgs[g][0] = fgamma_inject_rate_cgs[g][0];
      fgamma_inject_rate_cgs[g][1] = 0.0;
      data.fgamma_inject_rate_cgs[g][1] = fgamma_inject_rate_cgs[g][1];
      fgamma_inject_rate_cgs[g][2] = 0.0;
      data.fgamma_inject_rate_cgs[g][2] = fgamma_inject_rate_cgs[g][2];
    }
  }


  double abundances[rt_species_count];

  for (int spec = 0; spec < rt_species_count; spec++) {
    abundances[spec] = (double)(rpd->tchem.abundances[spec]);
    data.abundances[spec] = abundances[spec];
  }

  const double u = hydro_get_physical_internal_energy(p, xp, cosmo);

  double u_cgs = u * conv_factor_internal_energy_to_cgs;

  double T_cgs = rt_convert_u_to_temp(k_B_cgs, m_H_cgs, X_H, u_cgs, abundances);

  double T_min_cgs = hydro_props->minimal_temperature;

  double u_min_cgs =
      rt_convert_temp_to_u(k_B_cgs, m_H_cgs, T_min_cgs, X_H, abundances);

  u_cgs = fmax(u_cgs, u_min_cgs);

  data.u_cgs = u_cgs;

  data.u_min_cgs = u_min_cgs;

  /**************************/
  /* GET RATE COEFFICIENTS  */
  /**************************/
  int aindex[3];

  rt_get_index_to_species(aindex);
  for (int i = 0; i < 3; i++) {
    data.aindex[i] = aindex[i];
  }

  double alphalist[rt_species_count], betalist[rt_species_count],
      Gammalist[rt_species_count], sigmalist[3][3], epsilonlist[3][3];

  if (useparams == 1) {
    betalist[rt_sp_elec] = 0.0;
    betalist[rt_sp_HI] = rt_props->beta_cgs_H;
    betalist[rt_sp_HII] = 0.0;
    betalist[rt_sp_HeI] = 0.0;
    betalist[rt_sp_HeII] = 0.0;
    betalist[rt_sp_HeIII] = 0.0;
    alphalist[rt_sp_elec] = 0.0;
    alphalist[rt_sp_HI] = 0.0;
    alphalist[rt_sp_HeI] = 0.0;
    if (onthespot == 1) {
      alphalist[rt_sp_HII] = rt_props->alphaB_cgs_H;
      alphalist[rt_sp_HeII] = 0.0;
      alphalist[rt_sp_HeIII] = 0.0;
    } else {
      alphalist[rt_sp_HII] = rt_props->alphaA_cgs_H;
      alphalist[rt_sp_HeII] = 0.0;
      alphalist[rt_sp_HeIII] = 0.0;
    }
    sigmalist[0][0] = rt_props->sigma_cross_cgs_H[0];
    sigmalist[1][0] = rt_props->sigma_cross_cgs_H[1];
    sigmalist[2][0] = rt_props->sigma_cross_cgs_H[2];
    sigmalist[0][1] = 0.0;
    sigmalist[1][1] = 0.0;
    sigmalist[2][1] = 0.0;
    sigmalist[0][2] = 0.0;
    sigmalist[1][2] = 0.0;
    sigmalist[2][2] = 0.0;
    data.alphaA_cgs_H = rt_props->alphaA_cgs_H;
    data.alphaB_cgs_H = rt_props->alphaB_cgs_H;
    data.beta_cgs_H = rt_props->beta_cgs_H;
    data.sigma_cross_cgs_H[0] = rt_props->sigma_cross_cgs_H[0];
    data.sigma_cross_cgs_H[1] = rt_props->sigma_cross_cgs_H[1];
    data.sigma_cross_cgs_H[2] = rt_props->sigma_cross_cgs_H[2];
    for (int spec = 0; spec < rt_species_count; spec++) {
      Gammalist[spec] = 0.0;
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        epsilonlist[i][j] = 0.0;
      }
    }
  } else {
    rt_compute_rate_coefficients(T_cgs, onthespot, alphalist, betalist,
                                 Gammalist, sigmalist, epsilonlist);
  }

  data.useparams = rt_props->useparams;

  /**************************/
  /* SOLVING RATE EQUTAION  */
  /**************************/

  /* Try explicit solution */

  double new_abundances[rt_species_count], finish_abundances[rt_species_count],
      max_relative_change, new_ngamma_cgs[3], new_fgamma_cgs[3][3], u_new_cgs;

  max_relative_change = 0.0;
  /* compute net changes and cooling and heating for explicit solution */
  rt_compute_explicit_thermochemistry_solution(
      n_H_cgs, cred_cgs, dt_cgs, rho_cgs, u_cgs, u_min_cgs, abundances,
      ngamma_cgs, ngamma_inject_rate_cgs, fgamma_cgs, fgamma_inject_rate_cgs,
      alphalist, betalist, Gammalist, sigmalist, epsilonlist,
      aindex, &u_new_cgs, new_abundances, new_ngamma_cgs, new_fgamma_cgs,
      &max_relative_change);

  /* check whether xHI bigger than one */
  int errorHI = 0;
  if (new_abundances[rt_sp_HI] > 1.01) {
    errorHI = 1;
  } else {
    rt_enforce_constraint_equations(new_abundances, metal_mass_fraction,
                                    finish_abundances);
  }

  if ((max_relative_change > rt_props->explicitRelTolerance) ||
      (errorHI == 1)) {

    /**************************************
     * Explicit solution is insufficient. *
     * Use implicit solver.               *
     **************************************/
    realtype reltol, t;
    N_Vector abstol_vector, y;

    int maxsteps = 100000;
    int network_size, icount = 0;
    /* 3 for species;   */
    network_size = 3;
    /* 1 for thermal energy; */
    if (coolingon == 1) {
      network_size += 1;
    }
    /* 3 for radiation bins */
    if (fixphotondensity == 0) {
      network_size += 3;
      /* 9 for radiation flux */
      network_size += 9;
    }
   


    y = N_VNew_Serial(network_size);
    abstol_vector = N_VNew_Serial(network_size);
    for (int i = 0; i < 3; i++) {
      NV_Ith_S(y, icount) = (realtype)data.abundances[aindex[i]];
      NV_Ith_S(abstol_vector, icount) = (realtype)(rt_props->absoluteTolerance);
      icount += 1;
    }
    if (coolingon == 1) {
      NV_Ith_S(y, icount) = (realtype)u_cgs;
      NV_Ith_S(abstol_vector, icount) = (realtype)rt_props->absoluteTolerance;
      icount += 1;
    }
    if (fixphotondensity == 0) {
      for (int i = 0; i < 3; i++) {
        NV_Ith_S(y, icount) = (realtype)data.ngamma_cgs[i];
        NV_Ith_S(abstol_vector, icount) = (realtype)rt_props->absoluteTolerance;
        icount += 1;
      }
    }
    if (fixphotondensity == 0) {
      for (int i = 0; i < 3; i++) {
        NV_Ith_S(y, icount) = (realtype)data.fgamma_cgs[i][0];
        NV_Ith_S(abstol_vector, icount) = (realtype)rt_props->absoluteTolerance;
        icount += 1;
        NV_Ith_S(y, icount) = (realtype)data.fgamma_cgs[i][1];
        NV_Ith_S(abstol_vector, icount) = (realtype)rt_props->absoluteTolerance;
        icount += 1;
        NV_Ith_S(y, icount) = (realtype)data.fgamma_cgs[i][2];
        NV_Ith_S(abstol_vector, icount) = (realtype)rt_props->absoluteTolerance;
        icount += 1;
      }
    }

    data.network_size = network_size; 
    /* check if the number of inputs agrees with the number of equations */
    if (icount != network_size) {
      error("Error: at beginning: icount does not agree with network_size %i, %i", icount,network_size);
    }

    /* Set up the solver */
    /* Set the tolerances*/
    reltol = (realtype)rt_props->relativeTolerance;

    /* Use CVodeCreate to create the solver
     * memory and specify the Backward Differentiation
     * Formula. Note that CVODE now uses Newton iteration
     * iteration by default, so no need to specify this. */
    void* cvode_mem;
    cvode_mem = CVodeCreate(CV_BDF);

    /* Set the user data for CVode */
    CVodeSetUserData(cvode_mem, &data);

    /* Use CVodeSetMaxNumSteps to set the maximum number
     * of steps CVode takes. */
    CVodeSetMaxNumSteps(cvode_mem, maxsteps);

    /* Use CVodeInit to initialise the integrator
     * memory and specify the right hand side
     * function in y' = f(t,y) (i.e. the rate
     * equations), the initial time 0.0 and the
     * initial conditions, in y. */
    CVodeInit(cvode_mem, rt_frateeq, 0.0f, y);

    /* Use CVodeSVtolerances to specify the scalar
     * relative and absolute tolerances. */
    CVodeSVtolerances(cvode_mem, reltol, abstol_vector);

    /* Create a dense SUNMatrix to use in the
     * linear solver. */
    SUNMatrix A_sun;

    A_sun = SUNDenseMatrix(network_size, network_size);

    /* Create a denst SUNLinearSolver object
     * to use in CVode. */
    SUNLinearSolver LS_sun;
    LS_sun = SUNDenseLinearSolver(y, A_sun);

    /* Attach the matrix and linear
     * solver to CVode. */
    CVDlsSetLinearSolver(cvode_mem, LS_sun, A_sun);

    /* Specify the maximum number of convergence
     * test failures. */
    CVodeSetMaxConvFails(cvode_mem, 5000);

    /* Call CVode() to integrate the chemistry. */
    CVode(cvode_mem, (realtype)dt_cgs, y, &t, CV_NORMAL);

    /* Write the output abundances to the gas cell
     * Note that species not included in the reduced
     * network are kept constant in the GasVars struct. */
    icount = 0;
    for (int i = 0; i < 3; i++) {
      new_abundances[aindex[i]] = (double)NV_Ith_S(y, icount);
      icount += 1;
    }
    if (coolingon == 1) {
      u_cgs = (double)NV_Ith_S(y, icount);
      icount += 1;
    }

    if (fixphotondensity == 0) {
      for (int i = 0; i < 3; i++) {
        new_ngamma_cgs[i] = (double)NV_Ith_S(y, icount);
        icount += 1;
      }
    }

    if (fixphotondensity == 0) { 
      for (int i = 0; i < 3; i++) {
        new_fgamma_cgs[i][0] = (double)NV_Ith_S(y, icount);
        icount += 1;
        new_fgamma_cgs[i][1] = (double)NV_Ith_S(y, icount);
        icount += 1;
        new_fgamma_cgs[i][2] = (double)NV_Ith_S(y, icount);
        icount += 1;
      }
    }

    /* check if the number of outputs agrees with the number of equations */
    if (icount != network_size) {
      error("Error: at out: icount does not agree with network_size %i, %i", icount,network_size);
    }

    if (new_abundances[rt_sp_HI] > 1.01)
      error("HI fraction bigger than one after the CVODE solver");
    rt_enforce_constraint_equations(new_abundances, metal_mass_fraction,
                                    finish_abundances);
    SUNLinSolFree(LS_sun);
    SUNMatDestroy(A_sun);
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(abstol_vector);
    CVodeFree(&cvode_mem);
  }



  for (int spec = 0; spec < rt_species_count; spec++) {
    if (finish_abundances[spec] > 0.f) {
      if (finish_abundances[spec] < FLT_MAX) {
        rpd->tchem.abundances[spec] = (float)(finish_abundances[spec]);
      } else {
        error("finish_abundances larger than FLT_MAX");
      }
    } else {
      rpd->tchem.abundances[spec] = 0.f;
    }
  }
  if (coolingon == 1) {
    float u_new = 0.0f;
    if (u_new_cgs / conv_factor_internal_energy_to_cgs > 0.f) {
      if (u_new_cgs / conv_factor_internal_energy_to_cgs < FLT_MAX) {
        u_new = (float)(u_new_cgs / conv_factor_internal_energy_to_cgs);
      }
    }
    hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  }

  /* set radiation energy */
  float urad_new[RT_NGROUPS];
  urad_new[0] = 0.f;
  if (fixphotondensity == 0) {
    for (int i = 0; i < 3; i++) {
      urad_new[i + 1] = 0.f;
      if (new_ngamma_cgs[i] / rho_cgs / conv_factor_internal_energy_to_cgs *
              rt_props->ionizing_photon_energy_cgs[i] >
          0.f) {
        if (new_ngamma_cgs[i] / rho_cgs / conv_factor_internal_energy_to_cgs *
                rt_props->ionizing_photon_energy_cgs[i] <
            FLT_MAX) {
          urad_new[i + 1] = (float)(new_ngamma_cgs[i] / rho_cgs /
                                    conv_factor_internal_energy_to_cgs *
                                    rt_props->ionizing_photon_energy_cgs[i]);
        }
      }
    }
  } else {
    for (int i = 0; i < 3; i++) {
      urad_new[i + 1] = urad[i + 1];
    }
  }


  /* set radiation flux */
  if (fixphotondensity == 0) { 
    float frad_new[RT_NGROUPS][3];
    float frad_new_single[3];
    float fradinjmag;
    frad_new[0][0] = 0.f;
    frad_new[0][1] = 0.f;
    frad_new[0][2] = 0.f;      
    const char loc[30] = "rt_do_thermochemistry";
    for (int i = 0; i < 3; i++) {
      frad_new[i + 1][0] = 0.f;
      frad_new[i + 1][1] = 0.f;
      frad_new[i + 1][2] = 0.f; 
      if (new_ngamma_cgs[i] / rho_cgs / conv_factor_internal_energy_to_cgs *
              rt_props->ionizing_photon_energy_cgs[i] >
          0.f) {
        if (new_ngamma_cgs[i] / rho_cgs / conv_factor_internal_energy_to_cgs *
                rt_props->ionizing_photon_energy_cgs[i] <
            FLT_MAX) {
          fradinjmag = sqrtf(fradinj[i + 1][0]*fradinj[i + 1][0]
                    +fradinj[i + 1][1]*fradinj[i + 1][1]
                    +fradinj[i + 1][2]*fradinj[i + 1][2]);        
          if (fradinjmag > 0.f) {
            frad_new[i + 1][0] = urad_new[i + 1] * cred_phys * fradinj[i + 1][0] / fradinjmag; 
            frad_new[i + 1][1] = urad_new[i + 1] * cred_phys * fradinj[i + 1][1] / fradinjmag; 
            frad_new[i + 1][2] = urad_new[i + 1] * cred_phys * fradinj[i + 1][2] / fradinjmag; 
          } else {
            frad_new[i + 1][0] = (float)(new_fgamma_cgs[i][0] / rho_cgs /
                                  conv_factor_frad_to_cgs *
                                  rt_props->ionizing_photon_energy_cgs[i]);
            frad_new[i + 1][1] = (float)(new_fgamma_cgs[i][1] / rho_cgs /
                                  conv_factor_frad_to_cgs *
                                  rt_props->ionizing_photon_energy_cgs[i]);
            frad_new[i + 1][2] = (float)(new_fgamma_cgs[i][2] / rho_cgs /
                                  conv_factor_frad_to_cgs *
                                  rt_props->ionizing_photon_energy_cgs[i]);
          }
        }
      }
      frad_new_single[0] = frad_new[i+1][0]; 
      frad_new_single[1] = frad_new[i+1][1]; 
      frad_new_single[2] = frad_new[i+1][2];  
      rt_check_unphysical_state(&urad_new[i+1], frad_new_single,
                                0.0, cred_phys, loc);
      frad_new[i+1][0] = frad_new_single[0];
      frad_new[i+1][1] = frad_new_single[1];
      frad_new[i+1][2] = frad_new_single[2];      
    }
    rt_set_physical_radiation_flux_multifrequency(p, cosmo, frad_new);
  }
  rt_set_physical_urad_multifrequency(p, cosmo, urad_new);


  /* chi is in physical unit (L^2/M) */
  float chi_new[RT_NGROUPS];
  for (int i = 0; i < RT_NGROUPS; i++) {
    chi_new[i] = 0.0f;
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (finish_abundances[aindex[j]] * n_H_cgs / rho_cgs * sigmalist[i][j] *
              conv_factor_opacity_from_cgs >
          0.f) {
        if (finish_abundances[aindex[j]] * n_H_cgs / rho_cgs *
                sigmalist[i][j] * conv_factor_opacity_from_cgs <
            FLT_MAX) {
          chi_new[i + 1] +=
              (float)(finish_abundances[aindex[j]] * n_H_cgs / rho_cgs *
                      sigmalist[i][j] * conv_factor_opacity_from_cgs);
        }
      }
    }
  }
  rt_set_physical_radiation_opacity(p, cosmo, chi_new);

  rt_check_unphysical_elem_spec(p, rt_props);


}

/**
 * @brief Do the thermochemistry on a particle.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param dt The time-step of this particle.
 */
void rt_tchem(struct part* restrict p, struct xpart* restrict xp,
              struct rt_props* rt_props, const struct cosmology* restrict cosmo,
              const struct hydro_props* hydro_props,
              const struct phys_const* restrict phys_const,
              const struct unit_system* restrict us, const double dt) {
  rt_do_thermochemistry(p, xp, rt_props, cosmo, hydro_props, phys_const, us,
                        dt);

  /* reset injection rate after thermochemistry*/
  struct rt_part_data* rpd = &p->rt_data;
  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->dconserved_inj[g].urad = 0.0f;
    rpd->dconserved_inj[g].frad[0] = 0.0f;
    rpd->dconserved_inj[g].frad[1] = 0.0f;
    rpd->dconserved_inj[g].frad[2] = 0.0f;
  }
}
