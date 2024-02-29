#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"
#include "./cooling.h"

#ifdef GRACKLE
#include <grackle.h>
#define ENDRUNVAL 91234

//
// 'mode' -- tells the routine what to do
//
//     0 == solve chemistry and assign new abundances
//     1 == calculate and return cooling time
//     2 == calculate and return temperature
//     3 == calculate and return pressure
//     4 == calculate and return gamma (only valid when GRACKLE_CHEMISTRY>0)
//

static grackle_field_data my_fields;
static chemistry_data *my_grackle_data;

double CallGrackle(double u_old, double rho, double dt, double *ne_guess, int target, int mode)
{
    gr_float returnval=0.0,gamma,cooling_time,temperature,pressure;
    
    if(All.ComovingIntegrationOn) All.GrackleUnits.a_value = All.cf_atime;

    /* Set grid dimension and size.  grid_start and grid_end are used to ignore ghost zones. */
    my_fields.grid_rank = 3;
    my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation
    my_fields.grid_dimension = malloc(GRACKLE_NPART * sizeof(int));
    my_fields.grid_start = malloc(GRACKLE_NPART * sizeof(int));
    my_fields.grid_end = malloc(GRACKLE_NPART * sizeof(int));
  
    for (i = 0;i < 3;i++) {
      my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
      my_fields.grid_start[i] = 0;
      my_fields.grid_end[i] = 0;
    }
    my_fields.grid_dimension[0] = GRACKLE_NPART;
    my_fields.grid_end[0] = GRACKLE_NPART - 1;

    /* load physical info for particle to be cooled */
    my_fields.density[0]       = rho; 
    my_fields.internal_energy[0] = u_old; //+ SphP[target].DtInternalEnergy*dt;
    my_fields.x_velocity[0]    = SphP[target].VelPred[0];
    my_fields.y_velocity[0]    = SphP[target].VelPred[1];
    my_fields.z_velocity[0]    = SphP[target].VelPred[2];
#ifdef METALS
    double Zmet = P[target].Metallicity[0];
#ifdef GALSF_DUST
    Zmet = (P[target].Metallicity[0]*(P[target].Mass-P[target].dust_Mass) + P[target].dust_Metallicity[0]*P[target].dust_Mass) / P[target].Mass;  // include dust metals when cooling, in lieu of a separate dust cooling module
#endif
    my_fields.metal_density[0] = my_fields.density[0] * DMAX(Zmet,1.e-8);
#else
    my_fields.metal_density[0] = my_fields.density[0] * 0.02; // if METALS not tracked, assume solar-ish metal cooling
#endif
    my_fields.specific_heating_rate[0] = SphP[target].DtInternalEnergy * All.HubbleParam * All.UnitEnergy_in_cgs / (All.UnitMass_in_g * All.UnitTime_in_s); 	// convert du/dt to cgs
    gamma         = GAMMA;
    
#if (GRACKLE_CHEMISTRY >  0) // non-tabular
    gr_float tiny = 1.0e-20;
    
    // Atomic
    my_fields.e_density[0]    = my_fields.density[0] * *ne_guess;
    my_fields.HI_density[0]    = my_fields.density[0] * SphP[target].grHI;  //initialized with HYDROGEN_MASSFRAC
    my_fields.HII_density[0]   = my_fields.density[0] * SphP[target].grHII;
    my_fields.HM_density[0]    = my_fields.density[0] * SphP[target].grHM;
    my_fields.HeI_density[0]   = my_fields.density[0] * SphP[target].grHeI;
    my_fields.HeII_density[0]  = my_fields.density[0] * SphP[target].grHeII;
    my_fields.HeIII_density[0] = my_fields.density[0] * SphP[target].grHeIII;
#endif
    
#if (GRACKLE_CHEMISTRY >= 1) // Atomic+(H2+H2I+H2II)

    my_fields.H2I_density[0]  = my_fields.density[0] * tiny;
    my_fields.H2II_density[0] = my_fields.density[0] * tiny;
    my_fields.DI_density[0]   = my_fields.density[0] * tiny;
    my_fields.DII_density[0]  = my_fields.density[0] * tiny;
    my_fields.HDI_density[0]  = my_fields.density[0] * tiny;
    
#if (GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
    my_fields.H2I_density[0]  = my_fields.density[0] * SphP[target].grH2I;
    my_fields.H2II_density[0] = my_fields.density[0] * SphP[target].grH2II;
#endif
    
#if (GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
    my_fields.DI_density[0]   = my_fields.density[0] * SphP[target].grDI;
    my_fields.DII_density[0]  = my_fields.density[0] * SphP[target].grDII;
    my_fields.HDI_density[0]  = my_fields.density[0] * SphP[target].grHDI;
#endif

#ifdef ARTIST
    my_fields.RT_heating_rate[0]  = 0.;  // 0 for now, until I understand what values it wants
    my_fields.RT_HI_ionization_rate[0]  = SphP[target].RT_gammaHI;  // HI photoionisation rate in s^-1
#endif
    
    switch(mode) {
        case 0:  //solve chemistry & update values
            if(solve_chemistry(&All.GrackleUnits,&my_fields,dt) == 0) {
                fprintf(stderr, "Error in Grackle solve_chemistry.\n");
                endrun(ENDRUNVAL);
            }
            
            // Assign variables back
            *ne_guess            = my_fields.e_density[0]    / my_fields.density[0];
            
            SphP[target].grHI    = my_fields.HI_density[0]    / my_fields.density[0];
            SphP[target].grHII   = my_fields.HII_density[0]   / my_fields.density[0];
            SphP[target].grHM    = my_fields.HM_density[0]    / my_fields.density[0];
            
            SphP[target].grHeI   = my_fields.HeI_density[0]   / my_fields.density[0];
            SphP[target].grHeII  = my_fields.HeII_density[0]  / my_fields.density[0];
            SphP[target].grHeIII = my_fields.HeIII_density[0] / my_fields.density[0];
            
#if (GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
            SphP[target].grH2I   = my_fields.H2I_density[0]   / my_fields.density[0];
            SphP[target].grH2II  = my_fields.H2II_density[0]  / my_fields.density[0];
#endif
            
#if (GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
            SphP[target].grDI    = my_fields.DI_density[0]    / my_fields.density[0];
            SphP[target].grDII   = my_fields.DII_density[0]   / my_fields.density[0];
            SphP[target].grHDI   = my_fields.HDI_density[0]   / my_fields.density[0];
#endif
            returnval = my_fields.internal_energy[0];
    //if(ThisTask==0) printf("G3: dt=%g lognH=%g dudt_adb=%g dudt_cgs=%g u_cgs=%g u_old=%g u_grackle=%g u_adb=%g\n",dt*All.GrackleUnits.time_units,log10(my_fields.density[0]*All.GrackleUnits.density_units/PROTONMASS),SphP[target].DtInternalEnergy,my_fields.specific_heating_rate[0],SphP[target].InternalEnergy*All.UnitEnergy_in_cgs/All.UnitMass_in_g,u_old,returnval,u_old+SphP[target].DtInternalEnergy*dt);
            break;
            
        case 1:  //cooling time
            if(calculate_cooling_time(&All.GrackleUnits, &my_fields, &cooling_time) == 0) {
                fprintf(stderr, "Error in Grackle calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = cooling_time;
            break;
        case 2:  //calculate temperature
            if(calculate_temperature(&All.GrackleUnits, &my_fields, &temperature) == 0) {
                fprintf(stderr, "Error in Grackle calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = temperature;
            break;
        case 3:  //calculate pressure
            if(calculate_pressure(&All.GrackleUnits, &my_fields, &pressure) == 0) {
                fprintf(stderr, "Error in Grackle calculate_pressure.\n");
                endrun(ENDRUNVAL);
            }
            returnval = pressure;
            break;
        case 4:  //calculate gamma
            if(calculate_gamma(&All.GrackleUnits, &my_fields, &gamma) == 0) {
                fprintf(stderr, "Error in Grackle calculate_gamma.\n");
                endrun(ENDRUNVAL);
            }
            returnval = gamma;
            break;
    } //end switch
    
#else // tabular
    
    switch(mode){
        case 0:  //solve chemistry & update values (table)
            if(solve_chemistry_table(&All.GrackleUnits,&my_fields,dt) == 0){
                fprintf(stderr, "Error in Grackle solve_chemistry_table.\n");
                endrun(ENDRUNVAL);
            }
            convert_u_to_temp(energy, rho, ne_guess, target); //need to update *ne_guess for tabular!!, this may be wrong
            returnval = energy;
            break;
        case 1:  //cooling time (table)
            if(calculate_cooling_time_table(&All.GrackleUnits, &my_fields, &cooling_time) == 0){
                fprintf(stderr, "Error in Grackle calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = cooling_time;
            break;
        case 2:  //calculate temperature (table)
            if(calculate_temperature_table(&All.GrackleUnits, &my_fields, &temperature) == 0){
                fprintf(stderr, "Error in Grackle calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = temperature;
            break;
        case 3:  //calculate pressure (table)
            if(calculate_pressure_table(&All.GrackleUnits, &my_fields &pressure) == 0){
                fprintf(stderr, "Error in Grackle calculate_pressure.\n");
                endrun(ENDRUNVAL);
            }
            returnval = pressure;
            break;
    } //end switch
    
#endif // GRACKLE_CHEMISTRY
    
    return returnval;
}




//Initialize Grackle
void InitGrackle(void)
{
    int i;

    // Enable output for Task 0
    // grackle_verbose = 0;
    // if(ThisTask == 0) grackle_verbose = 1;
   
    // Create a chemistry object for parameters and rate data.
    my_grackle_data = (chemistry_data *) malloc(sizeof(chemistry_data));
    if (set_default_chemistry_parameters(my_grackle_data) == 0) {
        fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
        exit(ENDRUNVAL);
    }

    // Set parameter values for chemistry & cooling

    // Flag to activate the grackle machinery:
    grackle_data->use_grackle            = 1;                   // grackle on (duh)
    // Flag to include radiative cooling and actually update the thermal energy during the
    // chemistry solver. If off, the chemistry species will still be updated. The most
    // common reason to set this to off is to iterate the chemistry network to an equilibrium state. Default: 1.
    grackle_data->with_radiative_cooling = 1;                   // cooling on
    // Flag to control which primordial chemistry network is used (set by Config file)
#ifdef GRACKLE_CHEMISTRY
    grackle_data->primordial_chemistry = GRACKLE_CHEMISTRY;
#else
    grackle_data->primordial_chemistry = 0;                     // fully tabulated cooling
#endif
    // Flag to enable H2 formation on dust grains, dust cooling, and dust-gas heat transfer follow Omukai (2000). This assumes that the dust to gas ratio scales with the metallicity. Default: 0.
    grackle_data->h2_on_dust             = 0;                   // dust cooling/chemistry on
    // Flag to enable metal cooling using the Cloudy tables. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
#ifdef METALS
    grackle_data->metal_cooling          = 1;                   // metal cooling on
#else
    grackle_data->metal_cooling          = 0;                   // metal cooling on
#endif
    // Flag to enable an effective CMB temperature floor. This is implemented by subtracting the value of the cooling rate at TCMB from the total cooling rate. Default: 1.
    grackle_data->cmb_temperature_floor  = 1;
    // Flag to enable a UV background. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
    grackle_data->UVbackground           = 1;                  // UV background on
    // Path to the data file containing the metal cooling and UV background tables:
    grackle_data->grackle_data_file      = All.GrackleDataFile; // data file
    // The ratio of specific heats for an ideal gas. A direct calculation for the molecular component is used if primordial_chemistry > 1. Default: 5/3.
    grackle_data->Gamma                  = GAMMA;              // our eos set in Config.sh
    // Flag to control which three-body H2 formation rate is used.
    grackle_data->three_body_rate        = 0;
    // Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel (2004). Default: 0.
    grackle_data->cie_cooling                      = 0;
    // Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004). Default: 0
    grackle_data->h2_optical_depth_approximation   = 0;
    // Flag to enable a spatially uniform heating term approximating photo-electric heating from dust from Tasker & Bryan (2008). Default: 0.
    grackle_data->photoelectric_heating            = 0;         // photo-electric on [but not adjusted to local background, beware!]
    grackle_data->photoelectric_heating_rate       = 8.5e-26;
    // Flag to enable Compton heating from an X-ray background following Madau & Efstathiou (1999). Default: 0.
    grackle_data->Compton_xray_heating   = 0;
    // Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field, in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
    grackle_data->LWbackground_intensity           = 0;
    // Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption
    //    (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000). Default: 0.
    grackle_data->LWbackground_sawtooth_suppression = 0;
    // volumetric heating rates is being provided  in the volumetric_heating_rate field of grackle_field_data 
    grackle_data->use_volumetric_heating_rate	= 0;
    // specific heating rates is being provided  in the specific_heating_rate field of grackle_field_data 
    grackle_data->use_specific_heating_rate	= 1;
#ifdef ARTIST
    // arrays of ionization and heating rates from radiative transfer solutions are being provided
    grackle_data->use_radiative_transfer		= 1;
    // must be enabled to couple the passed radiative transfer fields to the chemistry solver
    grackle_data->radiative_transfer_coupled_rate_solver	= 1;
    // enable intermediate stepping in applying radiative transfer fields to chemistry solver. 
    grackle_data->radiative_transfer_intermediate_step	= 1;
    // only use hydrogen ionization and heating rates from the radiative transfer solutions.
    grackle_data->radiative_transfer_hydrogen_only	= 1;
#else
    // arrays of ionization and heating rates from radiative transfer solutions are being provided
    grackle_data->use_radiative_transfer		= 0;
    // must be enabled to couple the passed radiative transfer fields to the chemistry solver
    grackle_data->radiative_transfer_coupled_rate_solver	= 0;
    // enable intermediate stepping in applying radiative transfer fields to chemistry solver. 
    grackle_data->radiative_transfer_intermediate_step	= 0;
    // only use hydrogen ionization and heating rates from the radiative transfer solutions.
    grackle_data->radiative_transfer_hydrogen_only	= 0;
#endif
    // Use Rahmati+13 self-shielding; 0=none, 1=HI only, 2=HI+HeI, 3=HI+HeI but set HeII rates to 0
    grackle_data->self_shielding_method	= 0;
#ifdef GRACKLE_SELFSHIELD
    grackle_data->self_shielding_method	= GRACKLE_SELFSHIELD;
#endif
    
    // Finally, initialize the chemistry object.
    if (initialize_chemistry_data(&All.GrackleUnits) == 0) {
        fprintf(stderr, "Error in initialize_chemistry_data.\n");
        exit(ENDRUNVAL);
    }
    
}

#endif  //GRACKLE
