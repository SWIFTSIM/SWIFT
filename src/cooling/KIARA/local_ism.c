
/*
Code to calculate the local isrf incident upon each gas particle.
This is used in Grackle's calculation of the dust growth and evolution.
Author: Ewan Jones
Date: December 2021
*/

//* Function to find the number of, and indices, of all neighbouring gas particles.
void cooling_gas_neighbours(int target, int *num_neighbours, int *ind_neighbours, double smoothingFactor)
{
    /* Initialise required variables. */
    //     Constants
    int startNode=All.MaxPart, ignore=0;
    //     Softening length 
    double h = smoothingFactor * P[target].Hsml;
    if(h<=0){
        h=All.SofteningTable[0];
    }

    // Obtain the indeces and number of neighbouring gas particles.
    *num_neighbours = ngb_treefind_pairs_threads(P[target].Pos, h, -1, &startNode, 0, &ignore, &ignore, &ignore, ind_neighbours);
}

//* Function to find the number of, and indeces, of all neighbouring non-gas particles.
void cooling_nongas_neighbours(int target, int *num_neighbours, int *ind_neighbours, double smoothingFactor)
{
    /* Initialise required variables. */
    //     Constants
    int startNode=All.MaxPart, ignore=0;
    //     Softening length 
    double h = smoothingFactor * P[target].Hsml;
    if(h<=0){
        h=All.SofteningTable[0];
    }

    // Obtain the indeces and number of neighbouring gas particles.
    *num_neighbours = ngb_treefind_variable_threads_nongas(P[target].Pos, h, -1, &startNode, 0, &ignore, &ignore, &ignore,
                                                             ind_neighbours);
}

//* Function to calculate the local (average) sSFR for a given gas particle *//
void cooling_calc_local_ssfr(int target, double smoothingFactor)
{
    // Variables to hold sums within kernel.
    MyFloat total_sfr=0., total_stellar_mass=0.;

    // Find all neighbouring gas particles.
    int num_ngb_gas, ind_ngb_gas[NumPart];
    gas_neighbours(target, &num_ngb_gas, ind_ngb_gas, smoothingFactor);

    // Find all neighboring nongas particles;
    int num_ngb_nongas, ind_ngb_nongas[NumPart];
    nongas_neighbours(target, &num_ngb_nongas, ind_ngb_nongas, smoothingFactor);
    
    // Loop through the neighbouring gas particles and add their SFR to total.
    int n, n_ind;
    for (n=0; n < num_ngb_gas; n++){
        // Get index of neighbouring particle.
        n_ind = ind_ngb_gas[n];
        // Add its sfr to the total
        total_sfr += SphP[n_ind].Sfr;
    }

    // Loop through the neighbouring non-gas particles, if it is a star add stellar mass to total.
    for (n=0; n < num_ngb_nongas; n++){
        // Get index of neighbouring particle.
        n_ind = ind_ngb_nongas[n];
        // If its a star particle, add its stellar mass to the total.
        if (P[n_ind].Type == 4){
            total_stellar_mass += P[n_ind].Mass;
        }
    }

    // If the total stellar mass is zero then estimate it using the SFR and simulation timestep.
    if ( total_stellar_mass == 0 ) {
        total_stellar_mass = total_sfr * All.TimeStep;
    } 
        
    // Caclulate the local sSFR around the particle.
    if ( total_stellar_mass == 0 ) {
        P[target].ssfr_local = 0;
    } else {
        P[target].ssfr_local = total_sfr / total_stellar_mass;
    }
    
    //! debug
    if (P[target].ssfr_local != P[target].ssfr_local) {
        fprintf(stdout, "P[target].ssfr_local = %g, total_stellar_mass = %g, total_sfr = %g, All.TimeStep = %g \n", P[target].ssfr_local, total_stellar_mass, total_sfr, All.TimeStep);
    }

    // Save the timestep that this was calculated.
    P[target].ssfr_local_lastCalcTi = All.NumCurrentTiStep;
}

//* Function to smooth the local sSFR for a given gas particle, this helps alleviate non-physical boundaries *//
double cooling_smooth_local_ssfr(int target, double smoothingFactor)
{
    // INPUTS:
    //   target (int) --> array index of the particle, same convention as CallGrackle().
    // OUTPUTS:
    //   ssfr (double) --> will smooth ssfr with neighbouring particles and return value

    // Find all neighbouring gas particles.
    int num_ngb_gas, ind_ngb_gas[NumPart];
    gas_neighbours(target, &num_ngb_gas, ind_ngb_gas, smoothingFactor);

    // Check that neighbouring gas particles have calculated their local sSFR this timestep.
    int n, n_ind;
    for (n = 0; n < num_ngb_gas; n++) {
        n_ind = ind_ngb_gas[n];
        // If local sSFR not calculated this timestep then do it now.
        if (All.NumCurrentTiStep != P[n_ind].ssfr_local_lastCalcTi) {
            _calc_local_ssfr(n_ind, smoothingFactor);
        }
    }

    //* Compute the kernel-weighted average of the local sSFR values of neighbouring particles. *//
    // Variable definitions.
    double dx, dy, dz, r2, h2, hinv3, hinv4, u, wk, dwk;
    double sum_wk=0, sum_wk_sSFR=0, ssfr_smooth;
    // Softening length 
    double h = smoothingFactor * P[target].Hsml;
    if(h<=0){
        h=All.SofteningTable[0];
    }
    // Calculate required kernel values.
    h2 = h*h;
    hinv4 = 1/(h2*h2);
    hinv3 = hinv4*h;

    // Kernel weights for target particle.
    kernel_main(0, hinv3, hinv4, &wk, &dwk, -1);
    sum_wk += wk;
    sum_wk_sSFR += wk * P[target].ssfr_local;

    // Kernel weights for neighbouring particles.
    for (n=0; n < num_ngb_gas; n++){
        n_ind = ind_ngb_gas[n];

        // Calculate displacement between target particle and its neighbour.
        dx = P[target].Pos[0] - P[n_ind].Pos[0];
        dy = P[target].Pos[1] - P[n_ind].Pos[1];
	    dz = P[target].Pos[2] - P[n_ind].Pos[2];
	    r2 = dx * dx + dy * dy + dz * dz;
	    
        // Calculate kernel weights and weighted local sSFR for neighbouring particles.
	    if( r2 < h2 ) {
            u = sqrt(r2/h2);
            kernel_main(u, hinv3, hinv4, &wk, &dwk, -1);
            sum_wk += wk;
            sum_wk_sSFR += wk * P[n_ind].ssfr_local;
        }  
    }
    
    // Return smoothed value.
    if (sum_wk > 0.){
        ssfr_smooth = sum_wk_sSFR / sum_wk;
    } else {
        ssfr_smooth = 0.;
    }

    return ssfr_smooth;
}

//* Function to calculate local interstellar radiation field from ssfr * //
void cooling_calc_G0_local(int target, double smoothingFactor)
{
    // Smooth the local ssfr by calling appropriate function.
    double ssfr_smooth = _smooth_local_ssfr(target, smoothingFactor);

    // Calculate local ISRF and save in particle struct.
    // If ssfr is zero then set G0 to a default background value of 1.
    if (ssfr_smooth == 0.) {
        P[target].G0_local = 0.;
    } else {
        // The specific star formation rate in the MW (Timothy C. Licquia and Jeffrey A. Newman, 2015), in code units.
        double ssfr_mw = 2.71e-11 * (All.UnitTime_in_s / SEC_PER_YEAR);
        // G0 in the MW, taken from the SIGAME paper (Karen Olsen, Thomas R. Greve, 2017) in Habing units.
        double G0_mw = 0.6;
        // Scaling of the local SFR with MW properties to get G0.
        P[target].G0_local = G0_mw * ssfr_smooth / ssfr_mw;

        //! debug
        if (P[target].G0_local != P[target].G0_local) {
            fprintf(stdout, "G0 = %g, ssfr_mw = %g, ssfr_smooth = %g, G0_mw = %g, target = %d \n", P[target].G0_local, ssfr_mw, ssfr_smooth, G0_mw, target);
        }
    }

    // Save the timestep that this was calculated.
    P[target].G0_local_lastCalcTi = All.NumCurrentTiStep;
}

//* Function to calculate local sSFR of a given gas particle, smoothing to alleviate non-phsyical *//
//* boundaries. *//
int cooling_calculate_local_G0(int target, double smoothingFactor)
{
    // Check particle is gaseous and massive.
    if ( (P[target].Type != 0) || (P[target].Mass <= 0) ) {
        return 0; // Return fail
    }

    // Calculate target particle's G0
    if (All.NumCurrentTiStep != P[target].G0_local_lastCalcTi){
        // If target particle's sSFR has not been updated do it now
        if (All.NumCurrentTiStep != P[target].ssfr_local_lastCalcTi){
            _calc_local_ssfr(target, smoothingFactor);
        }
        _calc_G0_local(target, smoothingFactor);
    }

    // Check particle's values have been calculated and return pass/fail signal.
    if ( (All.NumCurrentTiStep == P[target].ssfr_local_lastCalcTi) && 
         (All.NumCurrentTiStep == P[target].G0_local_lastCalcTi) ) {
        return 1; // Return success
    } else {
        return 0; // Return fail
    }
}
