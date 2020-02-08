
// This code is used for simulating the extended MDS model described in Liu et al., 2020, Computational Brain & Behavior. 
// The complied version of the C code can be used together with MDS_test.R
// Installation of fftw3 package is required before using this code
// Instructions for installation of fftw3 can be found at http://www.fftw.org/
// Questions can be sent to Qingfang Liu, psychliuqf@gmail.com
// Feb, 2020


#include <errno.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <stddef.h>
#include <time.h>
#include <fftw3.h> // to use FFTW3 package for FFT convolution

// column major index
// column and row indexes start at 0, but nrows is physical
#define idx(r,c,nrows) (r+(c)*nrows)

void MDS_two_choice(double *pL,int *n_region,int *n_time,int *n_trial,int *n_rep,
    double *thresh1,double *thresh2,double *sigma1,double *sigma2,
    double *c0,double *c1,double *c2,double *a1,double *a2,double *a3,
    double *d1,double *d2,double *tau,int *L,int *ds_factor,
    double *A12,double *A34,double *A5,double *A6,
    double *Resp,double *t1,double *t0,double *y_BOLD
    ){
  
	// parameters
	// pL: probability condition of moving leftward - a vector
	// n_region: number of brain regions in the model
	// n_time: iteration time steps
	// n_trial: number of trials in a given experimental stimuli series
	// n_rep: number of repetitions for parameter recovery use
	// thresh1: first threshold for the evidence difference
	// thresh2: 2nd threshold for the response accumulation
	// sigma1: noise term std (control SNR) before responses
  // sigma2: noise term std after responses for visual cortex
	// c0,c1,c2: self-connection strength parameters
	// a1,a2,a3: inter-region connection strength parameters
	// d1,d2: exogenous connectivity matrix coefficients
  // tau: non-decision time
  // L: length of HRF (in msec)
  // ds_factor: downsampling factor from convolved HRF to (hypothetical) BOLD signals
  // Resp: response made 
	// t1: reaction time (before adding the non-decision time)
	// t0: intermediate time
  // y_BOLD: to store all convolved BOLD activity of all ROIs in all trials of all repititions AFTER downsampling
    
    double *SS = malloc(((*n_region) * (*n_time)* (*n_trial)) * sizeof(double)); //trajectories of all ROIs in all trials at current repetition
    double *y = malloc(((*n_region) * (*n_time)* (*n_trial)) * sizeof(double)); //convolved BOLD activity of all ROIs in all trials at current repetition

    int n_elements = (*n_region) * (*n_region); //number of elements in C and D matrices
    double C[n_elements],D[n_elements]; // declare C, D matrix as a vector
    double omega; // declare omega (the noise term)
    double U[(*n_region) * (*n_time)]; // declare the input 
    double S[(*n_region) * (*n_time)]; // declare S matrix for storing the current trial
    double end_state_regions[*n_region]; // decalre a vector to store the end state and to connect with next trial

    int ds_size = ((*n_time)*(*n_trial)*(*n_region))/(*ds_factor); // downsample size for one repetition
	  double ds_samples[ds_size]; // downsampled samples in one repetition

    // some index to use
    int t,k,r,j,i,z;
    double temp_evi_differ,evi_R6;
    int t_temp,t_resp;
    double term1,term2;
    
	  // initialize all Resps to be 9
    for (j=0; j<((*n_trial)*(*n_rep)); j++) {Resp[j]=9;}

    // ################### HRF setup ##############

     double alpha1 = 6,alpha2 = 16,beta1 = 1,beta2 = 1,c = 1.0/6; //fixed parameters in canonical HRF function
     double phi[(*L)*(*n_region)]; //declare phi as a matrix to store HRF basis for each brain region
     double t_hrf; //a temp var to transfer time unit from msec to sec
     double gamma_alpha1 = 120.0;
     double gamma_alpha2_1;gamma_alpha2_1 = 15*14*13*12*11.0;
     double gamma_alpha2_2;gamma_alpha2_2 = 10*9*8*7*6*5*4*3*2.0;

    for (j = 0; j < *L; j++){
          t_hrf = .001 * (j+1);
          phi[idx(0, j, *n_region)] = (*A12) * ((pow(t_hrf,(alpha1-1)) * pow(beta1,alpha1) * exp(-beta1 * t_hrf))/gamma_alpha1 -
                   c * ((pow(t_hrf,(alpha2-1)) * pow(beta2,alpha2) * exp(-beta2 * t_hrf))/(gamma_alpha2_1*gamma_alpha2_2)));
          phi[idx(1, j, *n_region)] = (*A12) * ((pow(t_hrf,(alpha1-1)) * pow(beta1,alpha1) * exp(-beta1 * t_hrf))/gamma_alpha1 -
                   c * ((pow(t_hrf,(alpha2-1)) * pow(beta2,alpha2) * exp(-beta2 * t_hrf))/(gamma_alpha2_1*gamma_alpha2_2)));
          phi[idx(2, j, *n_region)] = (*A34) * ((pow(t_hrf,(alpha1-1)) * pow(beta1,alpha1) * exp(-beta1 * t_hrf))/gamma_alpha1 -
                   c * ((pow(t_hrf,(alpha2-1)) * pow(beta2,alpha2) * exp(-beta2 * t_hrf))/(gamma_alpha2_1*gamma_alpha2_2)));
          phi[idx(3, j, *n_region)] = (*A34) * ((pow(t_hrf,(alpha1-1)) * pow(beta1,alpha1) * exp(-beta1 * t_hrf))/gamma_alpha1 -
                   c * ((pow(t_hrf,(alpha2-1)) * pow(beta2,alpha2) * exp(-beta2 * t_hrf))/(gamma_alpha2_1*gamma_alpha2_2)));
          phi[idx(4, j, *n_region)] = (*A5) * ((pow(t_hrf,(alpha1-1)) * pow(beta1,alpha1) * exp(-beta1 * t_hrf))/gamma_alpha1 -
                   c * ((pow(t_hrf,(alpha2-1)) * pow(beta2,alpha2) * exp(-beta2 * t_hrf))/(gamma_alpha2_1*gamma_alpha2_2)));
          phi[idx(5, j, *n_region)] = (*A6) * ((pow(t_hrf,(alpha1-1)) * pow(beta1,alpha1) * exp(-beta1 * t_hrf))/gamma_alpha1 -
                   c * ((pow(t_hrf,(alpha2-1)) * pow(beta2,alpha2) * exp(-beta2 * t_hrf))/(gamma_alpha2_1*gamma_alpha2_2)));
     }

     // define neural noise in each brain region
	  double nsig[*n_region]; for (j = 0; j < *n_region; j++) {nsig[j] = .05;}

	// #################### end of HRF setup ##############
 
// #############  FFTW convolution setup ###############

    int SIZE = (*n_time)*(*n_trial)+(*L); // size of the linear convolution result
    int SIZE_all; SIZE_all = SIZE * (*n_region); // multiplying n_region on SIZE to store fft result of h for all regions

    double *h,*x,*ifft_result; //real-valued input and output
    x = fftw_malloc ( sizeof ( double ) * SIZE );
    h = fftw_malloc ( sizeof ( double ) * SIZE );
    ifft_result = fftw_malloc ( sizeof ( double ) * SIZE );

    fftw_complex    *fft_result_x, *fft_result_h, *fft_result_prod, *fft_result_h_all; // declare complex variables for fft and ifft
    fftw_plan       plan_forward_x, plan_forward_h,plan_backward; // create fftw plans and plans can be executed many times
 
    fft_result_h  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * SIZE);
    fft_result_x  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * SIZE);
    fft_result_prod  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * SIZE);
    fft_result_h_all = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * SIZE_all);

    plan_forward_x  = fftw_plan_dft_r2c_1d ( SIZE, x, fft_result_x, FFTW_ESTIMATE );
    plan_forward_h  = fftw_plan_dft_r2c_1d ( SIZE, h, fft_result_h, FFTW_ESTIMATE );
    plan_backward = fftw_plan_dft_c2r_1d (SIZE, fft_result_prod, ifft_result, FFTW_ESTIMATE);

	//fft of h for each region
	for (r = 0; r < (*n_region); ++r) {
	// define h
    for(j = 0 ; j < (*L) ; j++) {h[j] = phi[idx(r, j, *n_region)]; }
    for(j = (*L); j < SIZE; j++){h[j] = 0.0;}

	fftw_execute(plan_forward_h);

	for (j = 0; j < SIZE; j++) {
		fft_result_h_all[idx(r,j,*n_region)][0] = fft_result_h[j][0];
		fft_result_h_all[idx(r,j,*n_region)][1] = fft_result_h[j][1];
	}
    }

//####################### end of FFTW convolution setup #################

    // define D matrix elements that never change across simulations
    for (j=0; j<n_elements; j++) {D[j]=0;}
         D[idx(0, 0, *n_region)] = *d1;
         D[idx(1, 1, *n_region)] = *d1;
         D[idx(4, 4, *n_region)] = *d2;

    // define C matrix (elements change from C1 to C2 and back to C1 within a trial)
    for (j=0; j<n_elements; j++) {C[j]=0;} // initialize C matrix at the start of each simulation

    // assume unequal self-connection strength within different brain regions
    C[idx(0, 0, *n_region)] = *c1; // self-connection of R1 and R2 is set to be c1
    C[idx(1, 1, *n_region)] = *c1;
    C[idx(2, 2, *n_region)] = *c2; // self-connection of R3 and R4 is set to be c2
    C[idx(3, 3, *n_region)] = *c2;
    C[idx(4, 4, *n_region)] = *c0; // self-connection of R5 and R6 is set to be c0
    C[idx(5, 5, *n_region)] = *c0;

    C[idx(2, 0, *n_region)] = *a1; // connection strength from R1 to R3
    C[idx(3, 1, *n_region)] = *a1; // connection strength from R2 to R4
    C[idx(5, 4, *n_region)] = *a3; // connection strength from R5 to R6

	// define the external inputs U
	for (j=0; j<((*n_time) * (*n_region)); j++) {U[j]=0;}
	
    for (z = 0; z < (*n_rep); ++z){ //loop through repetitions

    for (j=0; j<((*n_time) * (*n_region) * (*n_trial)); j++) {SS[j]=9;} // initialize SS to be all 9s

    for (i=0; i<(*n_trial); i++) { //loop through total trials
        
        // define state matrix S (changes on each trial)
        for (j=0; j<((*n_time) * (*n_region)); j++) {S[j]=0;} // initialize S to be all zeros
        
        // to make previous trials affect the next one
        if (i>0){
        	for (j = 0; j < (*n_region); ++j){
	 			S[idx(j,0,*n_region)] = end_state_regions[j];}
        }
		
        // define the external inputs U
		
        // input to R1 and R2, allow randomness at trial level
        for (j=0; j<(*n_time); j++) {
            U[idx(0, j, *n_region)] = rbinom(1000,pL[i]);
            U[idx(1, j, *n_region)] = 1000-U[idx(0, j, *n_region)];
            U[idx(0, j, *n_region)] = U[idx(0, j, *n_region)]/100.0;
            U[idx(1, j, *n_region)] = U[idx(1, j, *n_region)]/100.0;
        }
        
        // input to R5 (fixed input but is reset to zero as well after response is made)
    	for (j=0; j<(*n_time); j++) {U[idx(4, j, *n_region)]=1.0;}

        temp_evi_differ = 0; evi_R6 = 0; // some index to track the accumulated evidences
        t_temp = 0; t_resp = 0;

        // loop through time steps
        for (t=1; t<(*n_time); t++) {
            // loop through regions
            for (r=0; r<(*n_region); r++) {
                term1=0; term2=0; 
                for (k=0; k<(*n_region); k++) {term1 += C[idx(r, k, *n_region)] * S[idx(k, t-1, *n_region)];}
                term2 = D[idx(r,r,*n_region)] * U[idx(r,t,*n_region)];
                omega = rnorm(0.0,*sigma1); 
                S[idx(r, t, *n_region)] = term1 + term2 + omega; 
            }
            
            // note: term1 is the adding to region r from endogenous connectivity, and term2 is the adding to region r from external inputs
			      // When defining the noise term, assume no correlation b/t regions, and noises in all regions share the same distribution
            // define the noise term in every region so their noise ossilations can be different

            temp_evi_differ = S[idx(2, t, *n_region)]-S[idx(3, t, *n_region)];
            evi_R6 += S[idx(5, t, *n_region)];

            //check if the absolute accumulated evidence difference reaches the threshold
            if (temp_evi_differ > (*thresh1) || temp_evi_differ < (-1)* (*thresh1)) { 
                C[idx(4, 2, *n_region)] = *a2;
                C[idx(4, 3, *n_region)] = *a2;
                t_temp = t; // note down the current time
                if (temp_evi_differ>0) { // note down which is the winner accumulator
                    Resp[idx(z,i,*n_rep)] = 1.0;
                } else {
                    Resp[idx(z,i,*n_rep)] = 2.0;
                }
                break; //break the loop
            }
        }
        
        // start the second loop after changing matrix C
        for (t = t_temp+1; t<(*n_time); t++) {
            if (evi_R6 > *thresh2) {
              t_resp = t;
            	break;}
            
            for (r=0; r<(*n_region); r++) {
                term1=0;term2=0; 
                for (k=0; k<(*n_region); k++) {term1 += C[idx(r, k, *n_region)] * S[idx(k, t-1, *n_region)];}
                term2 = D[idx(r,r,*n_region)] * U[idx(r,t,*n_region)];
                omega = rnorm(0.0,*sigma1); 
                S[idx(r, t, *n_region)] = term1 + term2 + omega; 
            }
           
            evi_R6 += S[idx(5, t, *n_region)]; // calculate the evidence in R6 up to time t
        }

		   t1[idx(z,i,*n_rep)] = t_resp; //total RT in this trial
		   t0[idx(z,i,*n_rep)] = t_temp; // intermediate time in this trial

		////////////// let the evoluation not stop after making a response ////////////////////

        // reset all inputs to be zero
        for (j = t_resp + 1; j<(*n_time); j++) {
            U[idx(0, j, *n_region)] = 0.0;
            U[idx(1, j, *n_region)] = 0.0;
            U[idx(4, j, *n_region)] = 0.0;   
        }

        // change back engodenous connectivity matrix from C2 to C1

        C[idx(4, 2, *n_region)] = 0.0;
        C[idx(4, 3, *n_region)] = 0.0;

        // after the response is made, the noise variances of R1 and R2 (only) decrease so that fluctuation goes back to around zero
        
        for (t = t_resp+1; t<(*n_time); t++) {
             for (r=0; r<2; r++) {
               term1 = 0; term2 = 0; 
                for (k=0; k<(*n_region); k++) {term1 += C[idx(r, k, *n_region)] * S[idx(k, t-1, *n_region)];}
               		term2 = D[idx(r,r,*n_region)] * U[idx(r,t,*n_region)];
               		omega = rnorm(0.0,*sigma2); 
               		S[idx(r, t, *n_region)] = term1 + term2 + omega; 
            }

            for (r=2; r<(*n_region); r++) {
               term1 = 0; term2 = 0; 
                for (k=0; k<(*n_region); k++) {term1 += C[idx(r, k, *n_region)] * S[idx(k, t-1, *n_region)];}
               		term2 = D[idx(r,r,*n_region)] * U[idx(r,t,*n_region)];
               		omega = rnorm(0.0,*sigma1); 
               		S[idx(r, t, *n_region)] = term1 + term2 + omega; 
            }
        }

		// keep records of trajectories of all ROIs in all trials as SS
        for (j = 0; j < ((*n_time) * (*n_region)); j++){
        	SS[j + i * ((*n_region) * (*n_time)) ] = S[j];
        }

    // keep records of the end states across regions at current trial
        for (j = 0; j < (*n_region); j++){
        	end_state_regions[j] = S[idx(j,*n_time-1,*n_region)]; 
        }
        
    } // end of the current trial

// //////////// start of FFTW convolution for the current repetition (not on the trial level) ////////////////////////////

	for (r=0; r<(*n_region); r++){ // loop through regions

	for( i = 0 ; i < (*n_time)*(*n_trial) ; i++ ) {
	    x[i] = SS[idx(r, i, *n_region)]; 
	    }
	for (i = (*n_time)*(*n_trial); i < SIZE; ++i){
	    x[i] = 0.0; 
	}
	    
	fftw_execute( plan_forward_x );

	// multiply in the frequency domain

	for (i = 0; i < SIZE; ++i){
	  fft_result_prod[i][0] = fft_result_x[i][0] * fft_result_h_all[idx(r, i, *n_region)][0] - fft_result_x[i][1] * fft_result_h_all[idx(r, i, *n_region)][1];
	  fft_result_prod[i][1] = fft_result_x[i][0] * fft_result_h_all[idx(r, i, *n_region)][1] + fft_result_x[i][1] * fft_result_h_all[idx(r, i, *n_region)][0];
	}

	fftw_execute( plan_backward );

	// normalize and discard the samples after n_time*n_trial in the total length of SIZE
	for (i = 0; i < (*n_time) * (*n_trial); i++) {
	    y[i + r * (*n_time) * (*n_trial)] = ifft_result[i]/SIZE;
	  }
	} // end of loop of region

	//################### downsample the convolved neural signal #################################

	for (j = 0; j < ds_size; ++j){ds_samples[j] = y[j*(*ds_factor)];}

	// add noise
	for (r = 0; r < *n_region; ++r){
		for (j = 0; j< (*n_time) * (*n_trial) / (*ds_factor); ++j){
			ds_samples[idx(j,r,(*n_time) * (*n_trial) / (*ds_factor))] += rnorm(0.0,nsig[r]); 
	}
	}
	//###############  keep records of Neural activation in all repetitions as BOLD and downsample

	for (j = 0; j < ds_size; j++){y_BOLD[j + z * ds_size ] = ds_samples[j];}

	  } // end of the current repetition

	fftw_destroy_plan( plan_forward_x );
	fftw_destroy_plan( plan_forward_h );
	fftw_destroy_plan( plan_backward );

// free all variables to avoid "memory leak"  (very important!)
	fftw_free( x );
	fftw_free( h );
	fftw_free( fft_result_x );
	fftw_free( fft_result_h );
	fftw_free( fft_result_prod);
	fftw_free( ifft_result );
	fftw_free( fft_result_h_all);

  free(SS);
  free(y);
  // return RNG
  PutRNGstate();
}

