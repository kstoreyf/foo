/* This file is auto-generated from weight_functions.h.src */
#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif
// # -*- mode: c -*-
#pragma once

#include "defs.h"


#include <stdint.h>

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

// Info about a particle pair that we will pass to the weight function
typedef struct
{
    int64_t nprojbins;
    int64_t nsbins;
    double sqr_s;
    double *supp_sqr; //proper way to have array in struct?

} proj_struct_double;

//typedef double (*weight_func_t_double)(const pair_struct_double*);


//////////////////////////////////
// Projection functions
//////////////////////////////////

static inline void tophat_double(const proj_struct_double *proj, double *u){

    int ins = -1;
    for(int p=0;p<proj->nsbins;p++){
        u[p] = 0;
        if (proj->sqr_s >= proj->supp_sqr[p] && proj->sqr_s < proj->supp_sqr[p+1]){
            ins = p;
        }
    }
    if (ins>=0 && ins<proj->nprojbins){
        u[ins] = 1.0;
    }
}


//////////////////////////////////
// Utility functions
//////////////////////////////////

static inline void compute_amplitudes(int nprojbins, int nd1, int nd2, int nr1, int nr2,
            double *dd, double *dr, double *rd, double *rr, double *qq, double *amps){

    printf("Computing amps\n");
    printf("qq:\n");
    for(int i=0;i<nprojbins;i++){
        for (int j=0; j<nprojbins; j++){
            printf(" %f", qq[i*nprojbins+j]);
        }
        printf("\n");
    }

    // Computer numerator of estimator
    double numerator[nprojbins];
    double qqnorm[nprojbins*nprojbins];
    for (int i=0; i<nprojbins; i++){
        double ddnorm = dd[i]/(nd1*nd2);
        double drnorm = dr[i]/(nd1*nr2);
        double rdnorm = rd[i]/(nr1*nd2);
        double rrnorm = rr[i]/(nr1*nr2);
        numerator[i] = ddnorm - drnorm - rdnorm + rrnorm;
        for (int j=0; j<nprojbins; j++){
            qqnorm[i*nprojbins+j] = qq[i*nprojbins+j]/(nr1*nr2);
        }
    }

	int s;
	// Define all the used matrices
	gsl_matrix *qq_mat = gsl_matrix_alloc(nprojbins, nprojbins);
	gsl_matrix *qq_mat_inv = gsl_matrix_alloc(nprojbins, nprojbins);
	gsl_permutation *perm = gsl_permutation_alloc(nprojbins);
	// Fill the matrix m
	for (int i=0; i<nprojbins; i++){
        for (int j=0; j<nprojbins; j++){
            gsl_matrix_set(qq_mat, i, j, qqnorm[i*nprojbins+j]);
        }
    }
	// Make LU decomposition of matrix m
	gsl_linalg_LU_decomp(qq_mat, perm, &s);
	// Invert the matrix m
	gsl_linalg_LU_invert(qq_mat, perm, qq_mat_inv);

    printf("qqinv:\n");
    for(int i=0;i<nprojbins;i++){
        for (int j=0; j<nprojbins; j++){
            printf(" %f", gsl_matrix_get(qq_mat_inv, i, j));
        }
        printf("\n");
    }

    // Take inner product of qqinv * numerator, get amplitude vector
    // TODO: double check
    for (int i=0; i<nprojbins; i++){
        double aval = 0;
        for (int j=0; j<nprojbins; j++){
            aval += gsl_matrix_get(qq_mat_inv, i, j) * numerator[j];
        }
        amps[i] = aval;
    }
}

static inline void evaluate_xi(int nprojbins, double *amps, int nsvals, double *svals,
                                                int nsbins, double *sbins, double *xi){

    // will need to generalize, projbins won't always be rela   ted to sbins
    double supp_sqr[nsbins];
    //plus 1 because one more edge than number of bins
    for (int i=0; i<nsbins+1; i++){
        supp_sqr[i] = pow(sbins[i], 2);
    }

    // ns: number of s values at which to evaluate xi
    for (int i=0; i<nsvals; i++){
        //get basis function u for given value of s
        double u[nprojbins];
        double sqr_s = pow(svals[i], 2);
        proj_struct_double projdata = {.nprojbins=nprojbins, .nsbins=nsbins, .sqr_s=sqr_s, .supp_sqr=supp_sqr};
        tophat_double(&projdata, u);

        //multiply u by the amplitudes to get xi in that s bin (xi is vector of length ns_toeval)
        double xival = 0;
        for (int j=0; j<nprojbins; j++){
            xival += amps[j]*u[j];

        }
        xi[i] = xival;
    }

}


/* Gives a pointer to the weight function for the given weighting method
 * and instruction set.
 */
//static inline weight_func_t_double get_weight_func_by_method_double(const weight_method_t method){
//    switch(method){
//        case PAIR_PRODUCT:
//            return &pair_product_double;
//        default:
//        case NONE:
//            return NULL;
//    }
//}