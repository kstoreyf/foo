/* This file is auto-generated from countpairs_s_mu_mocks_impl.h.src */
#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif
// # -*- mode: c -*-
/* File: countpairs_s_mu_mocks_impl.h.src */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "defs.h" //for struct config_options
#include "weight_defs_double.h"
#include <inttypes.h> //for uint64_t

#include "countpairs_s_mu_mocks.h" //for definition of results_countpairs_mocks

    extern void interrupt_handler_countpairs_s_mu_mocks_double(int signo);

    typedef int (*countpairs_mocks_func_ptr_double)(const int64_t N0, double *x0, double *y0, double *z0, double *d0, const weight_struct_double *weights0,
                                                    const int64_t N1, double *x1, double *y1, double *z1, double *d1, const weight_struct_double *weights1,
                                                    const int same_cell,
                                                    const int fast_divide,
                                                    const double smax, const double smin, const int nsbin,
                                                    const int nmu_bins, const double *supp_sqr,
                                                    const double mu_max,
                                                    double *src_savg, uint64_t *src_npairs, double *src_projpairs,
                                                    double *src_projpairs_tensor,
                                                    double *src_weightavg, const weight_method_t weight_method);

    extern countpairs_mocks_func_ptr_double countpairs_s_mu_mocks_driver_double(const struct config_options *options) __attribute__((warn_unused_result));

    extern int countpairs_mocks_s_mu_double(const int64_t ND1, double *theta1, double *phi1, double *czD1,
                                            const int64_t ND2, double *theta2, double *phi2, double *czD2,
                                            const int numthreads,
                                            const int autocorr,
                                            const char *sbinfile,
                                            const double mu_max,
                                            const int nmu_bins,
                                            const int cosmology,
                                            results_countpairs_mocks_s_mu *results,
                                            struct config_options *options, struct extra_options *extra);

    extern int check_ra_dec_cz_s_mu_double(const int64_t N, double *phi, double *theta, double *cz);

#ifdef __cplusplus
}
#endif
