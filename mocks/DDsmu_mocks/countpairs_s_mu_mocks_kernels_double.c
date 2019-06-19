/* This file is auto-generated from countpairs_s_mu_mocks_kernels.c.src */
#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif
// # -*- mode: c -*-
/* File: countpairs_s_mu_mocks_kernels.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "defs.h"
#include "function_precision.h"
#include "utils.h"

#include "weight_functions_double.h"
#include "proj_functions_double.h"

#if defined(__AVX__)
#include "avx_calls.h"

static inline int countpairs_s_mu_mocks_avx_intrinsics_double(const int64_t N0, double *x0, double *y0, double *z0, double *d0, const weight_struct_double *weights0,
                                                              const int64_t N1, double *x1, double *y1, double *z1, double *d1, const weight_struct_double *weights1,
                                                              const int same_cell,
                                                              const int fast_divide,
                                                              const double smax, const double smin, const int nsbin,const int nmu_bins,
                                                              const double *supp_sqr, const double mu_max,
                                                              double *src_savg, uint64_t *src_npairs, double *src_projpairs,
                                                              double *src_projpairs_tensor, double *src_weightavg, const weight_method_t weight_method)
{
    if(N0 == 0 || N1 == 0) {
        return EXIT_SUCCESS;
    }

    if(src_npairs == NULL) {
        return EXIT_FAILURE;
    }

    const int32_t need_savg = src_savg != NULL;
    const int32_t need_weightavg = src_weightavg != NULL;

    const int64_t totnbins = (nmu_bins+1)*(nsbin+1);
    //TODO: make nprojbins parameter
    const int64_t nprojbins = nsbin-1;
    const double sqr_mumax = mu_max*mu_max;
    const double sqr_smax  = smax*smax;
    const double sqr_smin  = smin*smin;

    AVX_FLOATS m_supp_sqr[nsbin];
    AVX_FLOATS m_kbin[nsbin];
    for(int i=0;i<nsbin;i++) {
        m_supp_sqr[i] = AVX_SET_FLOAT(supp_sqr[i]);
        m_kbin[i] = AVX_SET_FLOAT((double) i);
    }

    uint64_t npairs[totnbins];
    const double dmu = mu_max/(double) nmu_bins;
    const double inv_dmu = 1.0/dmu;
    double savg[totnbins], weightavg[totnbins], projpairs[nprojbins];
    //double projpairs_tensor[nprojbins][nprojbins];
    double projpairs_tensor[nprojbins*nprojbins];
    for(int i=0;i<totnbins;i++) {
        npairs[i] = ZERO;
        if(need_savg) {
            savg[i] = ZERO;
        }
        if(need_weightavg){
            weightavg[i] = ZERO;
        }
    }
    for(int i=0;i<nprojbins;i++) {
        projpairs[i] = ZERO;
        for(int j=0;j<nprojbins;j++) {
            projpairs_tensor[i*nprojbins+j] = ZERO;
        }
    }



    // A copy whose pointers we can advance
    weight_struct_double local_w0 = {.weights={NULL}, .num_weights=0},
                         local_w1 = {.weights={NULL}, .num_weights=0};
    pair_struct_double pair = {.num_weights=0};
    avx_weight_func_t_double avx_weight_func = NULL;
    weight_func_t_double fallback_weight_func = NULL;
    if(need_weightavg){
        // Same particle list, new copy of num_weights pointers into that list
        local_w0 = *weights0;
        local_w1 = *weights1;

        pair.num_weights = local_w0.num_weights;

        avx_weight_func = get_avx_weight_func_by_method_double(weight_method);
        fallback_weight_func = get_weight_func_by_method_double(weight_method);
    }

    int64_t prev_j = 0, n_off = 0;
    for(int64_t i=0;i<N0;i++) {
        const double xpos = *x0++;
        const double ypos = *y0++;
        const double zpos = *z0++;
        const double dpos = *d0++;
        for(int w = 0; w < pair.num_weights; w++){
            // local_w0.weights[w] is a pointer to a float in the particle list of weights,
            // just as x0 is a pointer into the list of x-positions.
            // The advancement of the local_w0.weights[w] pointer should always mirror x0.
            pair.weights0[w].a = AVX_SET_FLOAT(*(local_w0.weights[w])++);
        }

        int64_t j;
        if(same_cell == 1) {
            d1++; n_off++;
            j = i+1;
        } else {
            for(;prev_j<N1;prev_j++) {
                const double dz = *d1 - dpos;
                if(dz > -smax) break;
                d1++; n_off++;
            }
            if(prev_j == N1) {
                break;
            }
            j = prev_j;
        }
        double *locald1 = d1;
        double *localx1 = x1 + n_off;
        double *localy1 = y1 + n_off;
        double *localz1 = z1 + n_off;
        for(int w = 0; w < local_w1.num_weights; w++){
            local_w1.weights[w] = weights1->weights[w] + n_off;
        }

        AVX_FLOATS m_xpos = AVX_SET_FLOAT(xpos);
        AVX_FLOATS m_ypos = AVX_SET_FLOAT(ypos);
        AVX_FLOATS m_zpos = AVX_SET_FLOAT(zpos);
        AVX_FLOATS m_dpos = AVX_SET_FLOAT(dpos);
        union int8 {
            AVX_INTS m_ibin;
            int ibin[AVX_NVEC];
        };

        union float8{
            AVX_FLOATS m_sep;
            double sep[AVX_NVEC];
        };

        const AVX_FLOATS m_sqr_smax = AVX_SET_FLOAT(sqr_smax);
        const AVX_FLOATS m_sqr_smin = AVX_SET_FLOAT(sqr_smin);
        const AVX_FLOATS m_sqr_mumax = AVX_SET_FLOAT(sqr_mumax);
        const AVX_FLOATS m_inv_dmu = AVX_SET_FLOAT(inv_dmu);
        const AVX_FLOATS m_nmu_bins = AVX_SET_FLOAT((double) nmu_bins);
        const AVX_FLOATS m_zero = AVX_SET_FLOAT(ZERO);
        const AVX_FLOATS m_one = AVX_SET_FLOAT((double) 1);

        for(;j<=(N1-AVX_NVEC);j+=AVX_NVEC){
            const AVX_FLOATS m_x2 = AVX_LOAD_FLOATS_UNALIGNED(localx1);
            const AVX_FLOATS m_y2 = AVX_LOAD_FLOATS_UNALIGNED(localy1);
            const AVX_FLOATS m_z2 = AVX_LOAD_FLOATS_UNALIGNED(localz1);
            const AVX_FLOATS m_d2 = AVX_LOAD_FLOATS_UNALIGNED(locald1);

            localx1 += AVX_NVEC;
            localy1 += AVX_NVEC;
            localz1 += AVX_NVEC;
            locald1 += AVX_NVEC;

            for(int w = 0; w < pair.num_weights; w++){
                pair.weights1[w].a = AVX_LOAD_FLOATS_UNALIGNED(local_w1.weights[w]);
                local_w1.weights[w] += AVX_NVEC;
            }

            union float8_weights{
                AVX_FLOATS m_weights;
                double weights[NVEC];
            };
            union float8_weights union_mweight;

            const AVX_FLOATS m_perpx = AVX_SUBTRACT_FLOATS(m_xpos, m_x2);
            const AVX_FLOATS m_perpy = AVX_SUBTRACT_FLOATS(m_ypos, m_y2);
            const AVX_FLOATS m_perpz = AVX_SUBTRACT_FLOATS(m_zpos, m_z2);

            const AVX_FLOATS m_parx = AVX_ADD_FLOATS(m_x2, m_xpos);
            const AVX_FLOATS m_pary = AVX_ADD_FLOATS(m_y2, m_ypos);
            const AVX_FLOATS m_parz = AVX_ADD_FLOATS(m_z2, m_zpos);

            AVX_FLOATS m_sqr_mu, m_sqr_s;
            {
                /*
                  //Technically l := 1/2 (v1 + v2) but the factor of 1/2 occurs both in numerator and denominator
                  and cancels out.

                  s \dot l := (parx*perpx + pary*perpy + parz*perp)
                           := (x1 + x2)*(x1 - x2) + (y1 + y2)*(y1 - y2) + (z1 + z2)*(z1 - z2)
                           := (x1^2 + y1^2 + z1^2) - (x2^2 + y2^2 + z2^2)
                           := d1^2 - d2^2
                */
                const AVX_FLOATS m_s_dot_l =  AVX_SUBTRACT_FLOATS(AVX_SQUARE_FLOAT(m_d2), AVX_SQUARE_FLOAT(m_dpos));
                const AVX_FLOATS m_sqr_s_dot_l = AVX_SQUARE_FLOAT(m_s_dot_l);// numerator := |s.l|^2
                const AVX_FLOATS m_sqr_perpx = AVX_SQUARE_FLOAT(m_perpx);
                const AVX_FLOATS m_sqr_perpy = AVX_SQUARE_FLOAT(m_perpy);
                const AVX_FLOATS m_sqr_perpz = AVX_SQUARE_FLOAT(m_perpz);
                m_sqr_s = AVX_ADD_FLOATS(m_sqr_perpx, AVX_ADD_FLOATS(m_sqr_perpy, m_sqr_perpz));//3-d separation

                //Create a mask where s^2 < smax^2
                const AVX_FLOATS m_mask_3d_sep = AVX_COMPARE_FLOATS(m_sqr_s, m_sqr_smax, _CMP_LT_OQ);
                if(AVX_TEST_COMPARISON(m_mask_3d_sep) == 0) {
                    continue;
                }
                const AVX_FLOATS m_sqr_norm_l = AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_parx),
                                                               AVX_ADD_FLOATS(AVX_SQUARE_FLOAT(m_pary),
                                                                              AVX_SQUARE_FLOAT(m_parz)));

                // \mu^2 := cos^2(\theta_between_s_and_l) = |s.l|^2 / (|s|^2 * |l|^2)
                const AVX_FLOATS m_sqr_norm_l_norm_s = AVX_MULTIPLY_FLOATS(m_sqr_norm_l, m_sqr_s);
                if (fast_divide == 0) {
                    m_sqr_mu = AVX_DIVIDE_FLOATS(m_sqr_s_dot_l, m_sqr_norm_l_norm_s);
                    //The divide is the actual operation we need
                    // but divides are about 10x slower than multiplies. So, I am replacing it
                    //with a approximate reciprocal in floating point
                    // + 2 iterations of newton-raphson in case of double
                } else {
                    //following blocks do an approximate reciprocal followed by two iterations of Newton-Raphson

#ifndef DOUBLE_PREC
                    //Taken from Intel's site: https://software.intel.com/en-us/articles/wiener-filtering-using-intel-advanced-vector-extensions
                    // (which has bugs in it, just FYI). Plus, https://techblog.lankes.org/2014/06/16/avx-isnt-always-faster-then-see/
                    __m256 rc  = _mm256_rcp_ps(m_sqr_norm_l_norm_s);
#else
                    //we have to do this for doubles now.
                    //if the vrcpps instruction is not generated, there will
                    //be a ~70 cycle performance hit from switching between
                    //AVX and SSE modes.
                    __m128 float_tmp1 =  _mm256_cvtpd_ps(m_sqr_norm_l_norm_s);
                    __m128 float_inv_tmp1 = _mm_rcp_ps(float_tmp1);
                    AVX_FLOATS rc = _mm256_cvtps_pd(float_inv_tmp1);
#endif//DOUBLE_PREC

                  //We have the double->float->approx. reciprocal->double process done.
                  //Now improve the accuracy of the divide with newton-raphson.

                  //Ist iteration of NewtonRaphson
                  AVX_FLOATS two = AVX_SET_FLOAT((double) 2.0);
                  AVX_FLOATS rc1 = AVX_MULTIPLY_FLOATS(rc,
                                                       AVX_SUBTRACT_FLOATS(two,
                                                                           AVX_MULTIPLY_FLOATS(m_sqr_norm_l_norm_s,rc)));
                  //2nd iteration of NewtonRaphson
                  AVX_FLOATS rc2 = AVX_MULTIPLY_FLOATS(rc1,
                                                       AVX_SUBTRACT_FLOATS(two,
                                                                           AVX_MULTIPLY_FLOATS(m_sqr_norm_l_norm_s,rc1)));
                  m_sqr_mu = AVX_MULTIPLY_FLOATS(m_sqr_s_dot_l,rc2);
                } //end of FAST_DIVIDE
            }

            const AVX_FLOATS m_mu = AVX_SQRT_FLOAT(m_sqr_mu);

            AVX_FLOATS m_mask_left;
            //Do the mask filters in a separate scope
            {
                const AVX_FLOATS m_mask_mumax = AVX_COMPARE_FLOATS(m_sqr_mu,m_sqr_mumax,_CMP_LT_OQ);
                const AVX_FLOATS m_smax_mask = AVX_COMPARE_FLOATS(m_sqr_s, m_sqr_smax, _CMP_LT_OQ);
                const AVX_FLOATS m_smin_mask = AVX_COMPARE_FLOATS(m_sqr_s, m_sqr_smin, _CMP_GE_OQ);
                const AVX_FLOATS m_s_mask = AVX_BITWISE_AND(m_smax_mask, m_smin_mask);

                m_mask_left = AVX_BITWISE_AND(m_mask_mumax, m_s_mask);
                if(AVX_TEST_COMPARISON(m_mask_left)==0) {
                    continue;
                }
                m_sqr_s = AVX_BLEND_FLOATS_WITH_MASK(m_zero,m_sqr_s,m_mask_left);
                m_sqr_mu  = AVX_BLEND_FLOATS_WITH_MASK(m_sqr_mumax,m_sqr_mu,m_mask_left);
            }

            union float8 union_msep;
            if(need_savg) {
                union_msep.m_sep = AVX_SQRT_FLOAT(m_sqr_s);
            }
            if(need_weightavg){
                pair.dx.a = m_perpx;
                pair.dy.a = m_perpy;
                pair.dz.a = m_perpz;

                pair.parx.a = m_parx;
                pair.pary.a = m_pary;
                pair.parz.a = m_parz;

                union_mweight.m_weights = avx_weight_func(&pair);
            }

            const AVX_FLOATS m_mask = m_mask_left;
            AVX_FLOATS m_sbin = AVX_SET_FLOAT((double) 0);
            for(int kbin=nsbin-1;kbin>=1;kbin--) {
                const AVX_FLOATS m_mask_low = AVX_COMPARE_FLOATS(m_sqr_s,m_supp_sqr[kbin-1],_CMP_GE_OQ);
                const AVX_FLOATS m_bin_mask = AVX_BITWISE_AND(m_mask_low,m_mask_left);
                m_sbin = AVX_BLEND_FLOATS_WITH_MASK(m_sbin,m_kbin[kbin], m_bin_mask);
                m_mask_left = AVX_COMPARE_FLOATS(m_sqr_s, m_supp_sqr[kbin-1],_CMP_LT_OQ);
                if(AVX_TEST_COMPARISON(m_mask_left) == 0) {
                    break;
                }
            }

            /* Compute the 1-D index to the [sbin, mubin] := sbin*(nmu_bins+1) + mubin */
            const AVX_FLOATS m_tmp2 = AVX_MULTIPLY_FLOATS(m_mu,m_inv_dmu);
            const AVX_FLOATS m_mubin = AVX_BLEND_FLOATS_WITH_MASK(m_nmu_bins, m_tmp2, m_mask);
            const AVX_FLOATS m_nmu_bins_p1 = AVX_ADD_FLOATS(m_nmu_bins,m_one);
            const AVX_FLOATS m_binproduct = AVX_ADD_FLOATS(AVX_MULTIPLY_FLOATS(m_sbin,m_nmu_bins_p1),m_mubin);
            union int8 union_finalbin;
            union_finalbin.m_ibin = AVX_TRUNCATE_FLOAT_TO_INT(m_binproduct);

            // proj index



#if  __INTEL_COMPILER
#pragma unroll(AVX_NVEC)
#endif
            for(int jj=0;jj<AVX_NVEC;jj++) {
                const int ibin=union_finalbin.ibin[jj];
                //TODO: handle weight better if need for projpairs
                const double weight = union_mweight.weights[jj];
                npairs[ibin]++;
                if(need_savg) {
                    savg[ibin] += union_msep.sep[jj];
                }
                if(need_weightavg){
                    //const double weight = union_mweight.weights[jj];
                    weightavg[ibin] += weight;
                }
                //TODO: separate when doing general proj
                //projpairs[ibin] += weight;
            }

//            //TODO: where should this go?? smarter way within avx shit?
//            int ins = -1;
//            for(int p=0;p<nprojbins-1;p++){
//                if (m_sqr_s>=m_supp_sqr[p] && m_sqr_s<m_supp_sqr[p+1]){
//                    ins = p;
//                    break;
//                }
//            }
//            if (ins>=0 && ins<nprojbins){
//                projpairs[ins] += 1.0;
//            }
//            //outer product to get tensor value
//            for(int p1=0;p1<nprojbins-1;p1++){
//                for(int p2=0;p2<nprojbins-1;p2++){
//                    projpairs_tensor[p1*nprojbins+p2] += projpairs[p1]*projpairs[p2];
//                }
//            }

        }//AVX j loop



        //Take care of the remainder
        for(;j<N1;j++) {
            const double parx = xpos + *localx1;
            const double pary = ypos + *localy1;
            const double parz = zpos + *localz1;

            const double perpx = xpos - *localx1;
            const double perpy = ypos - *localy1;
            const double perpz = zpos - *localz1;
            /*
              s := (perpx, perpy, perpz)
              l := 1/2 (parx, pary, parz)  //ignoring the factor 1/2 since it cancels out in both numerator and denominator
              
              s \dot l := (parx*perpx + pary*perpy + parz*perpz)
                       := (x1 + x2)*(x1 - x2) + (y1 + y2)*(y1 - y2) + (z1 + z2)*(z1 - z2)
                       := (x1^2 + y1^2 + z1^2) - (x2^2 + y2^2 + z2^2)
                       := d1^2 - d2^2
            */
            const double s_dot_l = dpos*dpos - (*locald1) * (*locald1);// s \dot l

            localx1++;localy1++;localz1++;locald1++;

            for(int w = 0; w < pair.num_weights; w++){
                pair.weights1[w].d = *local_w1.weights[w]++;
            }


            const double sqr_s = perpx*perpx + perpy*perpy + perpz*perpz;
            if(sqr_s >= sqr_smax || sqr_s < sqr_smin) continue;

            const double norm_l = (parx*parx + pary*pary + parz*parz);// := |l|^2
            const double sqr_s_dot_l = s_dot_l * s_dot_l;
            const double sqr_mu = sqr_s_dot_l/(norm_l * sqr_s);
            const int mubin  = (sqr_mu >= sqr_mumax) ? nmu_bins:(int) (SQRT(sqr_mu)*inv_dmu);
            double s, pairweight;
            if(need_savg) {
                s = SQRT(sqr_s);
            }
            if(need_weightavg){
                pair.dx.d = perpx;
                pair.dy.d = perpy;
                pair.dz.d = perpz;

                pair.parx.d = parx;
                pair.pary.d = pary;
                pair.parz.d = parz;

                pairweight = fallback_weight_func(&pair);
            }

            for(int kbin=nsbin-1;kbin>=1;kbin--) {
                if(sqr_s >= supp_sqr[kbin-1]) {
                    const int ibin = kbin*(nmu_bins+1) + mubin;
                    npairs[ibin]++;
                    if(need_savg) {
                        savg[ibin] += s;
                    }
                    if(need_weightavg){
                        weightavg[ibin] += pairweight;
                    }
                    //TODO: make general
                    //projpairs[ibin] += pairweight;
                    break;
                }
            }
            //TODO: what's the deal with remainder??
//            int ins = -1;
//            for(int p=0;p<nprojbins-1;p++){
//                if (sqr_s>=supp_sqr[p] && sqr_s<supp_sqr[p+1]){
//                    ins = p;
//                    break;
//                }
//            }
//            if (ins>=0 && ins<nprojbins){
//                projpairs[ins] += 1.0;
//            }
//            //outer product to get tensor value
//            for(int p1=0;p1<nprojbins-1;p1++){
//                for(int p2=0;p2<nprojbins-1;p2++){
//                    projpairs_tensor[p1*nprojbins+p2] += projpairs[p1]*projpairs[p2];
//                }
//            }
        }//remainder jloop
    }//i-loop

    for(int i=0;i<totnbins;i++) {
        src_npairs[i] += npairs[i];
        if(need_savg) {
            src_savg[i] += savg[i];
        }
        if(need_weightavg) {
            src_weightavg[i] += weightavg[i];
        }
    }
    for(int i=0;i<nprojbins;i++) {
        src_projpairs[i] += projpairs[i];
        for(int j=0;j<nprojbins;j++) {
            src_projpairs_tensor[i*nprojbins+j] += projpairs_tensor[i*nprojbins+j];
        }
    }
    return EXIT_SUCCESS;
}
#endif //AVX


#if defined(__SSE4_2__)
#include "sse_calls.h"

static inline int countpairs_s_mu_mocks_sse_intrinsics_double(const int64_t N0, double *x0, double *y0, double *z0, double *d0, const weight_struct_double *weights0,
                                                              const int64_t N1, double *x1, double *y1, double *z1, double *d1, const weight_struct_double *weights1,
                                                              const int same_cell,
                                                              const int fast_divide,
                                                              const double smax, const double smin, const int nsbin,
                                                              const int nmu_bins, const double *supp_sqr, const double mu_max,
                                                              double *src_savg, uint64_t *src_npairs, double *src_projpairs,
                                                              double *src_projpairs_tensor,
                                                              double *src_weightavg, const weight_method_t weight_method)
{
    if(N0 == 0 || N1 == 0) {
        return EXIT_SUCCESS;
    }
    if(src_npairs == NULL) {
        return EXIT_FAILURE;
    }

    const int32_t need_savg = src_savg != NULL;
    const int32_t need_weightavg = src_weightavg != NULL;
    (void) fast_divide; //unused

    const int64_t totnbins = (nmu_bins+1)*(nsbin+1);
    const int64_t nprojbins = nsbin-1;
    const double sqr_mumax = mu_max*mu_max;
    const double sqr_smax  = smax*smax;
    const double sqr_smin  = smin*smin;

    SSE_FLOATS m_supp_sqr[nsbin];
    SSE_FLOATS m_kbin[nsbin];
    for(int i=0;i<nsbin;i++) {
        m_supp_sqr[i] = SSE_SET_FLOAT(supp_sqr[i]);
        m_kbin[i] = SSE_SET_FLOAT((double) i);
    }

    uint64_t npairs[totnbins];
    const double dmu = mu_max/(double) nmu_bins;
    const double inv_dmu = 1.0/dmu;
    double savg[totnbins], weightavg[totnbins], projpairs[nprojbins];
    double projpairs_tensor[nprojbins*nprojbins];
    for(int64_t i=0;i<totnbins;i++) {
        npairs[i] = ZERO;
        if (need_savg) {
            savg[i] = ZERO;
        }
        if(need_weightavg){
            weightavg[i] = ZERO;
        }
    }
    for(int64_t i=0;i<nprojbins;i++) {
        projpairs[i] = ZERO;
        for(int j=0;j<nprojbins;j++) {
            projpairs_tensor[i*nprojbins+j] = ZERO;
        }
    }


    // A copy whose pointers we can advance
    weight_struct_double local_w0 = {.weights={NULL}, .num_weights=0},
                         local_w1 = {.weights={NULL}, .num_weights=0};
    pair_struct_double pair = {.num_weights=0};
    sse_weight_func_t_double sse_weight_func = NULL;
    weight_func_t_double fallback_weight_func = NULL;
    if(need_weightavg){
      // Same particle list, new copy of num_weights pointers into that list
      local_w0 = *weights0;
      local_w1 = *weights1;

      pair.num_weights = local_w0.num_weights;

      sse_weight_func = get_sse_weight_func_by_method_double(weight_method);
      fallback_weight_func = get_weight_func_by_method_double(weight_method);
    }

    int64_t prev_j=0, n_off = 0;
    for(int64_t i=0;i<N0;i++) {
        const double xpos = *x0++;
        const double ypos = *y0++;
        const double zpos = *z0++;
        const double dpos = *d0++;
        for(int w = 0; w < pair.num_weights; w++){
            // local_w0.weights[w] is a pointer to a float in the particle list of weights,
            // just as x0 is a pointer into the list of x-positions.
            // The advancement of the local_w0.weights[w] pointer should always mirror x0.
            pair.weights0[w].s = SSE_SET_FLOAT(*local_w0.weights[w]++);
        }

        int64_t j;
        if(same_cell == 1) {
            d1++; n_off++;
            j = i+1;
        } else {
            for(;prev_j<N1;prev_j++) {
                const double dz = *d1 - dpos;
                if(dz > -smax) break;
                d1++; n_off++;
            }
            if(prev_j == N1) {
                break;
            }
            j = prev_j;
        }
        double *locald1 = d1;
        double *localx1 = x1 + n_off;
        double *localy1 = y1 + n_off;
        double *localz1 = z1 + n_off;
        for(int w = 0; w < local_w1.num_weights; w++){
            local_w1.weights[w] = weights1->weights[w] + n_off;
        }

        const SSE_FLOATS m_xpos = SSE_SET_FLOAT(xpos);
        const SSE_FLOATS m_ypos = SSE_SET_FLOAT(ypos);
        const SSE_FLOATS m_zpos = SSE_SET_FLOAT(zpos);
        const SSE_FLOATS m_dpos = SSE_SET_FLOAT(dpos);

        union int8 {
            SSE_INTS m_ibin;
            int ibin[SSE_NVEC];
        };


        union float8{
            SSE_FLOATS m_sep;
            double sep[SSE_NVEC];
        };

        const SSE_FLOATS m_sqr_smax = SSE_SET_FLOAT(sqr_smax);
        const SSE_FLOATS m_sqr_smin = SSE_SET_FLOAT(sqr_smin);
        const SSE_FLOATS m_sqr_mumax = SSE_SET_FLOAT(sqr_mumax);
        const SSE_FLOATS m_inv_dmu = SSE_SET_FLOAT(inv_dmu);
        const SSE_FLOATS m_nmu_bins = SSE_SET_FLOAT((double) nmu_bins);
        const SSE_FLOATS m_zero = SSE_SET_FLOAT(ZERO);
        const SSE_FLOATS m_one = SSE_SET_FLOAT((double) 1);

        for(;j<=(N1-SSE_NVEC);j+=SSE_NVEC){
            const SSE_FLOATS m_x2 = SSE_LOAD_FLOATS_UNALIGNED(localx1);
            const SSE_FLOATS m_y2 = SSE_LOAD_FLOATS_UNALIGNED(localy1);
            const SSE_FLOATS m_z2 = SSE_LOAD_FLOATS_UNALIGNED(localz1);
            const SSE_FLOATS m_d2 = SSE_LOAD_FLOATS_UNALIGNED(locald1);

            localx1 += SSE_NVEC;
            localy1 += SSE_NVEC;
            localz1 += SSE_NVEC;
            locald1 += SSE_NVEC;

            for(int w = 0; w < pair.num_weights; w++){
                pair.weights1[w].s = SSE_LOAD_FLOATS_UNALIGNED(local_w1.weights[w]);
                local_w1.weights[w] += SSE_NVEC;
            }

            union float4_weights{
                SSE_FLOATS m_weights;
                double weights[SSE_NVEC];
            };
            union float4_weights union_mweight;

            const SSE_FLOATS m_perpx = SSE_SUBTRACT_FLOATS(m_xpos, m_x2);
            const SSE_FLOATS m_perpy = SSE_SUBTRACT_FLOATS(m_ypos, m_y2);
            const SSE_FLOATS m_perpz = SSE_SUBTRACT_FLOATS(m_zpos, m_z2);

            const SSE_FLOATS m_parx = SSE_ADD_FLOATS(m_x2, m_xpos);
            const SSE_FLOATS m_pary = SSE_ADD_FLOATS(m_y2, m_ypos);
            const SSE_FLOATS m_parz = SSE_ADD_FLOATS(m_z2, m_zpos);

            SSE_FLOATS m_sqr_s, m_sqr_mu;
            {
                const SSE_FLOATS m_s_dot_l =  SSE_SUBTRACT_FLOATS(SSE_SQUARE_FLOAT(m_d2), SSE_SQUARE_FLOAT(m_dpos));

                const SSE_FLOATS m_sqr_s_dot_l = SSE_SQUARE_FLOAT(m_s_dot_l);
                const SSE_FLOATS m_sqr_perpx = SSE_SQUARE_FLOAT(m_perpx);
                const SSE_FLOATS m_sqr_perpy = SSE_SQUARE_FLOAT(m_perpy);
                const SSE_FLOATS m_sqr_perpz = SSE_SQUARE_FLOAT(m_perpz);
                m_sqr_s = SSE_ADD_FLOATS(m_sqr_perpx, SSE_ADD_FLOATS(m_sqr_perpy, m_sqr_perpz));//3-d separation

                const SSE_FLOATS m_mask_3d_sep = SSE_COMPARE_FLOATS_LT(m_sqr_s, m_sqr_smax);
                const SSE_FLOATS m_sqr_norm_l = SSE_ADD_FLOATS(SSE_SQUARE_FLOAT(m_parx), SSE_ADD_FLOATS(SSE_SQUARE_FLOAT(m_pary), SSE_SQUARE_FLOAT(m_parz)));

                if(SSE_TEST_COMPARISON(m_mask_3d_sep)==0) {
                    continue;
                }

                // \mu^2 = \pi^2 / s^2
                const SSE_FLOATS m_sqr_norm_l_norm_s = SSE_MULTIPLY_FLOATS(m_sqr_norm_l, m_sqr_s);
                m_sqr_mu = SSE_DIVIDE_FLOATS(m_sqr_s_dot_l,m_sqr_norm_l_norm_s);
            }


            const SSE_FLOATS m_mu = SSE_SQRT_FLOAT(m_sqr_mu);

            SSE_FLOATS m_mask_left;
            //Do the mask filters in a separate scope
            {
                const SSE_FLOATS m_mask_mumax = SSE_COMPARE_FLOATS_LT(m_sqr_mu,m_sqr_mumax);
                const SSE_FLOATS m_smax_mask = SSE_COMPARE_FLOATS_LT(m_sqr_s, m_sqr_smax);
                const SSE_FLOATS m_smin_mask = SSE_COMPARE_FLOATS_GE(m_sqr_s, m_sqr_smin);
                const SSE_FLOATS m_s_mask = SSE_BITWISE_AND(m_smax_mask,m_smin_mask);

                m_mask_left = SSE_BITWISE_AND(m_mask_mumax, m_s_mask);
                if(SSE_TEST_COMPARISON(m_mask_left)==0) {
                    continue;
                }

                m_sqr_s = SSE_BLEND_FLOATS_WITH_MASK(m_zero,m_sqr_s,m_mask_left);
                m_sqr_mu  = SSE_BLEND_FLOATS_WITH_MASK(m_sqr_mumax,m_sqr_mu,m_mask_left);
            }
            union float8 union_msep;
            if(need_savg) {
                union_msep.m_sep = SSE_SQRT_FLOAT(m_sqr_s);
            }
            if(need_weightavg){
                pair.dx.s = m_perpx;
                pair.dy.s = m_perpy;
                pair.dz.s = m_perpz;

                pair.parx.s = m_parx;
                pair.pary.s = m_pary;
                pair.parz.s = m_parz;

                union_mweight.m_weights = sse_weight_func(&pair);
            }

            const SSE_FLOATS m_mask = m_mask_left;
            SSE_FLOATS m_sbin = SSE_SET_FLOAT((double) 0);
            for(int kbin=nsbin-1;kbin>=1;kbin--) {
                const SSE_FLOATS m_mask_low = SSE_COMPARE_FLOATS_GE(m_sqr_s,m_supp_sqr[kbin-1]);
                const SSE_FLOATS m_bin_mask = SSE_BITWISE_AND(m_mask_low,m_mask_left);
                m_sbin = SSE_BLEND_FLOATS_WITH_MASK(m_sbin,m_kbin[kbin], m_bin_mask);
                m_mask_left = SSE_COMPARE_FLOATS_LT(m_sqr_s, m_supp_sqr[kbin-1]);
                if(SSE_TEST_COMPARISON(m_mask_left) == 0) {
                    break;
                }
            }

            /* Compute the 1-D index to the [sbin, mubin] := sbin*(nmu_bins+1) + mubin */
            const SSE_FLOATS m_tmp2 = SSE_MULTIPLY_FLOATS(m_mu,m_inv_dmu);
            const SSE_FLOATS m_mubin = SSE_BLEND_FLOATS_WITH_MASK(m_nmu_bins, m_tmp2, m_mask);
            const SSE_FLOATS m_nmu_bins_p1 = SSE_ADD_FLOATS(m_nmu_bins,m_one);
            const SSE_FLOATS m_binproduct = SSE_ADD_FLOATS(SSE_MULTIPLY_FLOATS(m_sbin,m_nmu_bins_p1),m_mubin);
            union int8 union_finalbin;
            union_finalbin.m_ibin = SSE_TRUNCATE_FLOAT_TO_INT(m_binproduct);

#if  __INTEL_COMPILER
#pragma unroll(SSE_NVEC)
#endif
            for(int jj=0;jj<SSE_NVEC;jj++) {
                const int ibin=union_finalbin.ibin[jj];
                const double weight = union_mweight.weights[jj];

                npairs[ibin]++;
                if(need_savg) {
                    savg[ibin] += union_msep.sep[jj];
                }
                if(need_weightavg){
                    //const double weight = union_mweight.weights[jj];
                    weightavg[ibin] += weight;
                }
                //TODO: generalize
                //projpairs[ibin] += weight;
            }

//            //TODO: where should this go?? smarter way within sse shit?
//            int ins = -1;
//            for(int p=0;p<nprojbins-1;p++){
//                if (m_sqr_s>=m_supp_sqr[p] && m_sqr_s<m_supp_sqr[p+1]){
//                    ins = p;
//                    break;
//                }
//            }
//            if (ins>=0 && ins<nprojbins){
//                projpairs[ins] += 1.0;
//            }
//            //outer product to get tensor value
//            for(int p1=0;p1<nprojbins-1;p1++){
//                for(int p2=0;p2<nprojbins-1;p2++){
//                    projpairs_tensor[p1*nprojbins+p2] += projpairs[p1]*projpairs[p2];
//                }
//            }
        }//SSE j loop

        //Take care of the remainder
        for(;j<N1;j++) {
            const double parx = xpos + *localx1;
            const double pary = ypos + *localy1;
            const double parz = zpos + *localz1;

            const double perpx = xpos - *localx1;
            const double perpy = ypos - *localy1;
            const double perpz = zpos - *localz1;

            //parx*perpx + pary*perpy + parz*perpz == (x1^2 + y1^2 + z1^2) - (x2^2 + y2^2 + z2^2) == d1^2 - d2^2
            const double s_dot_l = dpos*dpos - (*locald1) * (*locald1);
            localx1++;localy1++;localz1++;locald1++;

            for(int w = 0; w < pair.num_weights; w++){
                pair.weights1[w].d = *local_w1.weights[w]++;
            }

            const double sqr_s = perpx*perpx + perpy*perpy + perpz*perpz;
            if(sqr_s >= sqr_smax || sqr_s < sqr_smin) continue;

            const double norm_l = (parx*parx + pary*pary + parz*parz);
            const double sqr_s_dot_l = s_dot_l * s_dot_l;
            const double sqr_mu = sqr_s_dot_l/(norm_l * sqr_s);
            const int mubin  = (sqr_mu >= sqr_mumax) ? nmu_bins:(int) (SQRT(sqr_mu)*inv_dmu);
            double s, pairweight;
            if(need_savg) {
                s = SQRT(sqr_s);
            }
            if(need_weightavg){
                pair.dx.d = perpx;
                pair.dy.d = perpy;
                pair.dz.d = perpz;

                pair.parx.d = parx;
                pair.pary.d = pary;
                pair.parz.d = parz;

                pairweight = fallback_weight_func(&pair);
            }


            for(int kbin=nsbin-1;kbin>=1;kbin--) {
                if(sqr_s >= supp_sqr[kbin-1]) {
                    const int ibin = kbin*(nmu_bins+1) + mubin;
                    npairs[ibin]++;
                    if(need_savg){
                        savg[ibin] += s;
                    }
                    if(need_weightavg){
                        weightavg[ibin] += pairweight;
                    }
                    //TODO: generalize
                    //projpairs[ibin] += pairweight;
                    break;
                }
            }
//            //TODO: where should this go?? smarter way within sse shit?
//            int ins = -1;
//            for(int p=0;p<nprojbins-1;p++){
//                if (sqr_s>=supp_sqr[p] && sqr_s<supp_sqr[p+1]){
//                    ins = p;
//                    break;
//                }
//            }
//            if (ins>=0 && ins<nprojbins){
//                projpairs[ins] += 1.0;
//            }
//            //outer product to get tensor value
//            for(int p1=0;p1<nprojbins-1;p1++){
//                for(int p2=0;p2<nprojbins-1;p2++){
//                    projpairs_tensor[p1*nprojbins+p2] += projpairs[p1]*projpairs[p2];
//                }
//            }
        }//remainder jloop
    }//i-loop

    for(int i=0;i<totnbins;i++) {
        src_npairs[i] += npairs[i];
        if(need_savg) {
            src_savg[i] += savg[i];
        }
        if(need_weightavg) {
            src_weightavg[i] += weightavg[i];
        }
    }
    for(int i=0;i<nprojbins;i++) {
        src_projpairs[i] += projpairs[i];
        for(int j=0;j<nprojbins;j++) {
            src_projpairs_tensor[i*nprojbins+j] += projpairs_tensor[i*nprojbins+j];
        }
    }


    return EXIT_SUCCESS;
}
#endif //SSE4.2 defined



static inline int countpairs_s_mu_mocks_fallback_double(const int64_t N0, double *x0, double *y0, double *z0, double *d0, const weight_struct_double *weights0,
                                                        const int64_t N1, double *x1, double *y1, double *z1, double *d1, const weight_struct_double *weights1,
                                                        const int same_cell,
                                                        const int fast_divide,
                                                        const double smax, const double smin, const int nsbin,
                                                        const int nmu_bins, const double *supp_sqr, const double mu_max,
                                                        double *src_savg, uint64_t *src_npairs, double *src_projpairs,
                                                        double *src_projpairs_tensor,
                                                        double *src_weightavg, const weight_method_t weight_method)
{
    if(N0 == 0 || N1 == 0) {
        return EXIT_SUCCESS;
    }

    if(src_npairs == NULL) {
        return EXIT_FAILURE;
    }

    const int32_t need_savg = src_savg != NULL;
    const int32_t need_weightavg = src_weightavg != NULL;

    (void) fast_divide;//unused parameter but required to keep the same function signature amongst the kernels

    const double sqr_smax  = smax*smax;
    const double sqr_smin  = smin*smin;
    const double sqr_mumax = mu_max*mu_max;

    /*----------------- FALLBACK CODE --------------------*/
    const int64_t totnbins = (nmu_bins+1)*(nsbin+1);
    //printf("%d\n", nsbin);
    const int64_t nprojbins = nsbin-1;
    uint64_t npairs[totnbins];
    double savg[totnbins], weightavg[totnbins], projpairs[nprojbins];
    double projpairs_tensor[nprojbins*nprojbins];
    for(int i=0;i<totnbins;i++) {
        npairs[i] = ZERO;
        if(need_savg) {
            savg[i]=ZERO;
        }
        if(need_weightavg){
            weightavg[i]=ZERO;
        }
    }
    for(int i=0;i<nprojbins;i++) {
        projpairs[i] = ZERO;
        for(int j=0;j<nprojbins;j++) {
            projpairs_tensor[i*nprojbins+j] = ZERO;
        }
    }


    // A copy whose pointers we can advance
    weight_struct_double local_w0 = {.weights={NULL}, .num_weights=0},
                         local_w1 = {.weights={NULL}, .num_weights=0};
    pair_struct_double pair = {.num_weights=0};
    weight_func_t_double weight_func = NULL;
    if(need_weightavg){
        // Same particle list, new copy of num_weights pointers into that list
        local_w0 = *weights0;
        local_w1 = *weights1;
        pair.num_weights = local_w0.num_weights;
        weight_func = get_weight_func_by_method_double(weight_method);
    }

    const double dmu = mu_max/(double) nmu_bins;
    const double inv_dmu = 1.0/dmu;

    int64_t nleft=N1, n_off = 0;
    for(int64_t i=0;i<N0;i++) {
        const double xpos = *x0++;
        const double ypos = *y0++;
        const double zpos = *z0++;
        const double dpos = *d0++;//d is the co-moving distance
        for(int w = 0; w < pair.num_weights; w++){
            pair.weights0[w].d = *local_w0.weights[w]++;
        }

        /* If in the same cell, unique pairs are guaranteed by not including the current particle */
        if(same_cell == 1) {
            d1++; n_off++;
            nleft--;
        } else {
            /* For a different cell, all pairs are unique pairs, since two cells are only opened for pairs once (accounted for in the assign_ngb_cells function)*/
            while(nleft > 0) {
                /*Particles are sorted on 'd', in increasing order */
                const double dz = *d1 - dpos;
                if(dz > -smax) break;
                d1++; n_off++;
                nleft--;
            }
            /*If no particle in the second cell satisfies distance constraints on 'dz' for the current 'i'th particle in first cell,
              then there can be no more pairs from any particles in the first cell (since the first cell is also sorted in increasing order in 'd')
             */
            if(nleft == 0) {
                i=N0;
                break;
            }
        }

        double *localx1 = x1 + n_off;
        double *localy1 = y1 + n_off;
        double *localz1 = z1 + n_off;
        double *locald1 = d1;
        for(int w = 0; w < pair.num_weights; w++){
            local_w1.weights[w] = weights1->weights[w] + n_off;
        }

        for(int64_t j=0;j<nleft;j++){
            const double parx = xpos + *localx1;
            const double pary = ypos + *localy1;
            const double parz = zpos + *localz1;

            const double perpx = xpos - *localx1;
            const double perpy = ypos - *localy1;
            const double perpz = zpos - *localz1;

            //parx*perpx + pary*perpy + parz*perpz == (x1^2 + y1^2 + z1^2) - (x2^2 + y2^2 + z2^2) == d1^2 - d2^2
            const double s_dot_l = dpos*dpos - (*locald1) * (*locald1);
            localx1++;localy1++;localz1++;locald1++;

            for(int w = 0; w < pair.num_weights; w++){
                pair.weights1[w].d = *local_w1.weights[w]++;
            }

            const double sqr_s = perpx*perpx + perpy*perpy + perpz*perpz;
            if(sqr_s >= sqr_smax || sqr_s < sqr_smin) continue;

            const double sqr_l = (parx*parx + pary*pary + parz*parz);
            const double sqr_s_dot_l = s_dot_l * s_dot_l;
            const double sqr_mu = sqr_s_dot_l/(sqr_l * sqr_s);
            const int mubin  = (sqr_mu >= sqr_mumax) ? nmu_bins:(int) (SQRT(sqr_mu)*inv_dmu);
            double s, pairweight;
            if(need_savg) {
                s = SQRT(sqr_s);
            }
            if(need_weightavg){
                pair.dx.d = perpx;
                pair.dy.d = perpy;
                pair.dz.d = perpz;
                
                pair.parx.d = parx;
                pair.pary.d = pary;
                pair.parz.d = parz;

                pairweight = weight_func(&pair);
            }

            for(int kbin=nsbin-1;kbin>=1;kbin--) {
                if(sqr_s >= supp_sqr[kbin-1]) {
                    const int ibin = kbin*(nmu_bins+1) + mubin;
                    npairs[ibin]++;
                    if(need_savg) {
                        savg[ibin]+=s;
                    }
                    if(need_weightavg){
                        weightavg[ibin] += pairweight;
                    }
                    // TODO: generalize
                    //projpairs[ibin] += pairweight;
                    break;
                }

            }//finding kbin

//            weight_struct_double local_w0 = {.weights={NULL}, .num_weights=0},
//                                 local_w1 = {.weights={NULL}, .num_weights=0};
//            pair_struct_double pair = {.num_weights=0};
//            weight_func_t_double weight_func = NULL;
//            if(need_weightavg){
//                // Same particle list, new copy of num_weights pointers into that list
//                local_w0 = *weights0;
//                local_w1 = *weights1;
//                pair.num_weights = local_w0.num_weights;
//                weight_func = get_weight_func_by_method_double(weight_method);
//            }

            //proj_struct_double projdata = {.nprojbins=nprojbins, .sqr_s=sqr_s, .supp_sqr=supp_sqr};
            //proj_func = tophat_double;

            //number of edges vs bins
            int nsbins = nsbin-1;
            proj_struct_double projdata = {.nprojbins=nprojbins, .nsbins=nsbins, .sqr_s=sqr_s, .supp_sqr=supp_sqr};

            double u[nprojbins];
            //u = proj_func(projdata);
            tophat_double(&projdata, u);

            //#magic - TODO, find better solution (was getting e-310s)
            double tiny = 1.0e-16;
            //outer product to get tensor value
            for(int p1=0;p1<nprojbins;p1++){
                projpairs[p1] += u[p1];
                for(int p2=0;p2<nprojbins;p2++){
                    double uval = u[p1]*u[p2];
                    // if value is basically zero, don't bother adding (and breaks ZERO)
                    if (uval>tiny || uval<-tiny){
                        projpairs_tensor[p1*nprojbins+p2] += uval;
                    }
                }
            }
//            printf("pt[0] %f\n", projpairs_tensor[0]);
//            printf("pt[last] %f\n", projpairs_tensor[nprojbins*nprojbins-1]);

        }//j loop over second set of particles
    }//i loop over first set of particles

    for(int i=0;i<totnbins;i++) {
        src_npairs[i] += npairs[i];
        if(need_savg) {
            src_savg[i] += savg[i];
        }
        if(need_weightavg){
            src_weightavg[i] += weightavg[i];
        }
    }
    for(int i=0;i<nprojbins;i++) {
        src_projpairs[i] += projpairs[i];
        for(int j=0;j<nprojbins;j++) {
            src_projpairs_tensor[i*nprojbins+j] += projpairs_tensor[i*nprojbins+j];
        }
    }

    return EXIT_SUCCESS;
}//end of fallback code
