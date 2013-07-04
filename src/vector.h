/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Have I already read this file? */
#ifndef VEC_MACRO

    /* Include the header file with the intrinsics. */
    #include <immintrin.h>
    
    /* Define the vector macro. */
    #define VEC_MACRO(elcount, type)  __attribute__((vector_size((elcount)*sizeof(type)))) type

    /* So what will the vector size be? */
    #ifdef NO__AVX__
        #define VECTORIZE
        #define VEC_SIZE 8
        #define VEC_FLOAT __m256
        #define VEC_DBL __m256d
        #define VEC_INT __m256i
        #define vec_load(a) _mm256_load_ps(a)
        #define vec_set1(a) _mm256_set1_ps(a)
        #define vec_set(a,b,c,d,e,f,g,h) _mm256_set_ps(h,g,f,e,d,c,b,a)
        #define vec_dbl_set(a,b,c,d) _mm256_set_pd(d,c,b,a)
        #define vec_sqrt(a) _mm256_sqrt_ps(a)
        #define vec_rcp(a) _mm256_rcp_ps(a)
        #define vec_rsqrt(a) _mm256_rsqrt_ps(a)
        #define vec_ftoi(a) _mm256_cvttps_epi32(a)
        #define vec_fmin(a,b) _mm256_min_ps(a,b)
        #define vec_fmax(a,b) _mm256_max_ps(a,b)
        #define vec_todbl_lo(a) _mm256_cvtps_pd(_mm256_extract128_ps(a,0))
        #define vec_todbl_hi(a) _mm256_cvtps_pd(_mm256_extract128_ps(a,1))
        #define vec_dbl_tofloat(a,b) _mm256_insertf128( _mm256_castps128_ps256(a) , b , 1 )
        #define vec_dbl_load(a) _mm256_load_pd(a)
        #define vec_dbl_set1(a) _mm256_set1_pd(a)
        #define vec_dbl_sqrt(a) _mm256_sqrt_pd(a)
        #define vec_dbl_rcp(a) _mm256_rcp_pd(a)
        #define vec_dbl_rsqrt(a) _mm256_rsqrt_pd(a)
        #define vec_dbl_ftoi(a) _mm256_cvttpd_epi32(a)
        #define vec_dbl_fmin(a,b) _mm256_min_pd(a,b)
        #define vec_dbl_fmax(a,b) _mm256_max_pd(a,b)
    #elif defined( NO__SSE2__ )
        #define VECTORIZE
        #define VEC_SIZE 4
        #define VEC_FLOAT __m128
        #define VEC_DBL __m128d
        #define VEC_INT __m128i
        #define vec_load(a) _mm_load_ps(a)
        #define vec_set1(a) _mm_set1_ps(a)
        #define vec_set(a,b,c,d) _mm_set_ps(d,c,b,a)
        #define vec_dbl_set(a,b) _mm_set_pd(b,a)
        #define vec_sqrt(a) _mm_sqrt_ps(a)
        #define vec_rcp(a) _mm_rcp_ps(a)
        #define vec_rsqrt(a) _mm_rsqrt_ps(a)
        #define vec_ftoi(a) _mm_cvttps_epi32(a)
        #define vec_fmin(a,b) _mm_min_ps(a,b)
        #define vec_fmax(a,b) _mm_max_ps(a,b)
        #define vec_todbl_lo(a) _mm_cvtps_pd(a)
        #define vec_todbl_hi(a) _mm_cvtps_pd(_mm_movehl_ps(a,a))
        #define vec_dbl_tofloat(a,b) _mm_movelh_ps( _mm_cvtpd_ps(a) , _mm_cvtpd_ps(b) )
        #define vec_dbl_load(a) _mm_load_pd(a)
        #define vec_dbl_set1(a) _mm_set1_pd(a)
        #define vec_dbl_sqrt(a) _mm_sqrt_pd(a)
        #define vec_dbl_rcp(a) _mm_rcp_pd(a)
        #define vec_dbl_rsqrt(a) _mm_rsqrt_pd(a)
        #define vec_dbl_ftoi(a) _mm_cvttpd_epi32(a)
        #define vec_dbl_fmin(a,b) _mm_min_pd(a,b)
        #define vec_dbl_fmax(a,b) _mm_max_pd(a,b)
    #else
        #define VEC_SIZE 4
    #endif

    /* Define the composite types for element access. */
    #ifdef VECTORIZE
    typedef union {
        VEC_FLOAT v;
        VEC_DBL vd;
        VEC_INT m;
        float f[VEC_SIZE];
        double d[VEC_SIZE/2];
        int i[VEC_SIZE];
        } vector;
    #endif

#endif
