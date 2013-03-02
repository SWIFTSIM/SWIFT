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
        #define VEC_INT __m256i
        #define vec_load(a) _mm256_load_ps(a)
        #define vec_set1(a) _mm256_set1_ps(a)
        #define vec_sqrt(a) _mm256_sqrt_ps(a)
        #define vec_rcp(a) _mm256_rcp_ps(a)
        #define vec_rsqrt(a) _mm256_rsqrt_ps(a)
        #define vec_ftoi(a) _mm256_cvttps_epi32(a)
        #define vec_fmin(a,b) _mm256_min_ps(a,b)
        #define vec_fmax(a,b) _mm256_max_ps(a,b)
    #elif defined( NO__SSE2__ )
        #define VECTORIZE
        #define VEC_SIZE 4
        #define VEC_FLOAT __m128
        #define VEC_INT __m128i
        #define vec_load(a) _mm_load_ps(a)
        #define vec_set1(a) _mm_set1_ps(a)
        #define vec_sqrt(a) _mm_sqrt_ps(a)
        #define vec_rcp(a) _mm_rcp_ps(a)
        #define vec_rsqrt(a) _mm_rsqrt_ps(a)
        #define vec_ftoi(a) _mm_cvttps_epi32(a)
        #define vec_fmin(a,b) _mm_min_ps(a,b)
        #define vec_fmax(a,b) _mm_max_ps(a,b)
    #else
        #define VEC_SIZE 4
    #endif
    // #ifdef __AVX__
    //     #define VEC_SIZE 8
    //     #define VEC_FLOAT VEC_MACRO(8,float)
    //     #define VEC_DOUBLE VEC_MACRO(4,double)
    //     #define VECTORIZE
    // #elif defined(__SSE2__)
    //     #define VEC_SIZE 4
    //     #define VEC_FLOAT VEC_MACRO(4,float)
    //     #define VEC_DOUBLE VEC_MACRO(2,double)
    //     #define VECTORIZE
    // #endif

    /* Define the composite types for element access. */
    #ifdef VECTORIZE
    typedef union {
        // VEC_MACRO(VEC_SIZE,float) v;
        // VEC_MACRO(VEC_SIZE,int) m;
        VEC_FLOAT v;
        VEC_INT m;
        float f[VEC_SIZE];
        int i[VEC_SIZE];
        } vector;
    #endif

#endif
