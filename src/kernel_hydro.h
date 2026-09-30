/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_KERNEL_HYDRO_H
#define SWIFT_KERNEL_HYDRO_H

/**
 * @file kernel_hydro.h
 * @brief Kernel functions for SPH (scalar and vector version).
 *
 * All constants and kernel coefficients are taken from table 1 of
 * Dehnen & Aly, MNRAS, 425, pp. 1062-1082 (2012).
 *
 * The numerical tables in this file are GENERATED (gen_kernel_hydro.py) from
 * exact rational expressions of the kernels; every constant is rounded once.
 *
 * kernel_coeffs: monomials in x = u / gamma on kernel_ivals uniform branches
 *   (to be multiplied by kernel_constant / gamma^d), as in the original
 *   file. Used by the hand-written SIMD code.
 * kernel_poly_coeffs: on kernel_poly_ivals uniform sub-intervals of
 *   [0, gamma), W(u) = sum_k c_k t^k with t = u - kernel_poly_origin[i] and
 *   the normalisation folded in (no rescaling of u, no final multiply).
 *   Origins are chosen per sub-interval to minimise the condition number;
 *   the last one is the support edge, so W -> 0 without cancellation.
 *   Used by the scalar kernel_deval(), kernel_eval(), kernel_eval_dWdx()
 *   and (double tables) kernel_eval_double().
 */

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "dimension.h"
#include "error.h"
#include "inline.h"
#include "minmax.h"
#include "vector.h"

/* ------------------------------------------------------------------------- */
#if defined(CUBIC_SPLINE_KERNEL)

#define kernel_name "Cubic spline (M4)"
#define kernel_degree 3 /*!< Degree of the polynomial (kernel_coeffs) */
#define kernel_ivals 2  /*!< Number of branches (kernel_coeffs) */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(1.8257418583505538)) /* sqrt(10/3) */
#define kernel_constant ((float)(16. * M_1_PI))
#define kernel_poly_degree 3
#define kernel_poly_ivals 2
#define kernel_poly_ivals_over_gamma 1.09544504f
#define kernel_poly_ivals_over_gamma_d 1.0954450977651105
#define kernel_poly_root 0.418429196f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    3.f, -3.f, 0.f, (float)(1. / 2.),
    -1.f, 3.f, -3.f, 1.f,
    0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.82574189f, 1.82574189f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.8257418870925903, 1.8257418870925903};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.412529588f, -0.753172517f, 0.f, 0.418429196f,
        -0.137509853f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.4125295735283098, -0.7531725220550779, 0.0, 0.4184291920938659,
        -0.13750985784276995, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(1.778001778002667)) /* sqrt(98/31) */
#define kernel_constant ((float)(80. * M_1_PI / 7.))
#define kernel_poly_degree 3
#define kernel_poly_ivals 2
#define kernel_poly_ivals_over_gamma 1.12485826f
#define kernel_poly_ivals_over_gamma_d 1.1248582631130086
#define kernel_poly_root 0.57537061f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    3.f, -3.f, 0.f, (float)(1. / 2.),
    -1.f, 3.f, -3.f, 1.f,
    0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.77800179f, 1.77800179f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.7780017852783203, 1.7780017852783203};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.614189446f, -1.09202993f, 0.f, 0.57537061f,
        -0.204729825f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.6141894687031897, -1.0920299718534143, 0.0, 0.5753706350402474,
        -0.20472982290106323, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_1D)
#define kernel_gamma ((float)(1.7320508075688772)) /* sqrt(3) */
#define kernel_constant ((float)(8. / 3.))
#define kernel_poly_degree 3
#define kernel_poly_ivals 2
#define kernel_poly_ivals_over_gamma 1.15470052f
#define kernel_poly_ivals_over_gamma_d 1.1547005591040844
#define kernel_poly_root 0.769800365f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    3.f, -3.f, 0.f, (float)(1. / 2.),
    -1.f, 3.f, -3.f, 1.f,
    0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.73205078f, 1.73205078f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.7320507764816284, 1.7320507764816284};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.888888955f, -1.53960085f, 0.f, 0.769800365f,
        -0.296296328f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.888888952704826, -1.5396008007383353, 0.0, 0.7698003727360563,
        -0.29629631756827535, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
};
#endif

/* ------------------------------------------------------------------------- */
#elif defined(QUARTIC_SPLINE_KERNEL)

#define kernel_name "Quartic spline (M5)"
#define kernel_degree 4 /*!< Degree of the polynomial (kernel_coeffs) */
#define kernel_ivals 5  /*!< Number of branches (kernel_coeffs) */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(2.0189321327181204)) /* sqrt(375/92) */
#define kernel_constant ((float)(15625. * M_1_PI / 512.))
#define kernel_poly_degree 4
#define kernel_poly_ivals 5
#define kernel_poly_ivals_over_gamma 2.47655678f
#define kernel_poly_ivals_over_gamma_d 2.4765567845593095
#define kernel_poly_root 0.434393048f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    6.f, 0.f, (float)(-12. / 5.), 0.f, (float)(46. / 125.),
    -4.f, 8.f, (float)(-24. / 5.), (float)(8. / 25.), (float)(44. / 125.),
    -4.f, 8.f, (float)(-24. / 5.), (float)(8. / 25.), (float)(44. / 125.),
    1.f, -4.f, 6.f, -4.f, 1.f,
    1.f, -4.f, 6.f, -4.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 0.807572842f, 1.21135926f, 1.61514568f, 2.0189321f, 2.0189321f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 0.8075728416442871, 1.2113592624664307, 1.6151456832885742, 2.0189321041107178, 2.0189321041107178};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.426284403f, 0.f, -0.695028901f, 0.f, 0.434393048f,
        -0.284189612f, 0.22950381f, 0.27801156f, -0.411610067f, 0.143538579f,
        -0.284189612f, -0.22950381f, 0.27801156f, -0.149676397f, 0.0302186478f,
        0.0710474029f, -0.114751905f, 0.0695028901f, -0.0187095497f, 0.00188866549f,
        0.0710474029f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.42628440457591743, 0.0, -0.6950289008076896, 0.0, 0.4343930506944791,
        -0.28418960305061164, 0.2295038053013444, 0.27801156032307583, -0.4116100739301256, 0.14353857327295833,
        -0.28418960305061164, -0.2295038053013444, 0.27801156032307583, -0.14967639052004567, 0.03021864700483333,
        0.07104740076265291, -0.1147519026506722, 0.06950289008076896, -0.018709548815005708, 0.001888665437802083,
        0.07104740076265291, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(1.9771727313078125)) /* sqrt(38150/9759) */
#define kernel_constant ((float)(46875. * M_1_PI / 2398.))
#define kernel_poly_degree 4
#define kernel_poly_ivals 5
#define kernel_poly_ivals_over_gamma 2.52886343f
#define kernel_poly_ivals_over_gamma_d 2.5288635222320974
#define kernel_poly_root 0.585734546f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    6.f, 0.f, (float)(-12. / 5.), 0.f, (float)(46. / 125.),
    -4.f, 8.f, (float)(-24. / 5.), (float)(8. / 25.), (float)(44. / 125.),
    -4.f, 8.f, (float)(-24. / 5.), (float)(8. / 25.), (float)(44. / 125.),
    1.f, -4.f, 6.f, -4.f, 1.f,
    1.f, -4.f, 6.f, -4.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 0.790869117f, 1.18630362f, 1.58173823f, 1.97717273f, 1.97717273f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 0.7908691167831421, 1.1863036155700684, 1.5817382335662842, 1.9771727323532104, 1.9771727323532104};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.624921978f, 0.f, -0.977181017f, 0.f, 0.585734546f,
        -0.416614652f, 0.329487622f, 0.390872419f, -0.566736281f, 0.19354704f,
        -0.416614652f, -0.329487622f, 0.390872419f, -0.20608595f, 0.0407467559f,
        0.104153663f, -0.164743811f, 0.0977180749f, -0.0257607326f, 0.00254667061f,
        0.104153663f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.6249219876340327, 0.0, -0.9771810166389667, 0.0, 0.5857345246825584,
        -0.41661465842268847, 0.3294876172813246, 0.39087243022237894, -0.5667363084921317, 0.19354704681784207,
        -0.41661465842268847, -0.3294876172813246, 0.39087243022237894, -0.20608595577656974, 0.040746754456606346,
        0.10415366460567212, -0.1647438086406623, 0.09771807809710438, -0.025760732823166793, 0.00254667061807822,
        0.10415366460567212, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_1D)
#define kernel_gamma ((float)(1.9364916731037085)) /* sqrt(15/4) */
#define kernel_constant ((float)(3125. / 768.))
#define kernel_poly_degree 4
#define kernel_poly_ivals 5
#define kernel_poly_ivals_over_gamma 2.58198881f
#define kernel_poly_ivals_over_gamma_d 2.581988824504585
#define kernel_poly_root 0.773251891f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    6.f, 0.f, (float)(-12. / 5.), 0.f, (float)(46. / 125.),
    -4.f, 8.f, (float)(-24. / 5.), (float)(8. / 25.), (float)(44. / 125.),
    -4.f, 8.f, (float)(-24. / 5.), (float)(8. / 25.), (float)(44. / 125.),
    1.f, -4.f, 6.f, -4.f, 1.f,
    1.f, -4.f, 6.f, -4.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 0.774596691f, 1.16189504f, 1.54919338f, 1.93649173f, 1.93649173f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 0.7745966911315918, 1.1618950366973877, 1.5491933822631836, 1.9364917278289795, 1.9364917278289795};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.896523774f, 0.f, -1.34478581f, 0.f, 0.773251891f,
        -0.597682536f, 0.462962925f, 0.537914336f, -0.763888836f, 0.255509317f,
        -0.597682536f, -0.462962925f, 0.537914336f, -0.277777761f, 0.0537914336f,
        0.149420634f, -0.231481463f, 0.134478584f, -0.0347222202f, 0.0033619646f,
        0.149420634f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.896523796054341, 0.0, -1.3447857700888226, 0.0, 0.7732518615052795,
        -0.5976825307028939, 0.4629629106296177, 0.5379143080355291, -0.7638888457138778, 0.2555093107582663,
        -0.5976825307028939, -0.4629629106296177, 0.5379143080355291, -0.27777776207777377, 0.05379143384384553,
        0.14942063267572347, -0.23148145531480885, 0.13447857700888227, -0.03472222025972172, 0.0033619646152403455,
        0.14942063267572347, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
};
#endif

/* ------------------------------------------------------------------------- */
#elif defined(QUINTIC_SPLINE_KERNEL)

#define kernel_name "Quintic spline (M6)"
#define kernel_degree 5 /*!< Degree of the polynomial (kernel_coeffs) */
#define kernel_ivals 3  /*!< Number of branches (kernel_coeffs) */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(2.1957751641342)) /* sqrt(135/28) */
#define kernel_constant ((float)(2187. * M_1_PI / 40.))
#define kernel_poly_degree 5
#define kernel_poly_ivals 3
#define kernel_poly_ivals_over_gamma 1.36626005f
#define kernel_poly_ivals_over_gamma_d 1.366260035968407
#define kernel_poly_root 0.446491212f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    -10.f, 10.f, 0.f, (float)(-20. / 9.), 0.f, (float)(22. / 81.),
    5.f, -15.f, (float)(50. / 3.), (float)(-70. / 9.), (float)(25. / 27.), (float)(17. / 81.),
    -1.f, 5.f, -10.f, 10.f, -5.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.46385014f, 2.19577527f, 2.19577527f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.4638501405715942, 2.195775270462036, 2.195775270462036};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        -0.322059274f, 0.707169771f, 0.f, -0.757681966f, 0.f, 0.446491212f,
        0.161029637f, 0.117861599f, -0.172531784f, 0.126280352f, -0.0462138802f, 0.00676502008f,
        -0.0322059281f, 0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        -0.3220592684918646, 0.7071697773775295, 0.0, -0.7576819777127839, 0.0, 0.44649120867951336,
        0.1610296342459323, 0.11786159756920776, -0.17253178642067926, 0.12628035018618788, -0.04621388085631803, 0.006765020149700395,
        -0.03220592684918646, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(2.158129829143054)) /* sqrt(12906/2771) */
#define kernel_constant ((float)(15309. * M_1_PI / 478.))
#define kernel_poly_degree 5
#define kernel_poly_ivals 3
#define kernel_poly_ivals_over_gamma 1.39009237f
#define kernel_poly_ivals_over_gamma_d 1.3900923932370532
#define kernel_poly_root 0.594499588f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    -10.f, 10.f, 0.f, (float)(-20. / 9.), 0.f, (float)(22. / 81.),
    5.f, -15.f, (float)(50. / 3.), (float)(-70. / 9.), (float)(25. / 27.), (float)(17. / 81.),
    -1.f, 5.f, -10.f, 10.f, -5.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.43875325f, 2.15812993f, 2.15812993f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.4387532472610474, 2.158129930496216, 2.158129930496216};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        -0.467547715f, 1.00902867f, 0.f, -1.04435027f, 0.f, 0.594499588f,
        0.233773857f, 0.168171406f, -0.241957262f, 0.174058408f, -0.0626067817f, 0.00900757127f,
        -0.04675477f, 0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        -0.4675477146338943, 1.0090287168865109, 0.0, -1.0443502821526107, 0.0, 0.5944995632618847,
        0.23377385731694714, 0.16817140636772607, -0.24195725724601458, 0.17405840920231958, -0.0626067805505772, 0.009007571628101254,
        -0.04675477146338943, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_1D)
#define kernel_gamma ((float)(2.1213203435596424)) /* sqrt(9/2) */
#define kernel_constant ((float)(243. / 40.))
#define kernel_poly_degree 5
#define kernel_poly_ivals 3
#define kernel_poly_ivals_over_gamma 1.41421366f
#define kernel_poly_ivals_over_gamma_d 1.414213626312762
#define kernel_poly_root 0.777817488f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    -10.f, 10.f, 0.f, (float)(-20. / 9.), 0.f, (float)(22. / 81.),
    5.f, -15.f, (float)(50. / 3.), (float)(-70. / 9.), (float)(25. / 27.), (float)(17. / 81.),
    -1.f, 5.f, -10.f, 10.f, -5.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.41421354f, 2.12132025f, 2.12132025f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.4142135381698608, 2.1213202476501465, 2.1213202476501465};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        -0.666666865f, 1.4142139f, 0.f, -1.41421378f, 0.f, 0.777817488f,
        0.333333433f, 0.235702381f, -0.333333343f, 0.235702246f, -0.0833333209f, 0.0117851105f,
        -0.066666685f, 0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        -0.6666668475153754, 1.4142138820714587, 0.0, -1.4142137541921045, 0.0, 0.7778174944720191,
        0.3333334237576877, 0.2357023799059775, -0.3333333561523545, 0.23570225262891595, -0.08333332213676188, 0.011785110241237268,
        -0.06666668475153754, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};
#endif

/* ------------------------------------------------------------------------- */
#elif defined(WENDLAND_C2_KERNEL)

#define kernel_name "Wendland C2"
#define kernel_degree 5 /*!< Degree of the polynomial (kernel_coeffs) */
#define kernel_ivals 1  /*!< Number of branches (kernel_coeffs) */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(1.9364916731037085)) /* sqrt(15/4) */
#define kernel_constant ((float)(21. * M_1_PI / 2.))
#define kernel_poly_degree 5
#define kernel_poly_ivals 4
#define kernel_poly_ivals_over_gamma 2.0655911f
#define kernel_poly_ivals_over_gamma_d 2.065591059603668
#define kernel_poly_root 0.460248619f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    4.f, -15.f, 20.f, -10.f, 0.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 0.968245864f, 1.45236874f, 1.93649173f, 1.93649173f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 0.9682458639144897, 1.4523687362670898, 1.9364917278289795, 1.9364917278289795};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.0676042885f, -0.490931809f, 1.26758051f, -1.22732961f, 0.f, 0.460248619f,
        0.0676042885f, -0.163643926f, 0.f, 0.306832403f, -0.297089189f, 0.086296618f,
        0.0676042885f, -2.01476489e-08f, -0.158447564f, 0.153416231f, -0.0557042435f, 0.00719138794f,
        0.0676042885f, 0.163643926f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.0676042888394359, -0.4909317978874823, 1.2675804873830907, -1.2273295640873905, 0.0, 0.4602486125460297,
        0.0676042888394359, -0.16364393262916077, 0.0, 0.3068323910218476, -0.2970891935218974, 0.08629661485238058,
        0.0676042888394359, -2.014764810783741e-08, -0.15844756092288392, 0.15341622384355558, -0.055704242073993125, 0.007191387891262731,
        0.0676042888394359, 0.16364393262916077, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(1.8973665961010275)) /* sqrt(18/5) */
#define kernel_constant ((float)(7. * M_1_PI))
#define kernel_poly_degree 5
#define kernel_poly_ivals 4
#define kernel_poly_ivals_over_gamma 2.10818505f
#define kernel_poly_ivals_over_gamma_d 2.1081850547223233
#define kernel_poly_root 0.618935883f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    4.f, -15.f, 20.f, -10.f, 0.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 0.948683321f, 1.42302501f, 1.89736664f, 1.89736664f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 0.9486833214759827, 1.4230250120162964, 1.8973666429519653, 1.8973666429519653};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.100681417f, -0.716360867f, 1.81226563f, -1.71926618f, 0.f, 0.618935883f,
        0.100681417f, -0.238786966f, 0.f, 0.429816544f, -0.407759786f, 0.116050474f,
        0.100681417f, 1.50027013e-08f, -0.226533204f, 0.214908257f, -0.0764549449f, 0.00967087038f,
        0.100681417f, 0.238786966f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.10068141970627255, -0.7163608774339807, 1.8122656442120484, -1.7192661907478977, 0.0, 0.6189358592355284,
        0.10068141970627255, -0.23878695914466025, 0.0, 0.4298165476869744, -0.40775979008501895, 0.11605047360666157,
        0.10068141970627255, 1.5002700642685973e-08, -0.22653320552650516, 0.21490825358984034, -0.07645494783141034, 0.009670870522019938,
        0.10068141970627255, 0.23878695914466025, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_1D)
#define kernel_gamma ((float)(1.620185174601965)) /* sqrt(21/8) */
#define kernel_constant ((float)(5. / 4.))
#define kernel_poly_degree 4
#define kernel_poly_ivals 4
#define kernel_poly_ivals_over_gamma 2.46885371f
#define kernel_poly_ivals_over_gamma_d 2.4688536570040185
#define kernel_poly_root 0.77151674f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    0.f, -3.f, 8.f, -6.f, 0.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 0.810092568f, 1.21513891f, 1.62018514f, 1.62018514f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 0.810092568397522, 1.2151389122009277, 1.620185136795044, 1.620185136795044};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        -0.335898489f, 1.45124733f, -1.76346695f, 0.f, 0.77151674f,
        -0.335898489f, 0.362811834f, 0.440866739f, -0.714285731f, 0.241098985f,
        -0.335898489f, -0.181405991f, 0.551083386f, -0.267857105f, 0.0391785689f,
        -0.335898489f, -0.725623667f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        -0.3358984880879215, 1.45124730099194, -1.7634669801607985, 0.0, 0.7715167678137558,
        -0.3358984880879215, 0.362811825247985, 0.4408667450401996, -0.7142857476213417, 0.24109898994179868,
        -0.3358984880879215, -0.18140599270843277, 0.5510833988623375, -0.2678570896637407, 0.03917856990001364,
        -0.3358984880879215, -0.72562365049597, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
};
#endif

/* ------------------------------------------------------------------------- */
#elif defined(WENDLAND_C4_KERNEL)

#define kernel_name "Wendland C4"
#define kernel_degree 8 /*!< Degree of the polynomial (kernel_coeffs) */
#define kernel_ivals 1  /*!< Number of branches (kernel_coeffs) */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(2.207940216581962)) /* sqrt(39/8) */
#define kernel_constant ((float)(495. * M_1_PI / 32.))
#define kernel_poly_degree 8
#define kernel_poly_ivals 4
#define kernel_poly_ivals_over_gamma 1.81164336f
#define kernel_poly_ivals_over_gamma_d 1.811643348956221
#define kernel_poly_root 0.457449853f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    (float)(35. / 3.), -64.f, 140.f, (float)(-448. / 3.), 70.f, 0.f, (float)(-28. / 3.), 0.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.10397005f, 1.65595508f, 2.2079401f, 2.2079401f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.1039700508117676, 1.6559550762176514, 2.207940101623535, 2.207940101623535};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.00944913365f, -0.114449121f, 0.552774251f, -1.30185854f, 1.34738708f, 0.f, -0.875801504f, 0.f, 0.457449853f,
        0.00944913365f, -0.0309966356f, -0.00921290368f, 0.142390788f, -0.140352815f, -0.123956248f, 0.314741164f, -0.211500332f, 0.0494379401f,
        0.00944913365f, 0.0107296044f, -0.0483677462f, 0.00254269247f, 0.0894749239f, -0.104588084f, 0.053668499f, -0.01345482f, 0.00134716521f,
        0.00944913365f, 0.0524558462f, 0.0737032294f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.009449133309073026, -0.11444911739623695, 0.5527742410193748, -1.3018585748172702, 1.3473870721788348, 0.0, -0.8758015057174228, 0.0, 0.45744984597623956,
        0.009449133309073026, -0.03099663596148084, -0.009212904016989579, 0.14239078162063892, -0.14035282001862862, -0.12395624787803225, 0.3147411661171988, -0.21150032591797502, 0.049437939083369645,
        0.009449133309073026, 0.010729604755897213, -0.048367746089195286, 0.002542692528939981, 0.08947492276187574, -0.1045880841470897, 0.053668500472429964, -0.013454819840764036, 0.001347165226339939,
        0.009449133309073026, 0.05245584547327527, 0.07370323213591663, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(2.1712405933672376)) /* sqrt(33/7) */
#define kernel_constant ((float)(9. * M_1_PI))
#define kernel_poly_degree 8
#define kernel_poly_ivals 4
#define kernel_poly_ivals_over_gamma 1.84226477f
#define kernel_poly_ivals_over_gamma_d 1.842264767274455
#define kernel_poly_root 0.607682526f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    (float)(35. / 3.), -64.f, 140.f, (float)(-448. / 3.), 70.f, 0.f, (float)(-28. / 3.), 0.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.08562028f, 1.62843037f, 2.17124057f, 2.17124057f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.0856202840805054, 1.6284303665161133, 2.1712405681610107, 2.1712405681610107};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.0143535715f, -0.170962602f, 0.812002003f, -1.88058853f, 1.91400468f, 0.f, -1.20308864f, 0.f, 0.607682526f,
        0.0143535715f, -0.0463023707f, -0.0135333668f, 0.205689371f, -0.199375495f, -0.173156857f, 0.432359993f, -0.285708815f, 0.0656740218f,
        0.0143535715f, 0.0160277374f, -0.071050182f, 0.00367304985f, 0.127101868f, -0.146101132f, 0.0737244561f, -0.0181756821f, 0.00178959349f,
        0.0143535715f, 0.0783578604f, 0.108266935f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.014353571515350723, -0.1709625971499692, 0.8120020268722712, -1.8805885249856453, 1.914004733187614, 0.0, -1.2030886614985747, 0.0, 0.6076825240965745,
        0.014353571515350723, -0.04630237006144999, -0.013533367114537854, 0.20568936992030495, -0.19937549304037644, -0.17315686351454743, 0.4323599877260503, -0.2857088181653382, 0.06567402278647876,
        0.014353571515350723, 0.01602773663849336, -0.071050184038618, 0.0036730498723872777, 0.1271018757185896, -0.14610113389384816, 0.0737244533798334, -0.018175681372817033, 0.0017895934776211103,
        0.014353571515350723, 0.07835785702706921, 0.10826693691630283, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_1D)
#error "Wendland C4 kernel not defined in 1D."
#endif

/* ------------------------------------------------------------------------- */
#elif defined(WENDLAND_C6_KERNEL)

#define kernel_name "Wendland C6"
#define kernel_degree 11 /*!< Degree of the polynomial (kernel_coeffs) */
#define kernel_ivals 1  /*!< Number of branches (kernel_coeffs) */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(2.449489742783178)) /* sqrt(6) */
#define kernel_constant ((float)(1365. * M_1_PI / 64.))
#define kernel_poly_degree 11
#define kernel_poly_ivals 4
#define kernel_poly_ivals_over_gamma 1.6329931f
#define kernel_poly_ivals_over_gamma_d 1.6329931024279474
#define kernel_poly_root 0.461929709f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    32.f, -231.f, 704.f, -1155.f, 1056.f, -462.f, 0.f, 66.f, 0.f, -11.f, 0.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.22474492f, 1.83711743f, 2.44948983f, 2.44948983f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.2247449159622192, 1.8371174335479736, 2.4494898319244385, 2.4494898319244385};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.000776057364f, -0.013722444f, 0.102439575f, -0.411673337f, 0.921956241f, -0.988016069f, 0.f, 0.846870959f, 0.f, -0.846871018f, 0.f, 0.461929709f,
        0.000776057364f, -0.00326724839f, -0.00160061836f, 0.0264647137f, -0.0288111325f, -0.0617510043f, 0.151258454f, -0.0396970771f, -0.194475174f, 0.261339098f, -0.137753263f, 0.0275172964f,
        0.000776057364f, 0.0019603495f, -0.00560216326f, -0.00808644388f, 0.0252097379f, -0.00385942729f, -0.0425414443f, 0.05975236f, -0.0395027548f, 0.0145633426f, -0.00289623532f, 0.000243613191f,
        0.000776057364f, 0.0071879467f, 0.0224086568f, 0.0235241912f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.0007760573368080279, -0.013722443510027344, 0.10243957591457366, -0.411673335263859, 0.9219562503343937, -0.9880160765445598, 0.0, 0.8468709843907398, 0.0, -0.8468710460290043, 0.0, 0.46192969509123993,
        0.0007760573368080279, -0.0032672484547684156, -0.0016006183736652134, 0.02646471440981951, -0.028811132822949802, -0.06175100478403499, 0.15125845832961105, -0.039697077393315926, -0.19447517486408972, 0.26133911186051306, -0.13775325888823206, 0.02751729628961488,
        0.0007760573368080279, 0.00196034958168389, -0.005602163139368994, -0.008086443519346325, 0.025209737364165045, -0.0038594272806792987, -0.042541442785443744, 0.05975235902305969, -0.03950275564819418, 0.01456334237597178, -0.0028962352097487228, 0.0002436131862360962,
        0.0007760573368080279, 0.007187946600490514, 0.022408657231312988, 0.02352419058650623, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(2.41522945769824)) /* sqrt(35/6) */
#define kernel_constant ((float)(78. * M_1_PI / 7.))
#define kernel_poly_degree 11
#define kernel_poly_ivals 4
#define kernel_poly_ivals_over_gamma 1.65615726f
#define kernel_poly_ivals_over_gamma_d 1.6561572729955074
#define kernel_poly_root 0.608036816f
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
    32.f, -231.f, 704.f, -1155.f, 1056.f, -462.f, 0.f, 66.f, 0.f, -11.f, 0.f, 1.f,
    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const float kernel_poly_origin[kernel_poly_ivals + 1] = {
    0.f, 1.20761478f, 1.81142211f, 2.41522956f, 2.41522956f};
static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {
    0.0, 1.207614779472351, 1.8114221096038818, 2.415229558944702, 2.415229558944702};
static const float
    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)]
    __attribute__((aligned(64))) = {
        0.00119271665f, -0.0207949411f, 0.153065309f, -0.606519163f, 1.33932161f, -1.41521144f, 0.f, 1.17934299f, 0.f, -1.14658356f, 0.f, 0.608036816f,
        0.00119271665f, -0.00495117623f, -0.00239164545f, 0.0389905162f, -0.0418538004f, -0.0884507149f, 0.213628784f, -0.0552817024f, -0.267036021f, 0.35382852f, -0.183896333f, 0.0362209417f,
        0.00119271665f, 0.00297070504f, -0.00837076083f, -0.0119137643f, 0.036622081f, -0.00552818505f, -0.060083095f, 0.0832104981f, -0.0542417094f, 0.0197174121f, -0.00386638916f, 0.000320667838f,
        0.00119271665f, 0.0108925877f, 0.0334830359f, 0.0346582383f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
};
static const double
    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {
        0.0011927166195616044, -0.02079494076522692, 0.1530653123434064, -0.6065191565030366, 1.3393215952935698, -1.4152114838251264, 0.0, 1.179343002063759, 0.0, -1.1465835703582563, 0.0, 0.6080367928344745,
        0.0011927166195616044, -0.004951176372673076, -0.002391645505365725, 0.03899051720376664, -0.041853799852924055, -0.0884507177390704, 0.21362878799327734, -0.0552817032217387, -0.26703600737990374, 0.35382852366524314, -0.18389633994454463, 0.036220941760647406,
        0.0011927166195616044, 0.002970705041597891, -0.008370761039458458, -0.011913764655169693, 0.0366220805522353, -0.00552818513861343, -0.060083094646078916, 0.08321049827637914, -0.05424170883796958, 0.019717412100809694, -0.0038663890933579065, 0.00032066782768889356,
        0.0011927166195616044, 0.010892588019880767, 0.033483037075120146, 0.034658237514459234, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};
#elif defined(HYDRO_DIMENSION_1D)
#error "Wendland C6 kernel not defined in 1D."
#endif

/* ------------------------------------------------------------------------- */
#else

#error "A kernel function must be chosen at configure time !!"

/* ------------------------------------------------------------------------- */
#endif

/* Ok, now comes the real deal. */

/* First some powers of gamma = H/h (evaluated in double, rounded once) */
#define kernel_gamma_inv ((float)(1. / (double)kernel_gamma))
#define kernel_gamma2 ((float)((double)kernel_gamma * (double)kernel_gamma))

/* define gamma^d, gamma^(d+1), 1/gamma^d and 1/gamma^(d+1) */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma_dim_d \
  ((double)kernel_gamma * (double)kernel_gamma * (double)kernel_gamma)
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma_dim_d ((double)kernel_gamma * (double)kernel_gamma)
#elif defined(HYDRO_DIMENSION_1D)
#define kernel_gamma_dim_d ((double)kernel_gamma)
#endif
#define kernel_gamma_dim ((float)(kernel_gamma_dim_d))
#define kernel_gamma_dim_plus_one \
  ((float)(kernel_gamma_dim_d * (double)kernel_gamma))
#define kernel_gamma_inv_dim ((float)(1. / kernel_gamma_dim_d))
#define kernel_gamma_inv_dim_plus_one \
  ((float)(1. / (kernel_gamma_dim_d * (double)kernel_gamma)))

/* The number of branches (floating point conversion) */
#define kernel_ivals_f ((float)(kernel_ivals))

/* Kernel self contribution (i.e. W(0,h)), identical to kernel_deval(0) */
#define kernel_root (kernel_poly_root)

/* Kernel normalisation constant (volume term) */
#define kernel_norm ((float)(hydro_dimension_unit_sphere * kernel_gamma_dim))

/* ------------------------------------------------------------------------- */

/**
 * @brief Select the sub-interval of the scalar kernel table.
 *
 * kernel_poly_ivals_over_gamma is rounded down at generation time such that
 * any u < kernel_gamma maps to a row < kernel_poly_ivals; u >= kernel_gamma
 * (and only those) map to the final all-zero row.
 */
__attribute__((always_inline, const)) INLINE static int kernel_poly_index(
    const float u) {
  const int temp = (int)(u * kernel_poly_ivals_over_gamma);
  return temp > kernel_poly_ivals ? kernel_poly_ivals : temp;
}

/**
 * @brief Computes the kernel function and its derivative.
 *
 * The kernel function needs to be mutliplied by \f$h^{-d}\f$ and the gradient
 * by \f$h^{-(d+1)}\f$, where \f$d\f$ is the dimensionality of the problem.
 *
 * Returns 0 if \f$u > \gamma = H/h\f$.
 *
 * @param u The ratio of the distance to the smoothing length \f$u = x/h\f$.
 * @param W (return) The value of the kernel function \f$W(x,h)\f$.
 * @param dW_dx (return) The norm of the gradient of \f$|\nabla W(x,h)|\f$.
 */
__attribute__((always_inline)) INLINE static void kernel_deval(
    float u, float *restrict W, float *restrict dW_dx) {

#ifdef KERNEL_HYDRO_EVAL_ALL_BRANCHES
  /* Gather-free variant: evaluate every sub-interval with constant
   * coefficients and select. Can be much faster when the compiler
   * auto-vectorises the neighbour loop (e.g. icx, clang -ffast-math). */
  const int ind = kernel_poly_index(u);
  float w = 0.f, dw_dx = 0.f;
  for (int i = 0; i < kernel_poly_ivals; i++) {
    const float t = u - kernel_poly_origin[i];
    const float *const c = &kernel_poly_coeffs[i * (kernel_poly_degree + 1)];
    float wi = c[0], di = 0.f;
    for (int k = 1; k <= kernel_poly_degree; k++) {
      di = di * t + wi;
      wi = wi * t + c[k];
    }
    w = (ind == i) ? wi : w;
    dw_dx = (ind == i) ? di : dw_dx;
  }
#else
  /* Pick the correct branch of the kernel */
  const int ind = kernel_poly_index(u);
  const float *const coeffs =
      &kernel_poly_coeffs[ind * (kernel_poly_degree + 1)];

  /* Local variable (exact near the edge of the support) */
  const float t = u - kernel_poly_origin[ind];

  /* Horner's scheme for the polynomial and its derivative */
  float w = coeffs[0];
  float dw_dx = 0.f;
  for (int k = 1; k <= kernel_poly_degree; k++) {
    dw_dx = dw_dx * t + w;
    w = w * t + coeffs[k];
  }
#endif

  /* Return everything (normalisation already in the coefficients) */
  *W = max(w, 0.f);
  *dW_dx = min(dw_dx, 0.f);
}

/**
 * @brief Computes the kernel function.
 *
 * The kernel function needs to be mutliplied by \f$h^{-d}\f$,
 * where \f$d\f$ is the dimensionality of the problem.
 *
 * Returns 0 if \f$u > \gamma = H/h\f$
 *
 * @param u The ratio of the distance to the smoothing length \f$u = x/h\f$.
 * @param W (return) The value of the kernel function \f$W(x,h)\f$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval(
    float u, float *restrict W) {

  const int ind = kernel_poly_index(u);
  const float *const coeffs =
      &kernel_poly_coeffs[ind * (kernel_poly_degree + 1)];
  const float t = u - kernel_poly_origin[ind];

  float w = coeffs[0];
  for (int k = 1; k <= kernel_poly_degree; k++) w = w * t + coeffs[k];

  *W = max(w, 0.f);
}

/**
 * @brief Computes the kernel function in double precision.
 *
 * Required for computing the projected kernel because rounding
 * error causes problems for the GSL integration function if
 * we evaluate in single precision. Uses double-precision coefficients
 * of the same kernel (same support kernel_gamma) as the float version.
 *
 * The kernel function needs to be mutliplied by \f$h^{-d}\f$,
 * where \f$d\f$ is the dimensionality of the problem.
 *
 * Returns 0 if \f$u > \gamma = H/h\f$
 *
 * @param u The ratio of the distance to the smoothing length \f$u = x/h\f$.
 * @param W (return) The value of the kernel function \f$W(x,h)\f$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_double(
    double u, double *restrict W) {

  const int temp = (int)(u * kernel_poly_ivals_over_gamma_d);
  const int ind = temp > kernel_poly_ivals ? kernel_poly_ivals : temp;
  const double *const coeffs =
      &kernel_poly_coeffs_d[ind * (kernel_poly_degree + 1)];
  const double t = u - kernel_poly_origin_d[ind];

  double w = coeffs[0];
  for (int k = 1; k <= kernel_poly_degree; k++) w = w * t + coeffs[k];

  *W = max(w, 0.);
}

/**
 * @brief Computes the kernel function derivative.
 *
 * The kernel function needs to be mutliplied by \f$h^{-d}\f$ and the gradient
 * by \f$h^{-(d+1)}\f$, where \f$d\f$ is the dimensionality of the problem.
 *
 * Returns 0 if \f$u > \gamma = H/h\f$.
 *
 * @param u The ratio of the distance to the smoothing length \f$u = x/h\f$.
 * @param dW_dx (return) The norm of the gradient of \f$|\nabla W(x,h)|\f$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_dWdx(
    float u, float *restrict dW_dx) {

  const int ind = kernel_poly_index(u);
  const float *const coeffs =
      &kernel_poly_coeffs[ind * (kernel_poly_degree + 1)];
  const float t = u - kernel_poly_origin[ind];

  float w = coeffs[0];
  float dw_dx = 0.f;
  for (int k = 1; k < kernel_poly_degree; k++) {
    dw_dx = dw_dx * t + w;
    w = w * t + coeffs[k];
  }
  dw_dx = dw_dx * t + w;

  *dW_dx = min(dw_dx, 0.f);
}


#ifdef WENDLAND_C2_KERNEL

/**
 * Computes dphi/dh for the chosen kernel.
 *
 * This corresponds to Appendix A of Price & Monaghan 2007, MNRAS, 374, 4
 * but for a Wendland-C2 kernel.
 *
 * Assumes r < H (i.e. u < kernel_gamma)
 *
 * @param u The ratio r / h.
 * @param h_inv 1 / h.
 */
__attribute__((always_inline, const)) INLINE static float potential_dh(
    const float u, const float h_inv) {

  /* Ratio of r to kernel support
   * Recall that in our definition the kernel edge is H = kernel_gamma * h. */
  const float q = fminf(u * kernel_gamma_inv, 1.f);

  /* -24 q^7 + 105 q^6 -168 q^5 + 105 q^4 - 21 q^2 */
  float dphi_dh = -24.f * q + 105.f;
  dphi_dh = dphi_dh * q - 168.f;
  dphi_dh = dphi_dh * q + 105.f;
  dphi_dh = dphi_dh * q;
  dphi_dh = dphi_dh * q - 21.f;
  dphi_dh = dphi_dh * q;
  dphi_dh = dphi_dh * q;

  return dphi_dh * h_inv * h_inv * kernel_gamma_inv * 0.25f;
}

#endif

/* -------------------------------------------------------------------------
 */

#ifdef WITH_OLD_VECTORIZATION
/**
 * @brief Computes the kernel function and its derivative (Vectorised version).
 *
 * Return 0 if $u > \\gamma = H/h$
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param w (return) The value of the kernel function $W(x,h)$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 */
__attribute__((always_inline)) INLINE static void kernel_deval_vec(
    vector *u, vector *w, vector *dw_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

  /* Load x and get the interval id. */
  vector ind;
  ind.m =
      vec_ftoi(vec_fmin(vec_mul(x.v, kernel_ivals_vec.v), kernel_ivals_vec.v));

  /* load the coefficients. */
  vector c[kernel_degree + 1];
  for (int k = 0; k < VEC_SIZE; k++)
    for (int j = 0; j < kernel_degree + 1; j++)
      c[j].f[k] = kernel_coeffs[ind.i[k] * (kernel_degree + 1) + j];

  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(c[0].v, x.v, c[1].v);
  dw_dx->v = c[0].v;

  /* And we're off! */
  for (int k = 2; k <= kernel_degree; k++) {
    dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
    w->v = vec_fma(x.v, w->v, c[k].v);
  }

  /* Return everything */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
}
#endif

#ifdef WITH_VECTORIZATION

static const vector kernel_gamma_inv_vec = FILL_VEC((float)kernel_gamma_inv);

static const vector kernel_ivals_vec = FILL_VEC((float)kernel_ivals);

static const vector kernel_constant_vec = FILL_VEC((float)kernel_constant);

static const vector kernel_gamma_inv_dim_vec =
    FILL_VEC((float)kernel_gamma_inv_dim);

static const vector kernel_gamma_inv_dim_plus_one_vec =
    FILL_VEC((float)kernel_gamma_inv_dim_plus_one);

/* Define constant vectors for the Wendland C2 and Cubic Spline kernel
 * coefficients. */
#ifdef WENDLAND_C2_KERNEL
static const vector wendland_const_c0 = FILL_VEC(4.f);
static const vector wendland_const_c1 = FILL_VEC(-15.f);
static const vector wendland_const_c2 = FILL_VEC(20.f);
static const vector wendland_const_c3 = FILL_VEC(-10.f);
static const vector wendland_const_c4 = FILL_VEC(0.f);
static const vector wendland_const_c5 = FILL_VEC(1.f);

static const vector wendland_dwdx_const_c0 = FILL_VEC(20.f);
static const vector wendland_dwdx_const_c1 = FILL_VEC(-60.f);
static const vector wendland_dwdx_const_c2 = FILL_VEC(60.f);
static const vector wendland_dwdx_const_c3 = FILL_VEC(-20.f);
#elif defined(CUBIC_SPLINE_KERNEL)
/* First region 0 < u < 0.5 */
static const vector cubic_1_const_c0 = FILL_VEC(3.f);
static const vector cubic_1_const_c1 = FILL_VEC(-3.f);
static const vector cubic_1_const_c2 = FILL_VEC(0.f);
static const vector cubic_1_const_c3 = FILL_VEC(0.5f);
static const vector cubic_1_dwdx_const_c0 = FILL_VEC(9.f);
static const vector cubic_1_dwdx_const_c1 = FILL_VEC(-6.f);
static const vector cubic_1_dwdx_const_c2 = FILL_VEC(0.f);

/* Second region 0.5 <= u < 1 */
static const vector cubic_2_const_c0 = FILL_VEC(-1.f);
static const vector cubic_2_const_c1 = FILL_VEC(3.f);
static const vector cubic_2_const_c2 = FILL_VEC(-3.f);
static const vector cubic_2_const_c3 = FILL_VEC(1.f);
static const vector cubic_2_dwdx_const_c0 = FILL_VEC(-3.f);
static const vector cubic_2_dwdx_const_c1 = FILL_VEC(6.f);
static const vector cubic_2_dwdx_const_c2 = FILL_VEC(-3.f);
static const vector cond = FILL_VEC(0.5f);
#endif

/**
 * @brief Computes the kernel function and its derivative for two particles
 * using vectors. The return value is undefined if $u > \\gamma = H/h$.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param w (return) The value of the kernel function $W(x,h)$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 */
__attribute__((always_inline)) INLINE static void kernel_deval_1_vec(
    vector *u, vector *w, vector *dw_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(wendland_const_c0.v, x.v, wendland_const_c1.v);
  dw_dx->v = wendland_const_c0.v;

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c3.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  w->v = vec_mul(x.v, w->v); /* wendland_const_c4 is zero. */

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c5.v);
#elif defined(CUBIC_SPLINE_KERNEL)
  vector w2, dw_dx2;
  mask_t mask_reg;

  /* Form a mask for one part of the kernel. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using a mask. */

  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(cubic_1_const_c0.v, x.v, cubic_1_const_c1.v);
  w2.v = vec_fma(cubic_2_const_c0.v, x.v, cubic_2_const_c1.v);
  dw_dx->v = cubic_1_const_c0.v;
  dw_dx2.v = cubic_2_const_c0.v;

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, w2.v);
  w->v = vec_mul(x.v, w->v); /* cubic_1_const_c2 is zero. */
  w2.v = vec_fma(x.v, w2.v, cubic_2_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, w2.v);
  w->v = vec_fma(x.v, w->v, cubic_1_const_c3.v);
  w2.v = vec_fma(x.v, w2.v, cubic_2_const_c3.v);

  /* Blend both kernel regions into one vector (mask out unneeded values). */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  w->v = vec_blend(mask_reg, w->v, w2.v);
  dw_dx->v = vec_blend(mask_reg, dw_dx->v, dw_dx2.v);

#else
#error \
    "Vectorisation not supported for this kernel!!! Choose a different one or configure with --disable-hand-vec."
#endif

  /* Return everyting */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
}

/**
 * @brief Computes the kernel function and its derivative for two particles
 * using interleaved vectors. The return value is undefined if $u > \\gamma =
 * H/h$.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param w (return) The value of the kernel function $W(x,h)$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 * @param u2 The ratio of the distance to the smoothing length $u = x/h$ for
 * second particle.
 * @param w2 (return) The value of the kernel function $W(x,h)$ for second
 * particle.
 * @param dw_dx2 (return) The norm of the gradient of $|\\nabla W(x,h)|$ for
 * second particle.
 */
__attribute__((always_inline)) INLINE static void kernel_deval_2_vec(
    vector *u, vector *w, vector *dw_dx, vector *u2, vector *w2,
    vector *dw_dx2) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x, x2;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);
  x2.v = vec_mul(u2->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(wendland_const_c0.v, x.v, wendland_const_c1.v);
  w2->v = vec_fma(wendland_const_c0.v, x2.v, wendland_const_c1.v);
  dw_dx->v = wendland_const_c0.v;
  dw_dx2->v = wendland_const_c0.v;

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c2.v);
  w2->v = vec_fma(x2.v, w2->v, wendland_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c3.v);
  w2->v = vec_fma(x2.v, w2->v, wendland_const_c3.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  w->v = vec_mul(x.v, w->v);    /* wendland_const_c4 is zero. */
  w2->v = vec_mul(x2.v, w2->v); /* wendland_const_c4 is zero. */

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c5.v);
  w2->v = vec_fma(x2.v, w2->v, wendland_const_c5.v);

  /* Return everything */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  w2->v = vec_mul(w2->v,
                  vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
  dw_dx2->v = vec_mul(dw_dx2->v, vec_mul(kernel_constant_vec.v,
                                         kernel_gamma_inv_dim_plus_one_vec.v));
#elif defined(CUBIC_SPLINE_KERNEL)
  vector w_2, dw_dx_2;
  vector w2_2, dw_dx2_2;
  mask_t mask_reg, mask_reg_v2;

  /* Form a mask for one part of the kernel for each vector. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v));     /* 0.5 < x < 1 */
  vec_create_mask(mask_reg_v2, vec_cmp_gte(x2.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(cubic_1_const_c0.v, x.v, cubic_1_const_c1.v);
  w2->v = vec_fma(cubic_1_const_c0.v, x2.v, cubic_1_const_c1.v);
  w_2.v = vec_fma(cubic_2_const_c0.v, x.v, cubic_2_const_c1.v);
  w2_2.v = vec_fma(cubic_2_const_c0.v, x2.v, cubic_2_const_c1.v);
  dw_dx->v = cubic_1_const_c0.v;
  dw_dx2->v = cubic_1_const_c0.v;
  dw_dx_2.v = cubic_2_const_c0.v;
  dw_dx2_2.v = cubic_2_const_c0.v;

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  dw_dx_2.v = vec_fma(dw_dx_2.v, x.v, w_2.v);
  dw_dx2_2.v = vec_fma(dw_dx2_2.v, x2.v, w2_2.v);
  w->v = vec_mul(x.v, w->v);    /* cubic_1_const_c2 is zero. */
  w2->v = vec_mul(x2.v, w2->v); /* cubic_1_const_c2 is zero. */
  w_2.v = vec_fma(x.v, w_2.v, cubic_2_const_c2.v);
  w2_2.v = vec_fma(x2.v, w2_2.v, cubic_2_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  dw_dx_2.v = vec_fma(dw_dx_2.v, x.v, w_2.v);
  dw_dx2_2.v = vec_fma(dw_dx2_2.v, x2.v, w2_2.v);
  w->v = vec_fma(x.v, w->v, cubic_1_const_c3.v);
  w2->v = vec_fma(x2.v, w2->v, cubic_1_const_c3.v);
  w_2.v = vec_fma(x.v, w_2.v, cubic_2_const_c3.v);
  w2_2.v = vec_fma(x2.v, w2_2.v, cubic_2_const_c3.v);

  /* Blend both kernel regions into one vector (mask out unneeded values). */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  w->v = vec_blend(mask_reg, w->v, w_2.v);
  w2->v = vec_blend(mask_reg_v2, w2->v, w2_2.v);
  dw_dx->v = vec_blend(mask_reg, dw_dx->v, dw_dx_2.v);
  dw_dx2->v = vec_blend(mask_reg_v2, dw_dx2->v, dw_dx2_2.v);

  /* Return everything */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  w2->v = vec_mul(w2->v,
                  vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
  dw_dx2->v = vec_mul(dw_dx2->v, vec_mul(kernel_constant_vec.v,
                                         kernel_gamma_inv_dim_plus_one_vec.v));

#endif
}

/**
 * @brief Computes the kernel function for two particles
 * using vectors. The return value is undefined if $u > \\gamma = H/h$.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param w (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_W_vec(vector *u,
                                                                    vector *w) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(wendland_const_c0.v, x.v, wendland_const_c1.v);

  /* Calculate the polynomial interleaving vector operations */
  w->v = vec_fma(x.v, w->v, wendland_const_c2.v);
  w->v = vec_fma(x.v, w->v, wendland_const_c3.v);
  w->v = vec_mul(x.v, w->v); /* wendland_const_c4 is zero.*/
  w->v = vec_fma(x.v, w->v, wendland_const_c5.v);
#elif defined(CUBIC_SPLINE_KERNEL)
  vector w2;
  mask_t mask_reg;

  /* Form a mask for each part of the kernel. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(cubic_1_const_c0.v, x.v, cubic_1_const_c1.v);
  w2.v = vec_fma(cubic_2_const_c0.v, x.v, cubic_2_const_c1.v);

  /* Calculate the polynomial interleaving vector operations. */
  w->v = vec_mul(x.v, w->v); /* cubic_1_const_c2 is zero */
  w2.v = vec_fma(x.v, w2.v, cubic_2_const_c2.v);

  w->v = vec_fma(x.v, w->v, cubic_1_const_c3.v);
  w2.v = vec_fma(x.v, w2.v, cubic_2_const_c3.v);

  /* Mask out unneeded values. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  w->v = vec_blend(mask_reg, w->v, w2.v);

#else
#error \
    "Vectorisation not supported for this kernel!!! Choose a different one or configure with --disable-hand-vec."
#endif

  /* Return everything */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
}

/**
 * @brief Computes the kernel function derivative for two particles
 * using vectors. The return value is undefined if $u > \\gamma = H/h$.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_dWdx_vec(
    vector *u, vector *dw_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(wendland_dwdx_const_c0.v, x.v, wendland_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c3.v);

  dw_dx->v = vec_mul(dw_dx->v, x.v);

#elif defined(CUBIC_SPLINE_KERNEL)
  vector dw_dx2;
  mask_t mask_reg1, mask_reg2;

  /* Form a mask for each part of the kernel. */
  vec_create_mask(mask_reg1, vec_cmp_lt(x.v, cond.v));  /* 0 < x < 0.5 */
  vec_create_mask(mask_reg2, vec_cmp_gte(x.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(cubic_1_dwdx_const_c0.v, x.v, cubic_1_dwdx_const_c1.v);
  dw_dx2.v = vec_fma(cubic_2_dwdx_const_c0.v, x.v, cubic_2_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_mul(dw_dx->v, x.v); /* cubic_1_dwdx_const_c2 is zero. */
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, cubic_2_dwdx_const_c2.v);

  /* Mask out unneeded values. */
  dw_dx->v = vec_and_mask(dw_dx->v, mask_reg1);
  dw_dx2.v = vec_and_mask(dw_dx2.v, mask_reg2);

  /* Added both dwdx and dwdx2 together to form complete result. */
  dw_dx->v = vec_add(dw_dx->v, dw_dx2.v);
#else
#error \
    "Vectorisation not supported for this kernel!!! Choose a different one or configure with --disable-hand-vec."
#endif

  /* Return everything */
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
}

/**
 * @brief Computes the kernel function derivative for two particles
 * using vectors.
 *
 * Return 0 if $u > \\gamma = H/h$
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_dWdx_force_vec(
    vector *u, vector *dw_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(wendland_dwdx_const_c0.v, x.v, wendland_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c3.v);

  dw_dx->v = vec_mul(dw_dx->v, x.v);

#elif defined(CUBIC_SPLINE_KERNEL)
  vector dw_dx2;
  mask_t mask_reg;

  /* Form a mask for each part of the kernel. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(cubic_1_dwdx_const_c0.v, x.v, cubic_1_dwdx_const_c1.v);
  dw_dx2.v = vec_fma(cubic_2_dwdx_const_c0.v, x.v, cubic_2_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_mul(dw_dx->v, x.v); /* cubic_1_dwdx_const_c2 is zero. */
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, cubic_2_dwdx_const_c2.v);

  /* Mask out unneeded values. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  dw_dx->v = vec_blend(mask_reg, dw_dx->v, dw_dx2.v);

#else
#error \
    "Vectorisation not supported for this kernel!!! Choose a different one or configure with --disable-hand-vec."
#endif

  /* Mask out result for particles that lie outside of the kernel function. */
  mask_t mask;
  vec_create_mask(mask, vec_cmp_lt(x.v, vec_set1(1.f))); /* x < 1 */

  dw_dx->v = vec_and_mask(dw_dx->v, mask);

  /* Return everything */
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
}

/**
 * @brief Computes the kernel function derivative for two particles
 * using interleaved vectors.
 *
 * Return 0 if $u > \\gamma = H/h$
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 * @param u_2 The ratio of the distance to the smoothing length $u = x/h$ for
 * second particle.
 * @param dw_dx_2 (return) The norm of the gradient of $|\\nabla W(x,h)|$ for
 * second particle.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_dWdx_force_2_vec(
    vector *u, vector *dw_dx, vector *u_2, vector *dw_dx_2) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x, x_2;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);
  x_2.v = vec_mul(u_2->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(wendland_dwdx_const_c0.v, x.v, wendland_dwdx_const_c1.v);
  dw_dx_2->v =
      vec_fma(wendland_dwdx_const_c0.v, x_2.v, wendland_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c2.v);
  dw_dx_2->v = vec_fma(dw_dx_2->v, x_2.v, wendland_dwdx_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c3.v);
  dw_dx_2->v = vec_fma(dw_dx_2->v, x_2.v, wendland_dwdx_const_c3.v);

  dw_dx->v = vec_mul(dw_dx->v, x.v);
  dw_dx_2->v = vec_mul(dw_dx_2->v, x_2.v);

#elif defined(CUBIC_SPLINE_KERNEL)
  vector dw_dx2, dw_dx2_2;
  mask_t mask_reg;
  mask_t mask_reg_v2;

  /* Form a mask for one part of the kernel. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v));      /* 0.5 < x < 1 */
  vec_create_mask(mask_reg_v2, vec_cmp_gte(x_2.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(cubic_1_dwdx_const_c0.v, x.v, cubic_1_dwdx_const_c1.v);
  dw_dx_2->v = vec_fma(cubic_1_dwdx_const_c0.v, x_2.v, cubic_1_dwdx_const_c1.v);
  dw_dx2.v = vec_fma(cubic_2_dwdx_const_c0.v, x.v, cubic_2_dwdx_const_c1.v);
  dw_dx2_2.v = vec_fma(cubic_2_dwdx_const_c0.v, x_2.v, cubic_2_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_mul(dw_dx->v, x.v);       /* cubic_1_dwdx_const_c2 is zero. */
  dw_dx_2->v = vec_mul(dw_dx_2->v, x_2.v); /* cubic_1_dwdx_const_c2 is zero. */
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, cubic_2_dwdx_const_c2.v);
  dw_dx2_2.v = vec_fma(dw_dx2_2.v, x_2.v, cubic_2_dwdx_const_c2.v);

  /* Mask out unneeded values. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  dw_dx->v = vec_blend(mask_reg, dw_dx->v, dw_dx2.v);
  dw_dx_2->v = vec_blend(mask_reg_v2, dw_dx_2->v, dw_dx2_2.v);

#else
#error \
    "Vectorisation not supported for this kernel!!! Choose a different one or configure with --disable-hand-vec."
#endif

  /* Mask out result for particles that lie outside of the kernel function. */
  mask_t mask, mask_2;
  vec_create_mask(mask, vec_cmp_lt(x.v, vec_set1(1.f)));     /* x < 1 */
  vec_create_mask(mask_2, vec_cmp_lt(x_2.v, vec_set1(1.f))); /* x < 1 */

  dw_dx->v = vec_and_mask(dw_dx->v, mask);
  dw_dx_2->v = vec_and_mask(dw_dx_2->v, mask_2);

  /* Return everything */
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
  dw_dx_2->v = vec_mul(
      dw_dx_2->v,
      vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_plus_one_vec.v));
}

#endif /* WITH_VECTORIZATION */

/* Some cross-check functions */
void hydro_kernel_dump(int N);

#endif  // SWIFT_KERNEL_HYDRO_H
