#
# SYNOPSIS
#
#   AX_ASM_ARM_CNTVCT
#
# DESCRIPTION
#
#   Check whether the CNTVCT_EL0 exists on this platform. Defines
#   HAVE_ARMV8_CNTVCT_EL0 if true.
#
# LICENSE
#
#   Copyright (c) 2019 Matthieu Schaller <schaller@strw.leidenuniv.nl>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_ASM_ARM_CNTVCT],
[AC_CACHE_CHECK([for CNTVCT_EL0 asm instruction on ARM v8.1a],
   [ax_cv_asm_arm_cntvct_works],
    [AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>

int
main()
{
   uint64_t cc = 0;
   __asm__ __volatile__("mrs %0,  CNTVCT_EL0" : "=r"(cc));
   return 0;
}
    ]])],
    [ax_cv_asm_arm_cntvct_works=yes],
    [ax_cv_asm_arm_cntvct_works=no],
    [ax_cv_asm_arm_cntvct_works=no])])
if test "$ax_cv_asm_arm_cntvct_works" = "yes" ; then
  AC_DEFINE([HAVE_ARMV8_CNTVCT_EL0], [1],
    [Define to 1 if the ARM v8.1a instruction CNTVCT_EL0 exists.])
fi
])
