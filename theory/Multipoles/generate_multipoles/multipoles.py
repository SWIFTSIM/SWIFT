import numpy as np
import sys


def factorial(x):
    if x == 0:
        return 1
    else:
        return x * factorial(x - 1)


SUFFIXES = {1: "st", 2: "nd", 3: "rd"}


def ordinal(num):
    suffix = SUFFIXES.get(num % 10, "th")
    return str(num) + suffix


# Get the order
order = int(sys.argv[1])

print("-------------------------------------------------")
print("Generating code for multipoles of order", order, "(only).")
print("-------------------------------------------------\n")

print("-------------------------------------------------")
print("Multipole structure:")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n" % (order - 1))

print("/* %s order terms */" % ordinal(order))

# Create all the terms relevent for this order
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                print("float M_%d%d%d;" % (i, j, k))

if order > 0:
    print("#endif")

print("")
print("-------------------------------------------------")

print("Field tensor structure:")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n" % (order - 1))

print("/* %s order terms */" % ordinal(order))

# Create all the terms relevent for this order
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                print("float F_%d%d%d;" % (i, j, k))
if order > 0:
    print("#endif")

print("")
print("-------------------------------------------------")

print("gravity_field_tensors_add():")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d" % (order - 1))

print("/* %s order terms */" % ordinal(order))

# Create all the terms relevent for this order
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                print("la->F_%d%d%d += lb->F_%d%d%d;" % (i, j, k, i, j, k))
if order > 0:
    print("#endif")

print("")
print("-------------------------------------------------")

print("gravity_multipole_add():")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d" % (order - 1))

print("/* %s order terms */" % ordinal(order))

# Create all the terms relevent for this order
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                print("ma->M_%d%d%d += mb->M_%d%d%d;" % (i, j, k, i, j, k))

if order > 0:
    print("#endif")

print("")
print("-------------------------------------------------")

print("gravity_multipole_equal():")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d" % (order - 1))

# Create all the terms relevent for this order
print("/* Manhattan Norm of %s order terms */" % ordinal(order))
print("const float order%d_norm = " % order, end=" ")
first = True
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                if first:
                    first = False
                else:
                    print("+", end=" ")
                print("fabsf(ma->M_%d%d%d)" % (i, j, k), end=" ")
                print("+ fabsf(mb->M_%d%d%d)" % (i, j, k), end=" ")
print(";\n")
print("/* Compare %s order terms above 1%% of norm */" % ordinal(order))
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                print(
                    "if (fabsf(ma->M_%d%d%d + mb->M_%d%d%d) > 0.01f * order%d_norm &&"
                    % (i, j, k, i, j, k, order)
                )
                print(
                    "    fabsf(ma->M_%d%d%d - mb->M_%d%d%d) / fabsf(ma->M_%d%d%d + mb->M_%d%d%d) > tolerance) {"
                    % (i, j, k, i, j, k, i, j, k, i, j, k)
                )
                print('  message("M_%d%d%d term different");' % (i, j, k))
                print("  return 0;")
                print("}")

if order > 0:
    print("#endif")


print("")
print("-------------------------------------------------")

print("gravity_multipole_compute_power():")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d" % (order - 1))

print("/* %s order terms */" % ordinal(order))

# Add the terms to the multipole power
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                fact1 = factorial(i) * factorial(j) * factorial(k)
                fact2 = float(factorial(i + j + k))
                frac = fact1 / fact2
                if frac == 1.0:
                    print(
                        "power[%d] += m->M_%d%d%d * m->M_%d%d%d;"
                        % (order, i, j, k, i, j, k)
                    )
                else:
                    print(
                        "power[%d] += %12.15e * m->M_%d%d%d * m->M_%d%d%d;"
                        % (order, frac, i, j, k, i, j, k)
                    )

print("")
print("m->power[%d] = sqrt(power[%d]);" % (order, order))

if order > 0:
    print("#endif")

print("")
print("-------------------------------------------------")

print("gravity_P2M(): (loop)")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d" % (order - 1))

print("/* %s order terms */" % ordinal(order))

# Create all the terms relevent for this order
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                if order % 2 == 0:
                    print("M_%d%d%d += m * X_%d%d%d(dx);" % (i, j, k, i, j, k))
                else:
                    print("M_%d%d%d += -m * X_%d%d%d(dx);" % (i, j, k, i, j, k))

if order > 0:
    print("#endif")

print("")
print("-------------------------------------------------")

print("gravity_P2M(): (storing)")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d" % (order - 1))

print("/* %s order terms */" % ordinal(order))

# Create all the terms relevent for this order
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                print("m->m_pole.M_%d%d%d = M_%d%d%d;" % (i, j, k, i, j, k))

if order > 0:
    print("#endif")


print("")
print("-------------------------------------------------")

print("gravity_M2M():")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d" % (order - 1))

print("/* Shift %s order terms (1st order mpole (all 0) commented out) */" % ordinal(order))

# Create all the terms relevent for this order
for i in range(order + 1):
    for j in range(order + 1):
        for k in range(order + 1):
            if i + j + k == order:
                print("m_a->M_%d%d%d = m_b->M_%d%d%d" % (i, j, k, i, j, k), end=" ")

                for ii in range(order + 1):
                    for jj in range(order + 1):
                        for kk in range(order + 1):

                            if not (ii == 0 and jj == 0 and kk == 0):
                                for iii in range(order + 1):
                                    for jjj in range(order + 1):
                                        for kkk in range(order + 1):
                                            if (
                                                ii + iii == i
                                                and jj + jjj == j
                                                and kk + kkk == k
                                            ):
                                                if iii + jjj + kkk == 1:
                                                    print(
                                                        "/* + X_%d%d%d(dx) * m_b->M_%d%d%d */"
                                                        % (ii, jj, kk, iii, jjj, kkk),
                                                        end=" ",
                                                    )
                                                else:
                                                    print(
                                                        "+ X_%d%d%d(dx) * m_b->M_%d%d%d"
                                                        % (ii, jj, kk, iii, jjj, kkk),
                                                        end=" ",
                                                    )

                print(";")

if order > 0:
    print("#endif")


print("")
print("-------------------------------------------------")

print("gravity_M2L():")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n" % (order - 1))

# Loop over LHS order
for l in range(order + 1):
    print(
        "/* Compute %s order field tensor terms (addition to rank %d) */"
        % (ordinal(order), l)
    )

    for i in range(l + 1):
        for j in range(l + 1):
            for k in range(l + 1):
                if i + j + k == l:
                    print("l_b->F_%d%d%d +=" % (i, j, k), end=" ")

                    first = True
                    for ii in range(order + 1):
                        for jj in range(order + 1):
                            for kk in range(order + 1):
                                if ii + jj + kk == order - l:
                                    if first:
                                        first = False
                                    else:
                                        print("+", end=" ")
                                    print(
                                        "m_a->M_%d%d%d * D_%d%d%d(dx, dy, dz, r_inv)"
                                        % (ii, jj, kk, i + ii, j + jj, k + kk),
                                        end=" ",
                                    )
                    print(";")
    print("")

if order > 0:
    print("#endif")


print("")
print("-------------------------------------------------")

print("gravity_L2L():")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n" % (order - 1))

# Loop over LHS order
for l in range(order + 1):
    print(
        "/* Shift %s order field tensor terms (addition to rank %d) */"
        % (ordinal(order), l)
    )

    for i in range(l + 1):
        for j in range(l + 1):
            for k in range(l + 1):
                if i + j + k == l:
                    print("la->F_%d%d%d +=" % (i, j, k), end=" ")

                    first = True
                    for ii in range(order + 1):
                        for jj in range(order + 1):
                            for kk in range(order + 1):
                                if ii + jj + kk == order - l:
                                    if first:
                                        first = False
                                    else:
                                        print("+", end=" ")
                                    print(
                                        "X_%d%d%d(dx) * lb->F_%d%d%d"
                                        % (ii, jj, kk, i + ii, j + jj, k + kk),
                                        end=" ",
                                    )
                    print(";")
    print("")

if order > 0:
    print("#endif")

print("")
print("-------------------------------------------------")

print("gravity_L2P():")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n" % (order - 1))

    print("/* %s order contributions */" % (ordinal(order - 1)))

    for r in range(3):
        print("gp->a_grav[%d] +=" % (r), end=" ")

        first = True
        for i in range(order + 1):
            for j in range(order + 1):
                for k in range(order + 1):
                    if i + j + k == order - 1:
                        if first:
                            first = False
                        else:
                            print("+", end=" ")
                        if r == 0:
                            ii = i + 1
                            jj = j
                            kk = k
                        if r == 1:
                            ii = i
                            jj = j + 1
                            kk = k
                        if r == 2:
                            ii = i
                            jj = j
                            kk = k + 1
                        print(
                            "X_%d%d%d(dx) * lb->F_%d%d%d" % (i, j, k, ii, jj, kk),
                            end=" ",
                        )
        print(";")

    print("")
if order > 0:
    print("#endif")

print("")
print("-------------------------------------------------")

print("gravity_M2P():")
print("-------------------------------------------------\n")

if order > 0:
    print("#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n" % (order - 1))

print("/* %s order contributions */" % (ordinal(order)))


for r in range(4):
    if r == 0:
        print("*f_x =", end=" ")
    if r == 1:
        print("*f_y =", end=" ")
    if r == 2:
        print("*f_z =", end=" ")
    if r == 3:
        print("*pot =", end=" ")

    first = True
    for i in range(order + 1):
        for j in range(order + 1):
            for k in range(order + 1):
                if i + j + k == order:
                    if first:
                        first = False
                    else:
                        print("+", end=" ")
                    if r == 0:
                        ii = i + 1
                        jj = j
                        kk = k
                    if r == 1:
                        ii = i
                        jj = j + 1
                        kk = k
                    if r == 2:
                        ii = i
                        jj = j
                        kk = k + 1
                    if r == 3:
                        ii = i
                        jj = j
                        kk = k
                    print("m->M_%d%d%d * d.D_%d%d%d" % (i, j, k, ii, jj, kk), end=" ")

    print(";")

print("")

if order > 0:
    print("#endif")

print("")
print("-------------------------------------------------")
