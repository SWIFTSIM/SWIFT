import numpy as np
import sys

def factorial(x):
    if x == 0:
        return 1
    else:
        return x * factorial(x-1)

SUFFIXES = {1: 'st', 2: 'nd', 3: 'rd'}
def ordinal(num):
    suffix = SUFFIXES.get(num % 10, 'th')
    return str(num) + suffix

# Get the order
order = int(sys.argv[1])

print "-------------------------------------------------"
print "Generating code for multipoles of order", order, "(only)."
print "-------------------------------------------------\n"

print "-------------------------------------------------"
print "Multipole structure:"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

print "/* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "float M_%d%d%d;"%(i,j,k)

if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "Field tensor structure:"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

print "/* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "float F_%d%d%d;"%(i,j,k)
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "gravity_field_tensors_add():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "/* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "la->F_%d%d%d += lb->F_%d%d%d;"%(i,j,k,i,j,k)
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "gravity_multipole_add():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "/* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "ma->M_%d%d%d += mb->M_%d%d%d;"%(i,j,k,i,j,k)
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "gravity_multipole_equal():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

# Create all the terms relevent for this order
print "/* Manhattan Norm of %s order terms */"%ordinal(order)
print "const float order%d_norm = "%order,
first = True
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                if first:
                    first = False
                else:
                    print "+",
                print "fabsf(ma->M_%d%d%d)"%(i,j,k),
                print "+ fabsf(mb->M_%d%d%d)"%(i,j,k),
print ";\n"
print "/* Compare %s order terms above 1%% of norm */"%ordinal(order)
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "if (fabsf(ma->M_%d%d%d + mb->M_%d%d%d) > 0.01f * order%d_norm &&"%(i,j,k,i,j,k,order)
                print "    fabsf(ma->M_%d%d%d - mb->M_%d%d%d) / fabsf(ma->M_%d%d%d + mb->M_%d%d%d) > tolerance) {"%(i,j,k,i,j,k,i,j,k,i,j,k)
                print "  message(\"M_%d%d%d term different\");"%(i,j,k)
                print "  return 0;"
                print "}"

if order > 0:
    print "#endif"


print ""
print "-------------------------------------------------"

print "gravity_P2M(): (loop)"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "/* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                if order % 2 == 0:
                    print "M_%d%d%d += m * X_%d%d%d(dx);"%(i,j,k,i,j,k)
                else:
                    print "M_%d%d%d += -m * X_%d%d%d(dx);"%(i,j,k,i,j,k)
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"
    
print "gravity_P2M(): (storing)"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "/* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "m->m_pole.M_%d%d%d = M_%d%d%d;"%(i,j,k,i,j,k)
if order > 0:
    print "#endif"


print ""
print "-------------------------------------------------"

print "gravity_M2M():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "/* Shift %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "m_a->M_%d%d%d = m_b->M_%d%d%d"%(i,j,k,i,j,k),

                for ii in range(order+1):
                    for jj in range(order+1):
                        for kk in range(order+1):

                            if not(ii == 0 and jj == 0 and kk == 0):
                                for iii in range(order+1):
                                    for jjj in range(order+1):
                                        for kkk in range(order+1):
                                            if ii+iii == i and jj+jjj == j and kk+kkk == k:
                                                print "+ X_%d%d%d(dx) * m_b->M_%d%d%d"%(ii, jj, kk, iii, jjj, kkk),

                                        
                print ";"
if order > 0:
    print "#endif"

    
print ""
print "-------------------------------------------------"

print "gravity_M2L():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

# Loop over LHS order
for l in range(order + 1):
    print "/* Compute %s order field tensor terms (addition to rank %d) */"%(ordinal(order), l)

    for i in range(l+1):
        for j in range(l+1):
            for k in range(l+1):
                if i + j + k == l:
                    print "l_b->F_%d%d%d +="%(i,j,k),

                    first = True
                    for ii in range(order+1):
                        for jj in range(order+1):
                            for kk in range(order+1):
                                if ii + jj + kk  == order - l:
                                    if first:
                                        first = False
                                    else:
                                        print "+",
                                    print "m_a->M_%d%d%d * D_%d%d%d(dx, dy, dz, r_inv)"%(ii,jj,kk,i+ii,j+jj,k+kk),
                    print ";"
    print ""
    
if order > 0:
    print "#endif"


print ""
print "-------------------------------------------------"

print "gravity_L2L():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

# Loop over LHS order
for l in range(order + 1):
    print "/* Shift %s order field tensor terms (addition to rank %d) */"%(ordinal(order), l)

    for i in range(l+1):
        for j in range(l+1):
            for k in range(l+1):
                if i + j + k == l:
                    print "la->F_%d%d%d +="%(i,j,k),

                    first = True
                    for ii in range(order+1):
                        for jj in range(order+1):
                            for kk in range(order+1):
                                if ii + jj + kk  == order - l:
                                    if first:
                                        first = False
                                    else:
                                        print "+",
                                    print "X_%d%d%d(dx) * lb->F_%d%d%d"%(ii,jj,kk,i+ii,j+jj,k+kk),
                    print ";"
    print ""
    
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "gravity_L2P():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

    print "/* %s order contributions */"%(ordinal(order-1))

    for r in range(3):
        print "gp->a_grav[%d] +="%(r),

        first = True
        for i in range(order + 1):
            for j in range(order + 1):
                for k in range(order + 1):
                    if i + j + k == order-1:
                        if first:
                            first = False
                        else:
                            print "+",
                        if r == 0:
                            ii = i+1
                            jj = j
                            kk = k
                        if r == 1:
                            ii = i
                            jj = j+1
                            kk = k
                        if r == 2:
                            ii = i
                            jj = j
                            kk = k+1
                        print "X_%d%d%d(dx) * lb->F_%d%d%d"%(i,j,k,ii,jj,kk),
        print ";"

    print ""
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

