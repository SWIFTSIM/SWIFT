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

print "----------------------------------------------------"
print "Generating code for vector powers of order", order, "(only)."
print "----------------------------------------------------\n"

print "/***************************/"
print "/* %s order vector powers */"%ordinal(order)
print "/***************************/\n"

# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                fact = factorial(i) * factorial(j) * factorial(k)
                print "/**"
                print "* @brief Compute \\f$ \\frac{1}{(%d,%d,%d)!}\\vec{v}^{(%d,%d,%d)} \\f$."%(i,j,k,i,j,k)
                print "*"
                print "* Note \\f$ \\vec{v}^{(%d,%d,%d)} ="%(i,j,k),
                if i > 0: print "v_x^%d"%i,
                if j > 0: print "v_y^%d"%j,
                if k > 0: print "v_z^%d"%k,
                print "\\f$"
                print "* and \\f$ \\frac{1}{(%d,%d,%d)!} = 1/(%d!*%d!*%d!) = 1/%d = %e \\f$"%(i,j,k,i,j,k, fact, 1./fact)
                print "*"
                print "* @param v vector (\\f$ v \\f$)."
                print "*/"
                print "__attribute__((always_inline, const)) INLINE static double X_%d%d%d(const double v[3]) {"%(i,j,k)
                print ""
                print "  return",
                if fact != 1:
                    print "%12.15e"%(1./fact),
                else:
                    print "1.",
                for ii in range(i):
                    print "* v[0]",
                for jj in range(j):
                    print "* v[1]",
                for kk in range(k):
                    print "* v[2]",
                print ";"
                print "}\n"


