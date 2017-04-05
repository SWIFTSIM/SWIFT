import numpy as np
import sys

SUFFIXES = {1: 'st', 2: 'nd', 3: 'rd'}
def ordinal(num):
    suffix = SUFFIXES.get(num % 10, 'th')
    return str(num) + suffix

# Get the order
n = int(sys.argv[1])
#verbose = bool(sys.argv[2])

def pochhammer(x):
    if x == 1:
        return "-(1. / 2.)"
    if x == 2:
        return "(3. / 4.)"
    if x == 3:
        return "-(15. / 8.)"
    if x == 4:
        return "(105. /16.)"
    if x == 5:
        return "-(945. /32.)"
    if x == 6:
        return "(10395. / 64.)"
    if x == 7:
        return "-(135135. / 128.)"
    if x == 8:
        return "(2027025. / 256.)"
    if x == 9:
        return "-(34459425. / 512.)"
    else:
        print "Invalid x"
        exit(-1)

# Returns nth derivative of 1/sqrt[u]
def derivative_phi(n):
    power = 2 * n + 1
    string = str(pochhammer(n))
    #string += " * r_inv^%d"%power
    string += " * r_inv"
    for i in range(power - 1):
        string += " * r_inv"
    return string

def u_derivative_is_non_zero(deriv):
    if len(deriv) > 2:
        return False
    elif len(deriv) == 2:
        if deriv == "xx" or deriv == "yy" or deriv == "zz":
            return True
        else:
            return False
    elif len(deriv) == 1:
        return True

# All non-zero derivatives of u(r_x,r_y,r_z) = r_x^2 + r_y^2 + r_z^2
def u_derivative(deriv):
    if len(deriv) > 2:
        return ""
    elif deriv == "x":
        return "2. * r_x"
    elif deriv == "y":
        return "2. * r_y"
    elif deriv == "z":
        return "2. * r_z"
    elif deriv == "xx":
        return "2."
    elif deriv == "yy":
        return "2."
    elif deriv == "zz":
        return "2."
    else:
        return ""
    
# Return the array of derivatives to compute from a multi-index
def map_derivatives(num_x, num_y, num_z):

    array = []
    for i in range(num_x):
        array.append('x')
    for i in range(num_y):
        array.append('y')
    for i in range(num_z):
        array.append('z')
    
    return array

    

    
# Generate all permutations of {1....n}
permutations = []
temp = []
for i in range(n):
    temp.append(1)
permutations.extend([list(temp)])

# Use permutations in lexicographic order to generate partitions
while temp != range(1,n+1):
    i = n-1
    while i >= 0 and temp[0] == 1:
        # Maximal value up to that point
        if i > 0:
            max_i = max(temp[0:i])
        else:
            max_i = 0
        # Can't go beyond the order not above the max+1
        if temp[i] != n and not max_i+1 < temp[i]+1:
            temp[i] += 1
            for j in range(i+1, n):
                temp[j] = 1
            i = n-1
            if temp[0] == 1:
                permutations.extend([list(temp)])
        else:
            i -= 1

                
# Print the partitions and compute cardinality of each partition and block
num_terms = len(permutations)
partitions = []
print "Paritions of {1...%d}"%n
print "Number of partitions:", num_terms
norm_pi = []
norm_b = []
for i in range(num_terms):
    print "%3d)"%(i+1), permutations[i], "  ---> pi =",
    cardinality = 0
    this_norm_b = []
    my_str = ""
    for current in range(1,n+1):
        count = 0
        for j in range(n):
            if permutations[i][j] == current:
                count += 1
        if count > 0:
            this_norm_b.append(count)
            my_str += "["
            print "[",
            cardinality += 1
            for j in range(n):
                if permutations[i][j] == current:
                    print range(n)[j]+1,
                    my_str += str(range(n)[j]+1)
            print "]",
            my_str += "]"
    partitions.append(my_str)
    for j in range(n-cardinality):
        print "   ",
    print "  ---> |pi| =", cardinality,
    print "  ---> [|B|] =", this_norm_b 
    norm_pi.append(cardinality)
    norm_b.extend([this_norm_b])

# Now we can start generating functions
print "----------------------------------------------------------"
print "Generating code for gravity derivatives of order", n, "(only)."
print "----------------------------------------------------------\n"

print "/*********************************/"
print "/* %s order gravity derivatives */"%ordinal(n)
print "/*********************************/\n"

# Create all the terms relevent for this order
for i in range(n+1):
    for j in range(n+1):
        for k in range(n+1):
            if i + j + k == n:

                # Wrtie the documentation and prototype
                print "/**"
                print " * @brief Compute \\f$ \\frac{\\partial^%d}{"%n,
                if i > 0: print "\\partial_x^%d"%i,
                if j > 0: print "\\partial_y^%d"%j,
                if k > 0: print "\\partial_z^%d"%k,
                print "}\\phi(x, y, z} \\f$."
                print " *"
                print " * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)"
                print " */"
                print "__attribute__((always_inline)) INLINE static double D_%d%d%d(double r_x, double r_y, double r_z, double r_inv) {"%(i,j,k)
                print "return",

                derivatives = map_derivatives(i, j, k)
                
                first_term = True
                for p in range(num_terms):
                    # Get the partition
                    part = ""
                    part += partitions[p]
                    
                    # Replace the numbers by derivatives
                    for b in range(len(derivatives)):
                        part = part.replace(str(b+1), derivatives[b])
                       
                    # Check whether the derivatives are non-zero
                    non_zero = True
                    for b in part:
                        if b=='[':
                            deriv = ""
                        elif b.isalnum():
                            deriv += b
                        elif b==']':
                            if not u_derivative_is_non_zero(deriv):
                                non_zero = False
                                break
                    # Print the derivative of phi(u)
                    if non_zero:

                        #print "/* %18s */"%part,
                        
                        # Get the first symbol
                        if first_term:
                            first_term = False
                        else:
                            print "+",

                        #Write the derivative of phi
                        #print "\\phi^(%d)(u) * ("%norm_pi[p],
                        print derivative_phi(norm_pi[p]), "* (",
                            
                        # Write out the derivatives of u
                        first = True
                        for b in part:
                            if b=='[':
                                deriv = ""
                                if first:
                                    first = False
                                else:
                                    print "*",
                            elif b.isalnum():
                                deriv += b
                            elif b==']':
                                print u_derivative(deriv),
                    
                        print ")",
                    
                # Finish the code fo the function
                print ";"
                print "}\n"
                #exit(-1)
