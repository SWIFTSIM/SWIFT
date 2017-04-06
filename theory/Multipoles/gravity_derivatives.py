import numpy as np
import sys

SUFFIXES = {1: 'st', 2: 'nd', 3: 'rd'}
def ordinal(num):
    suffix = SUFFIXES.get(num % 10, 'th')
    return str(num) + suffix

# Get the order
n = int(sys.argv[1])
verbose = int(sys.argv[2])

# Get the sorting indexes
def argsort(seq):
    return sorted(range(len(seq)), key=seq.__getitem__)

def pochhammer(x):
    if x == 1:
        return "- 1."
    if x == 2:
        return "+ 3."
    if x == 3:
        return "- 15."
    if x == 4:
        return "+ 105."
    if x == 5:
        return "- 945."
    if x == 6:
        return "+ 10395."
    if x == 7:
        return "- 135135."
    if x == 8:
        return "+ 2027025."
    if x == 9:
        return "- 34459425."
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

# All non-zero derivatives of u(r_x,r_y,r_z) = r_x^2 + r_y^2 + r_z^2 divided by two
def u_derivative(deriv):
    if len(deriv) > 2:
        return ""
    elif deriv == "x":
        return "r_x"
    elif deriv == "y":
        return "r_y"
    elif deriv == "z":
        return "r_z"
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
norm_pi = []
norm_b = []
for i in range(num_terms):
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
            cardinality += 1
            for j in range(n):
                if permutations[i][j] == current:
                    my_str += str(range(n)[j]+1)
            my_str += "]"
    partitions.append(my_str)
    norm_pi.append(cardinality)
    norm_b.extend([this_norm_b])

# Sort partitions in order of |pi|
index = argsort(norm_pi)
norm_pi = [ norm_pi[i] for i in index]
norm_b = [ norm_b[i] for i in index]
partitions = [ partitions[i] for i in index]
permutations = [ permutations[i] for i in index]

if verbose:
    print "Paritions of {1...%d}"%n
    print "Number of partitions:", num_terms
    for i in range(num_terms):
        print "%3d)"%(i+1), permutations[i], "  ---> pi =",
        print partitions[i],
        for j in range(n-norm_pi[i]):
            print " ",
        print "  ---> |pi| =", norm_pi[i],
        print "  ---> [|B|] =", norm_b[i]

# Now we can start generating functions
print "----------------------------------------------------------"
print "Generating code for gravity derivatives of order", n, "(only)."
print "----------------------------------------------------------\n"

print "/*********************************/"
print "/* %s order gravity derivatives */"%ordinal(n)
print "/*********************************/\n"

# Write out the functions
for i in range(n+1):
    for j in range(n+1):
        for k in range(n+1):
            if i + j + k == n:

                #if i != 4 or j != 2 or k != 0:
                #    continue
                
                derivatives = map_derivatives(i, j, k)

                terms = []
                count_terms = []
                norm_pi2 = []
                count = 0
                
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
                                count += 1
                                break
                            
                    if non_zero:
                        
                        # Sort the blocks in lexicographic order
                        part2 = []
                        temp = ""
                        for c in range(len(part)):
                            if part[c] == '[':
                                temp = ""
                            elif part[c] == ']':
                                part2.append(temp)
                            else:
                                temp += part[c]
                        part2 = sorted(part2)

                        terms.append(part2)
                        norm_pi2.append(norm_pi[p])

                # Sort list of blocks
                index = argsort(terms)
                terms = [ terms[u] for u in index]
                norm_pi2 = [ norm_pi2[u] for u in index]

                # Make the list unique
                terms2 = []
                norm_pi3 = []
                for b in range(len(terms)):
                    if terms2 == [] or terms2[-1] != terms[b]:
                        terms2.append(terms[b])
                        count_terms.append(1)
                        norm_pi3.append(norm_pi2[b])
                    else:
                        count_terms[-1] += 1
                    
                            
                # Write the documentation and prototype
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


                # Print the whole thing
                for b in range(len(terms2)):

                    part2 = terms2[b]
                    #print "/* %18s*/" %(part2),
                    
                    # Write the derivative of phi
                    print derivative_phi(norm_pi3[b]),

                    # Write out multiplicity
                    if count_terms[b] > 1:
                        print "* %.1f "%count_terms[b],
                        
                    # Write out the derivatives of u
                    first = True
                    for c in part2:
                        deriv = u_derivative(c)
                        if deriv != "":
                            if first :
                                print " * (",
                                first = False
                            else:
                                print "*",
                            print deriv,
                    if not first:
                        print ")",
                    
                # Finish the code fo the function
                print ";"
                print "/* %d zero-valued terms not written out */"%(count)
                print "}\n"
