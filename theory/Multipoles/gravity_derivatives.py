import numpy as np
import sys

SUFFIXES = {1: 'st', 2: 'nd', 3: 'rd'}
def ordinal(num):
    suffix = SUFFIXES.get(num % 10, 'th')
    return str(num) + suffix

# Get the order
n = int(sys.argv[1])

# All non-zero derivatives of u(r_x,r_y,r_z) = r_x^2 + r_y^2 + r_z^2
U_100 = "2. * r_x"
U_200 = "2."
U_010 = "2. * r_y"
U_020 = "2."
U_001 = "2. * r_z"
U_002 = "2."

# Generate all partitions of {1....n}
partitions = []
temp = []
for i in range(n):
    temp.append(1)
partitions.extend([list(temp)])

if n > 1:
    while temp != range(1,n+1):
        i = n-1
        while i >= 0 and temp[0] == 1:
            if i > 0:
                max_i = max(temp[0:i])
            else:
                max_i = 0
            if temp[i] != n and not max_i+1 < temp[i]+1:
                #print "ok"
                temp[i] += 1
                for j in range(i+1, n):
                    temp[j] = 1
                i = n-1
                if temp[0] == 1:
                    partitions.extend([list(temp)])
            else:
                i -= 1

                
# Print the partitions
print "Number of partitions:", len(partitions)
for i in range(len(partitions)):
    print "%3d)"%(i+1), partitions[i], "-->",
    for current in range(1,n+1):
        count = 0
        for j in range(n):
            if partitions[i][j] == current:
                count += 1
        if count > 0:
            print "[",
            for j in range(n):
                if partitions[i][j] == current:
                    print range(n)[j]+1,
            print "]",
    print ""
        
