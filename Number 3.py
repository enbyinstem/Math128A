import math

def maxPrimeFactors (x):
    maxPrime = -1;

    while (x % 2) == 0:
        maxPrime = 2
        x >>= 1

    
    for i in range(3, int(math.sqrt(x))+1, 2):
        while (x % i) == 0:
            maxPrime = i
            x = x/i


    if x > 2:
        maxPrime = x

    return int(maxPrime)

x=600851475143
print(maxPrimeFactors(x))
