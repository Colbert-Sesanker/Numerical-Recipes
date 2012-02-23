# Colbert Sesanker, 2012, used in Project Connjur-SE

import math

def complex_fft(x): 
    """ input: vector with 2^n points, outputs fft of same length"""       
    N = len(x)   
    M = N /2 
    if N == 1 :                 # recursion base, test if length x is one
        return x
    if (N % 2):                 # using numeric output from N % 2 as boolean
        return 0
    x_even = x[::2]      
    x_odd  = x[1::2]
    X_even = fft(x_even)
    X_odd  = fft(x_odd)
    X = [0]*N                   # Allocate Zeros    
    for k in range(M) :
        twiddle = math.e ** (2j * math.pi * k / N)
        X[k] = X_even[k] + X_odd[k] * twiddle
        X[M+k] = X_even[k] - X_odd[k] * twiddle   
    return X

def fft(x):
    complex_spectrum = complex_fft(x)
    return map(abs, complex_spectrum)
    
