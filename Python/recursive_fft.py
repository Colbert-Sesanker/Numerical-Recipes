# Colbert Sesanker, 2012, fft used in Project Connjur-SE
# minimal dependencies for GAE cloud use

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
    if (complex_fft(x) == 0):
        return "input data size not a power of two"       
    return map(abs, complex_fft(x))

def fft_test():
    return fft([math.sin(x) for x in range(16)])

def nextPow2(num):     
    for l in map(lambda x: num - x, [2**i for i in range(35)]):
        if (l <= 0):    
            return num - l  

def parse(sched):
    s = [[0]*nextPow2(max([x[0] for x in sched])) for l in range(nextPow2(max([x[1] for x in sched])))] 
    for i in sched:        
        s[i[1]-1] [i[0]-1] = 1   # point [x,y] = s[y][x]
    return s   

def fft2(x):
    return map(fft, map(list, zip(*map(fft, parse(x))) ) )




