#(C) Christian Benner 2020.

import numpy as np

def areSNPsIncluded( snps, n_snps, fname ):
    for snp in snps:
        if snp >= n_snps:
            print('Cannot read dosages for ' + str( snp ) + 'th SNP. File "' + fname + '" contains only ' + str( n_snps ) + ' SNPs!')
            return

def convertIntToFloat( x, n_bytes ):
    if type( x ) is np.ndarray:
        return convertIntToFloat_array( x, n_bytes )
    else:
        return convertIntToFloat_scalar( x, n_bytes )

def convertIntToFloat_array( x, n_bytes ):
    int_na = getIntNA( n_bytes )
    y = np.zeros( len( x ) )
    y[ x == int_na ] = np.nan
    y[ x != int_na ] = np.ldexp( x[ x != int_na ], -1 * ( 8 * n_bytes - 2 ) )
    
    return y

def convertIntToFloat_scalar( x, n_bytes ):
    if x == getIntNA( n_bytes ):
        return np.nan
    else:
        return np.ldexp( x, -1 * ( 8 * n_bytes - 2 ) )

def getIntNA( n_bytes ):
    if n_bytes == 2:
        return 53248
    elif n_bytes == 4:
        return 3489660928
    elif n_bytes == 8:
        return 14987979559889010688
    elif n_bytes == 1:
        return 208
    else:
        print('Only 1, 2, 4 and 8 bytes are supported!')
        return
