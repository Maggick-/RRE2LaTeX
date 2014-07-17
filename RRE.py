# RRE a function to put a matrix into [R]educe [R]ow [E]chelon form and write steps out LaTeX
#
# ABOUT:
#
# PUBLIC REPOSITORY:
#   - https://github.com/Maggick-/RRE2LaTeX
#
# HISTORY:
#   2014-07-17: Nick Maggio (nickmaggio [at] ymail [dot] com)
#   - Intial version
#
# TODO:
#   - Fix error handling
#   - Tidy up code
#   - Allow fractions to be printed out correctly
#   - Add function comments
#
# LICENSE:
#   The MIT License (MIT)
#
#   Copyright 2014 Nick Maggio
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy
#   of this software and associated documentation files (the "Software"), to 
#   deal in the Software without restriction, including without limitation the 
#   rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
#   sell copies of the Software, and to permit persons to whom the Software is 
#   furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in 
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
#   IN THE SOFTWARE.
#
# REFERENCES:
#   [*] https://github.com/pavdpr/RRE
#
#

__author__ = "Nick Maggio"
__copyright__ = "Copyright 2014, Maggick"
__credits__ = ["Nick Maggio"]

__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Nick Maggio"
__email__ = "nickmaggio@ymail.com"
__status__ = "Production"

import fractions
import decimal
import numpy as np

def RRE2LaTeX( A, b = [], fid = 1, sas = False, sen = False ):
    """
    Inputs:
        A   - an n x m matrix
        b   - augmented n x p matrix
            - defaults is []
        fid - file id to open
        sas - [show all step]
        sen - [supress equation numbers]

    """
    if isinstance(A,list):
        A = np.array(A)
    # fi

    if isinstance(b,list):
        b = np.array(b)
    # fi

    [ n, m ] = A.shape
    if b.size != 0:
        N = b.shape[0]
        if N != n: # TODO
            raise ValueError, 'message'
        # fi
    # fi

    i = 0

    fprintf( fid, '\\allowdisplaybreaks' )
    if not sen:
        fprintf( fid, '\\begin{flalign}' )
    else:
        fprintf( fid, '\\begin{flalign*}' )
    # fi
    fprintf( fid, RRELaTeXmat( A, b ) )
    fprintf( fid, '\t&&\\mbox{Original Matrix}\\\\' )

    while i < min( n, m ):
        if not RREhasPivot( A, i ):
            i += 1
        else:
            [ A, b ] = RREreduce( A, i, b, fid, sas )
            i += 1
        # fi
    if not sen:
        fprintf( fid, '\\end{flalign}\n' )
    else:
        fprintf( fid, '\\end{flalign*}\n' )
    # fi
# fed RRE

def RREreduce( A, i, b = [], fid = 1, sas = False):
    if isinstance(b,list):
        b = np.array(b)
    # fi
    n = A.shape[0]
    idx = RREgetPivot( A, i )
    if idx != i:
        [ A, b ] = RREswap( A, b, i, idx)
        fprintf( fid, RRELaTeXmat( A, b ) )
        var = '\t&& R_{%d' %(i+1) +'}\\leftrightarrow R_{%d' %(idx+1) +'}\\\\'
        fprintf( fid, var )
    # fi
    if A[ i, i ] == 0.0:
        raise ValueError, 'message'
    # fi
    [ A, b, s ] = RREone( A, b, i )
    if A[ i, i ] != 1.0:
        fprintf( fid, RRELaTeXmat( A, b ) )
        var = '\t&& R_{%d' %i + '}=%s' %writeFrac( s ) + '\\cdot R_{%d' %( i + 1 ) + '}\\\\' 
        fprintf( fid, var )
    # fi
    if sas:
        msg = '\\begin{array}{c}'
    # fi
    for j in range(n):
        if ( ( j != i ) and ( A[ j, i] != 0 ) ):
            [ A, b, s ] = RREzero( A, b, i, j)
            if not sas:
                fprintf( fid, RRELaTeXmat( A, b ) )
                var = '\t&& R_{%d' %( j + 1 ) + '}=R_{%d' %( j + 1 ) + '}+%s' %writeFrac( s ) + '\\cdot R_{%d' %( i + 1 )  +'}\\\\'
                fprintf( fid, var)
            else:
                msg += 'R_{' + str( j + 1 ) + '}=R_{' + str( j + 1 ) + \
                '}+' + writeFrac( s ) + '\\cdot R_{' + str( i + 1 ) + '}\\\\' 
            # fi ton
        # fi
    # rof
    if sas:
        msg += '\\end{array}\\\\'
        fprintf( fid, RRELaTeXmat( A, b ) )
        fprintf( fid, '\t&&' )
        fprintf( fid, msg )
    return A, b
# fed RREreduce

def RREhasPivot( A, i ):
    if sum( abs( A[ i: , i] ) ) == 0.0:
        return False
    else:
        return True
# fed RREhasPivot

def RREone( A, b, i):
    s = A[ i, i ]
    if s == 0: # TODO
        raise ValueError, 'message'
    # fi
    s = 1 / s
    A[ i, : ] = A[ i, : ] * s
    if b.size > 0:
        b[ i, : ] = b[ i, :] * s
    # fi
    return A, b, s
# fed RREone

def RREzero( A, b, i, j):
    if A[ i, i ] == 0: # TODO
        raise ValueError, 'message'
    elif i == j: # TODO
        raise ValueError, 'message'
    s = A[ j, i ] / A[ i, i]
    A[ j, :] = A[ j, :] - s * A[ i, : ]
    if b.size > 0:
        b[ j, : ] = b[ j, : ] - s * b[ i, :] 
    return A, b, s
# fed RREzero

def RREgetPivot( A, i ):
    return int(np.nonzero( abs( A[ i:, i ] ) == max( abs( A[ i:, i ] ) ) )[0]+i)

# fed RREgetPivot

def RREswap( A, b = [], i = 0, j = 0):
    if isinstance(b,list):
        b = np.array(b)
    # fi
    if i == j:
        return A, b
    # fi
    [ n, m ] = A.shape
    A[ [ j, i ], : ] = A[ [ i, j ], : ]
    if b.size > 0:
        b[ [ j, i ], : ] = b[ [ i, j ], : ]
    # fi
    return A, b
# fed RREswap

def RRELaTeXmat( A, b = [] ):
    """

    """
    if isinstance(b,list):
        b = np.array(b)
    # fi
    if b.size == 0:
        p = 0
    else:
        p = b.shape[1]
        A = np.concatenate( ( A , b ), 1 )
    # fi
    [ n, m ] = A.shape
    LaTeX = '\t\\left[\\begin{array}{'
    for i in range(m-p):
        LaTeX += 'l'
    # fi
    if p > 0:
        LaTeX += '|'
        for i in range(p):
            LaTeX += 'l'
        # rof
    # fi
    LaTeX += '}\n'
    for i in range(n):
        LaTeX += '\t\t'
        for j in range(m):
            LaTeX += writeFrac( A[ i, j ] )
            if j != m - 1:
                LaTeX += ' & '
            # fi
        # rof
        LaTeX += '\\\\\n'
    # rof
    LaTeX += '\t\\end{array}\\right]'
    return LaTeX
# fed RRELaTeXmat

def writeFrac(dec):
    var = dec2frac(dec)
    if var.denominator == 1:
        return str( var )
    else:
        return '%.4f' %dec
        #return frac2LaTeX( var )
    # fi
# fed writeFrac


def dec2frac(dec):
    return fractions.Fraction.from_decimal(decimal.Decimal(str(dec))).limit_denominator(100000000)
# fed dec2frac

def frac2LaTeX(frac):
    if isinstance(frac,fractions.Fraction):
        return '\\frac{' + str(frac.numerator) + '}{' + str(frac.denominator) + '}'
    # fi
# fed frac2LaTeX

def fprintf(fid,message):
    if isinstance(fid,file):
        fid.write(message)
    else:
        print message
    # fi
# fed fprintf

if __name__ == "__main__":

    S = np.array([[1.,3,1],[1,1,-1],[3,11,5]])
    v = np.array([[9.],[1],[35]])

    print RRE2LaTeX(S,v)

