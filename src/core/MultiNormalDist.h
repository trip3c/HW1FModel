//
//  MultiNormalDist.h
//  test
//
//

#ifndef MultiNormalDist_h
#define MultiNormalDist_h


# include <math.h>
# include <iostream>
# include <string>

using namespace std;

vector<double> multinormal_sample ( int m, int n, const vector<double> a, const vector<double>  mu, int seed );
vector<double> r8po_fa ( int n, const vector<double> a );
vector<double> r8vec_normal_01_new ( int n, int seed );
vector<double> r8vec_uniform_01_new ( int n, int seed );

//****************************************************************************80



vector <double> multinormal_sample ( int m, int n, const vector<double> a, const vector<double>  mu, int seed )

{
    int i;
    int j;
    int k;
    
    //
    //  Compute the upper triangular Cholesky factor R of the variance-covariance
    //  matrix.
    //
    vector<double> r = r8po_fa ( m, a );
    vector<double> x(m*n);
    
    if ( r.empty() )
    {
        cout << "\n";
        cout << "MULTINORMAL_SAMPLE - Fatal error!\n";
        cout << "  The variance-covariance matrix is not positive definite symmetric.\n";
//        exit ( 1 );
        return x;
    }
    //
    //  Y = MxN matrix of samples of the 1D normal distribution with mean 0
    //  and variance 1.
    //
    vector<double> y = r8vec_normal_01_new ( m*n, seed );
    //
    //  Compute X = MU + R' * Y.
    //
    
    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < m; i++ )
        {
            x[i+j*m] = mu[i];
            for ( k = 0; k < m; k++ )
            {
                x[i+j*m] = x[i+j*m] + r[k+i*m] * y[k+j*m];
            }
        }
    }
    
    return x;
}

vector <double> r8po_fa ( int n, const vector<double> a )

//
//  Purpose:
//
//    R8PO_FA factors a R8PO matrix.
//
//    The R8PO storage format is appropriate for a symmetric positive definite
//    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
//    upper triangular matrix, so it will be in R8GE storage format.)
//
//    Only the diagonal and upper triangle of the square array are used.
//    This same storage format is used when the matrix is factored by
//    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
//    is set to zero.
//
//    The positive definite symmetric matrix A has a Cholesky factorization
//    of the form:
//
//      A = R' * R
//
//    where R is an upper triangular matrix with positive elements on
//    its diagonal.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix in R8PO storage.
//
//    Output, double R8PO_FA[N*N], the Cholesky factor in SGE
//    storage, or NULL if there was an error.
//
{
    int i;
    int j;
    int k;
    double s;
    
    vector<double> b(n*n);
    
    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < n; i++ )
        {
            b[i+j*n] = a[i+j*n];
        }
    }
    
    for ( j = 0; j < n; j++ )
    {
        for ( k = 0; k <= j-1; k++ )
        {
            for ( i = 0; i <= k-1; i++ )
            {
                b[k+j*n] = b[k+j*n] - b[i+k*n] * b[i+j*n];
            }
            b[k+j*n] = b[k+j*n] / b[k+k*n];
        }
        
        s = b[j+j*n];
        for ( i = 0; i <= j-1; i++ )
        {
            s = s - b[i+j*n] * b[i+j*n];
        }
        
        if ( s <= 0.0 )
        {
            b.clear();
            return b;
        }
        
        b[j+j*n] = sqrt ( s );
    }
    //
    //  Since the Cholesky factor is in R8GE format, zero out the lower triangle.
    //
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < i; j++ )
        {
            b[i+j*n] = 0.0;
        }
    }
    
    return b;
}
//****************************************************************************80

vector<double> r8vec_normal_01_new ( int n, int seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
//
//
//  Parameters:
//
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, real R(N+1), is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int SAVED, is 0 or 1 depending on whether there is a
//    single saved value left over from the previous call.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.  This starts off as 1:N, but is adjusted
//    if we have a saved value that can be immediately stored in X(1),
//    and so on.
//
//    Local, real Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
# define PI 3.141592653589793
    
    int i;
    int m;
    static int made = 0;
    vector<double> r;
    static int saved = 0;
    int x_hi;
    int x_lo;
    static double y = 0.0;
    
    vector<double> x(n);
    //
    //  Record the range of X we need to fill in.
    //
    x_lo = 1;
    x_hi = n;
    //
    //  Use up the old value, if we have it.
    //
    if ( saved == 1 )
    {
        x[0] = y;
        saved = 0;
        x_lo = 2;
    }
    //
    //  Maybe we don't need any more values.
    //
    if ( x_hi - x_lo + 1 == 0 )
    {
    }
    //
    //  If we need just one new value, do that here to avoid null arrays.
    //
    else if ( x_hi - x_lo + 1 == 1 )
    {
        r = r8vec_uniform_01_new ( 2, seed );
        
        x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * PI * r[1] );
        y =         sqrt ( - 2.0 * log ( r[0] ) ) * sin ( 2.0 * PI * r[1] );
        
        saved = 1;
        
        made = made + 2;
    }
    //
    //  If we require an even number of values, that's easy.
    //
    else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
    {
        m = ( x_hi - x_lo + 1 ) / 2;
        
        r = r8vec_uniform_01_new ( 2*m, seed );
        
        for ( i = 0; i <= 2*m-2; i = i + 2 )
        {
            x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
            x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );
        }
        made = made + x_hi - x_lo + 1;
    }
    //
    //  If we require an odd number of values, we generate an even number,
    //  and handle the last pair specially, storing one in X(N), and
    //  saving the other for later.
    //
    else
    {
        x_hi = x_hi - 1;
        
        m = ( x_hi - x_lo + 1 ) / 2 + 1;
        
        r = r8vec_uniform_01_new ( 2*m, seed );
        
        for ( i = 0; i <= 2*m-4; i = i + 2 )
        {
            x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
            x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );
        }
        
        i = 2*m - 2;
        
        x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
        y           = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );
        
        saved = 1;
        
        made = made + x_hi - x_lo + 2;
    }
    
    return x;
# undef PI
}

vector<double> r8vec_uniform_01_new ( int n, int seed )
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
    int i;
    int i4_huge = 2147483647;
    int k;
    vector<double> r(n);
    
    for ( i = 0; i < n; i++ )
    {
        k = seed / 127773;
        
        seed = 16807 * ( seed - k * 127773 ) - k * 2836;
        
        if ( seed < 0 )
        {
            seed = seed + i4_huge;
        }
        
        r[i] = ( double ) ( seed ) * 4.656612875E-10;
    }
    
    return r;
}


#endif /* MultiNormalDist_h */
