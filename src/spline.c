# include <cstdlib>
# include <iostream>
# include <iomanip>
//# include <cmath>
# include <ctime>
# include <math.h>
# include <string.h>

using namespace std;

# include "spline.h"

//******************************************************************************

double basis_function_b_val ( double tdata[], double tval )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_FUNCTION_B_VAL evaluates the B spline basis function.
//
//  Discussion:
//
//    The B spline basis function is a piecewise cubic which
//    has the properties that:
//
//    * it equals 2/3 at TDATA(3), 1/6 at TDATA(2) and TDATA(4);
//    * it is 0 for TVAL <= TDATA(1) or TDATA(5) <= TVAL;
//    * it is strictly increasing from TDATA(1) to TDATA(3),
//      and strictly decreasing from TDATA(3) to TDATA(5);
//    * the function and its first two derivatives are continuous
//      at each node TDATA(I).
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Davies and Philip Samuels,
//    An Introduction to Computational Geometry for Curves and Surfaces,
//    Clarendon Press, 1996.
//
//  Parameters:
//
//    Input, double TDATA(5), the nodes associated with the basis function.
//    The entries of TDATA are assumed to be distinct and increasing.
//
//    Input, double TVAL, a point at which the B spline basis function is
//    to be evaluated.
//
//    Output, double BASIS_FUNCTION_B_VAL, the value of the function at TVAL.
//
{
# define NDATA 5

    int left;
    int right;
    double u;
    double yval = 0.0;
    //
    if ( tval <= tdata[0] || tdata[NDATA-1] <= tval ) {
        yval = 0.0E+00;
        return yval;
    }
    //
    //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
    //
    dvec_bracket ( NDATA, tdata, tval, &left, &right );
    //
    //  U is the normalized coordinate of TVAL in this interval.
    //
    u = ( tval - tdata[left-1] ) / ( tdata[right-1] - tdata[left-1] );
    //
    //  Now evaluate the function.
    //
    if ( tval < tdata[1] ) {
        yval = pow ( u, 3 ) / 6.0E+00;
    } else if ( tval < tdata[2] ) {
        yval = ( ( (     - 3.0E+00
                         * u + 3.0E+00 )
                   * u + 3.0E+00 )
                 * u + 1.0E+00 ) / 6.0E+00;
    } else if ( tval < tdata[3]) {
        yval = ( ( (     + 3.0E+00
                         * u - 6.0E+00 )
                   * u + 0.0E+00 )
                 * u + 4.0E+00 ) / 6.0E+00;
    } else if ( tval < tdata[4] ) {
        yval = pow ( ( 1.0E+00 - u ), 3 ) / 6.0E+00;
    }

    return yval;

# undef NDATA
}
//******************************************************************************

double basis_function_beta_val ( double beta1, double beta2, double tdata[], double tval )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_FUNCTION_BETA_VAL evaluates the beta spline basis function.
//
//  Discussion:
//
//    With BETA1 = 1 and BETA2 = 0, the beta spline basis function
//    equals the B spline basis function.
//
//    With BETA1 large, and BETA2 = 0, the beta spline basis function
//    skews to the right, that is, its maximum increases, and occurs
//    to the right of the center.
//
//    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
//    a linear basis function; that is, its value in the outer two intervals
//    goes to zero, and its behavior in the inner two intervals approaches
//    a piecewise linear function that is 0 at the second node, 1 at the
//    third, and 0 at the fourth.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Davies and Philip Samuels,
//    An Introduction to Computational Geometry for Curves and Surfaces,
//    Clarendon Press, 1996, page 129.
//
//  Parameters:
//
//    Input, double BETA1, the skew or bias parameter.
//    BETA1 = 1 for no skew or bias.
//
//    Input, double BETA2, the tension parameter.
//    BETA2 = 0 for no tension.
//
//    Input, double TDATA[5], the nodes associated with the basis function.
//    The entries of TDATA are assumed to be distinct and increasing.
//
//    Input, double TVAL, a point at which the B spline basis function is
//    to be evaluated.
//
//    Output, double BASIS_FUNCTION_BETA_VAL, the value of the function at TVAL.
//
{
# define NDATA 5

    double a;
    double b;
    double c;
    double d;
    int left;
    int right;
    double u;
    double yval = 0.0;
    //
    if ( tval <= tdata[0] || tdata[NDATA-1] <= tval ) {
        yval = 0.0E+00;
        return yval;
    }
    //
    //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
    //
    dvec_bracket ( NDATA, tdata, tval, &left, &right );
    //
    //  U is the normalized coordinate of TVAL in this interval.
    //
    u = ( tval - tdata[left-1] ) / ( tdata[right-1] - tdata[left-1] );
    //
    //  Now evaluate the function.
    //
    if ( tval < tdata[1] ) {
        yval = 2.0E+00 * u * u * u;
    } else if ( tval < tdata[2] ) {
        a = beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1 * beta1
            + 6.0E+00 * ( 1.0E+00 - beta1 * beta1 )
            - 3.0E+00 * ( 2.0E+00 + beta2 + 2.0E+00 * beta1 )
            + 2.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1 * beta1 );

        b = - 6.0E+00 * ( 1.0E+00 - beta1 * beta1 )
            + 6.0E+00 * ( 2.0E+00 + beta2 + 2.0E+00 * beta1 )
            - 6.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1 * beta1 );

        c = - 3.0E+00 * ( 2.0E+00 + beta2 + 2.0E+00 * beta1 )
            + 6.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1 * beta1 );

        d = - 2.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1 * beta1 );

        yval = a + b * u + c * u * u + d * u * u * u;
    } else if ( tval < tdata[3] ) {
        a = beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1 * beta1;

        b = - 6.0E+00 * beta1 * ( 1.0E+00 - beta1 * beta1 );

        c = - 3.0E+00 * ( beta2 + 2.0E+00 * beta1 * beta1
                          + 2.0E+00 * beta1 * beta1 * beta1 );

        d = 2.0E+00 * ( beta2 + beta1 + beta1 * beta1 + beta1 * beta1 * beta1 );

        yval = a + b * u + c * u * u + d * u * u * u;
    } else if ( tval < tdata[4] ) {
        yval = 2.0E+00 * pow ( beta1 * ( 1.0E+00 - u ), 3 );
    }

    yval = yval / ( 2.0E+00 + beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1 * beta1
                    + 2.0E+00 * beta1 * beta1 * beta1 );

    return yval;
# undef NDATA
}
//******************************************************************************

double *basis_matrix_b_uni ( void )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_B_UNI sets up the uniform B spline basis matrix.
//
//  Modified:
//
//    11 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Foley, van Dam, Feiner, Hughes,
//    Computer Graphics: Principles and Practice,
//    page 493.
//
//  Parameters:
//
//    Output, double BASIS_MATRIX_B_UNI[4*4], the basis matrix.
//
{
    int i;
    int j;
    double *mbasis;
    double mbasis_save[4*4] = {
        -1.0E+00 / 6.0E+00,
        3.0E+00 / 6.0E+00,
        -3.0E+00 / 6.0E+00,
        1.0E+00 / 6.0E+00,
        3.0E+00 / 6.0E+00,
        -6.0E+00 / 6.0E+00,
        0.0E+00,
        4.0E+00 / 6.0E+00,
        -3.0E+00 / 6.0E+00,
        3.0E+00 / 6.0E+00,
        3.0E+00 / 6.0E+00,
        1.0E+00 / 6.0E+00,
        1.0E+00 / 6.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00
    };

    mbasis = new double[4*4];

    for ( j = 0; j < 4; j++ ) {
        for ( i = 0; i < 4; i++ ) {
            mbasis[i+j*4] = mbasis_save[i+j*4];
        }
    }

    return mbasis;
}
//******************************************************************************

double *basis_matrix_beta_uni ( double beta1, double beta2 )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_BETA_UNI sets up the uniform beta spline basis matrix.
//
//  Discussion:
//
//    If BETA1 = 1 and BETA2 = 0, then the beta spline reduces to
//    the B spline.
//
//  Modified:
//
//    12 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Foley, van Dam, Feiner, Hughes,
//    Computer Graphics: Principles and Practice,
//    page 505.
//
//  Parameters:
//
//    Input, double BETA1, the skew or bias parameter.
//    BETA1 = 1 for no skew or bias.
//
//    Input, double BETA2, the tension parameter.
//    BETA2 = 0 for no tension.
//
//    Output, double BASIS_MATRIX_BETA_UNI[4*4], the basis matrix.
//
{
    double delta;
    int i;
    int j;
    double *mbasis;

    mbasis = new double[4*4];

    mbasis[0+0*4] = - 2.0E+00 * beta1 * beta1 * beta1;
    mbasis[0+1*4] =   2.0E+00 * beta2
                      + 2.0E+00 * beta1 * ( beta1 * beta1 + beta1 + 1.0E+00 );
    mbasis[0+2*4] = - 2.0E+00 * ( beta2 + beta1 * beta1 + beta1 + 1.0E+00 );
    mbasis[0+3*4] =   2.0E+00;

    mbasis[1+0*4] =   6.0E+00 * beta1 * beta1 * beta1;
    mbasis[1+1*4] = - 3.0E+00 * beta2
                    - 6.0E+00 * beta1 * beta1 * ( beta1 + 1.0E+00 );
    mbasis[1+2*4] =   3.0E+00 * beta2 + 6.0E+00 * beta1 * beta1;
    mbasis[1+3*4] =   0.0E+00;

    mbasis[2+0*4] = - 6.0E+00 * beta1 * beta1 * beta1;
    mbasis[2+1*4] =   6.0E+00 * beta1 * ( beta1 - 1.0E+00 ) * ( beta1 + 1.0E+00 );
    mbasis[2+2*4] =   6.0E+00 * beta1;
    mbasis[2+3*4] =   0.0E+00;

    mbasis[3+0*4] =   2.0E+00 * beta1 * beta1 * beta1;
    mbasis[3+1*4] =   4.0E+00 * beta1 * ( beta1 + 1.0E+00 ) + beta2;
    mbasis[3+2*4] =   2.0E+00;
    mbasis[3+3*4] =   0.0E+00;

    delta = ( ( 2.0E+00
                * beta1 + 4.0E+00 )
              * beta1 + 4.0E+00 )
            * beta1 + 2.0E+00 + beta2;

    for ( j = 0; j < 4; j++ ) {
        for ( i = 0; i < 4; i++ ) {
            mbasis[i+j*4] = mbasis[i+j*4] / delta;
        }
    }

    return mbasis;
}
//******************************************************************************

double *basis_matrix_bezier ( void )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_BEZIER_UNI sets up the cubic Bezier spline basis matrix.
//
//  Discussion:
//
//    This basis matrix assumes that the data points are stored as
//    ( P1, P2, P3, P4 ).  P1 is the function value at T = 0, while
//    P2 is used to approximate the derivative at T = 0 by
//    dP/dt = 3 * ( P2 - P1 ).  Similarly, P4 is the function value
//    at T = 1, and P3 is used to approximate the derivative at T = 1
//    by dP/dT = 3 * ( P4 - P3 ).
//
//  Modified:
//
//    13 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Foley, van Dam, Feiner, Hughes,
//    Computer Graphics: Principles and Practice,
//    page 489.
//
//  Parameters:
//
//    Output, double BASIS_MATRIX_BEZIER[4*4], the basis matrix.
//
{
    double *mbasis;

    mbasis = new double[4*4];

    mbasis[0+0*4] = -1.0E+00;
    mbasis[0+1*4] =  3.0E+00;
    mbasis[0+2*4] = -3.0E+00;
    mbasis[0+3*4] =  1.0E+00;

    mbasis[1+0*4] =  3.0E+00;
    mbasis[1+1*4] = -6.0E+00;
    mbasis[1+2*4] =  3.0E+00;
    mbasis[1+3*4] =  0.0E+00;

    mbasis[2+0*4] = -3.0E+00;
    mbasis[2+1*4] =  3.0E+00;
    mbasis[2+2*4] =  0.0E+00;
    mbasis[2+3*4] =  0.0E+00;

    mbasis[3+0*4] =  1.0E+00;
    mbasis[3+1*4] =  0.0E+00;
    mbasis[3+2*4] =  0.0E+00;
    mbasis[3+3*4] =  0.0E+00;

    return mbasis;
}
//******************************************************************************

double *basis_matrix_hermite ( void )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_HERMITE sets up the Hermite spline basis matrix.
//
//  Discussion:
//
//    This basis matrix assumes that the data points are stored as
//    ( P1, P2, P1', P2' ), with P1 and P1' being the data value and
//    the derivative dP/dT at T = 0, while P2 and P2' apply at T = 1.
//
//  Modified:
//
//    13 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Foley, van Dam, Feiner, Hughes,
//    Computer Graphics: Principles and Practice,
//    page 484.
//
//  Parameters:
//
//    Output, double BASIS_MATRIX_HERMITE[4*4], the basis matrix.
//
{
    double *mbasis;

    mbasis = new double[4*4];

    mbasis[0+0*4] =  2.0E+00;
    mbasis[0+1*4] = -2.0E+00;
    mbasis[0+2*4] =  1.0E+00;
    mbasis[0+3*4] =  1.0E+00;

    mbasis[1+0*4] = -3.0E+00;
    mbasis[1+1*4] =  3.0E+00;
    mbasis[1+2*4] = -2.0E+00;
    mbasis[1+3*4] = -1.0E+00;

    mbasis[2+0*4] =  0.0E+00;
    mbasis[2+1*4] =  0.0E+00;
    mbasis[2+2*4] =  1.0E+00;
    mbasis[2+3*4] =  0.0E+00;

    mbasis[3+0*4] =  1.0E+00;
    mbasis[3+1*4] =  0.0E+00;
    mbasis[3+2*4] =  0.0E+00;
    mbasis[3+3*4] =  0.0E+00;

    return mbasis;
}
//******************************************************************************

double *basis_matrix_overhauser_nonuni ( double alpha, double beta )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_OVERHAUSER_NONUNI sets up the nonuniform Overhauser spline basis matrix.
//
//  Discussion:
//
//    This basis matrix assumes that the data points P1, P2, P3 and
//    P4 are not uniformly spaced in T, and that P2 corresponds to T = 0,
//    and P3 to T = 1.
//
//  Modified:
//
//    13 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, BETA.
//    ALPHA = | P2 - P1 | / ( | P3 - P2 | + | P2 - P1 | )
//    BETA  = | P3 - P2 | / ( | P4 - P3 | + | P3 - P2 | ).
//
//    Output, double BASIS_MATRIX_OVERHAUSER_NONUNI[4*4], the basis matrix.
//
{
    double *mbasis;

    mbasis = new double[4*4];

    mbasis[0+0*4] = - ( 1.0E+00 - alpha ) * ( 1.0E+00 - alpha ) / alpha;
    mbasis[0+1*4] =   beta + ( 1.0E+00 - alpha ) / alpha;
    mbasis[0+2*4] =   alpha - 1.0E+00 / ( 1.0E+00 - beta );
    mbasis[0+3*4] =   beta * beta / ( 1.0E+00 - beta );

    mbasis[1+0*4] =   2.0E+00 * ( 1.0E+00 - alpha ) * ( 1.0E+00 - alpha ) / alpha;
    mbasis[1+1*4] = ( - 2.0E+00 * ( 1.0E+00 - alpha ) - alpha * beta ) / alpha;
    mbasis[1+2*4] = ( 2.0E+00 * ( 1.0E+00 - alpha )
                      - beta * ( 1.0E+00 - 2.0E+00 * alpha ) ) / ( 1.0E+00 - beta );
    mbasis[1+3*4] = - beta * beta / ( 1.0E+00 - beta );

    mbasis[2+0*4] = - ( 1.0E+00 - alpha ) * ( 1.0E+00 - alpha ) / alpha;
    mbasis[2+1*4] =   ( 1.0E+00 - 2.0E+00 * alpha ) / alpha;
    mbasis[2+2*4] =   alpha;
    mbasis[2+3*4] =   0.0E+00;

    mbasis[3+0*4] =   0.0E+00;
    mbasis[3+1*4] =   1.0E+00;
    mbasis[3+2*4] =   0.0E+00;
    mbasis[3+3*4] =   0.0E+00;

    return mbasis;
}
//******************************************************************************

double *basis_matrix_overhauser_nul ( double alpha )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_OVERHAUSER_NUL sets up the nonuniform left Overhauser spline basis matrix.
//
//  Discussion:
//
//    This basis matrix assumes that the data points P1, P2, and
//    P3 are not uniformly spaced in T, and that P1 corresponds to T = 0,
//    and P2 to T = 1. (???)
//
//  Modified:
//
//    13 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA.
//    ALPHA = | P2 - P1 | / ( | P3 - P2 | + | P2 - P1 | )
//
//    Output, double BASIS_MATRIX_OVERHAUSER_NUL[3*3], the basis matrix.
//
{
    double *mbasis;

    mbasis = new double[3*3];

    mbasis[0+0*3] =   1.0E+00 / alpha;
    mbasis[0+1*3] = - 1.0E+00 / ( alpha * ( 1.0E+00 - alpha ) );
    mbasis[0+2*3] =   1.0E+00 / ( 1.0E+00 - alpha );

    mbasis[1+0*3] = - ( 1.0E+00 + alpha ) / alpha;
    mbasis[1+1*3] =   1.0E+00 / ( alpha * ( 1.0E+00 - alpha ) );
    mbasis[1+2*3] = - alpha / ( 1.0E+00 - alpha );

    mbasis[2+0*3] =   1.0E+00;
    mbasis[2+1*3] =   0.0E+00;
    mbasis[2+2*3] =   0.0E+00;

    return mbasis;
}
//******************************************************************************

double *basis_matrix_overhauser_nur ( double beta )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_OVERHAUSER_NUR sets up the nonuniform right Overhauser spline basis matrix.
//
//  Discussion:
//
//    This basis matrix assumes that the data points PN-2, PN-1, and
//    PN are not uniformly spaced in T, and that PN-1 corresponds to T = 0,
//    and PN to T = 1. (???)
//
//  Modified:
//
//    14 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double BETA.
//    BETA = | P(N) - P(N-1) | / ( | P(N) - P(N-1) | + | P(N-1) - P(N-2) | )
//
//    Output, double BASIS_MATRIX_OVERHAUSER_NUR[3*3], the basis matrix.
//
{
    double *mbasis;

    mbasis = new double[3*3];

    mbasis[0+0*3] =   1.0E+00 / beta;
    mbasis[0+1*3] = - 1.0E+00 / ( beta * ( 1.0E+00 - beta ) );
    mbasis[0+2*3] =   1.0E+00 / ( 1.0E+00 - beta );

    mbasis[1+0*3] = - ( 1.0E+00 + beta ) / beta;
    mbasis[1+1*3] =   1.0E+00 / ( beta * ( 1.0E+00 - beta ) );
    mbasis[1+2*3] = - beta / ( 1.0E+00 - beta );

    mbasis[2+0*3] =   1.0E+00;
    mbasis[2+1*3] =   0.0E+00;
    mbasis[2+2*3] =   0.0E+00;

    return mbasis;
}
//******************************************************************************

double *basis_matrix_overhauser_uni ( void)

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_OVERHAUSER_UNI sets up the uniform Overhauser spline basis matrix.
//
//  Discussion:
//
//    This basis matrix assumes that the data points P1, P2, P3 and
//    P4 are uniformly spaced in T, and that P2 corresponds to T = 0,
//    and P3 to T = 1.
//
//  Modified:
//
//    14 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Foley, van Dam, Feiner, Hughes,
//    Computer Graphics: Principles and Practice,
//    page 505.
//
//  Parameters:
//
//    Output, double BASIS_MATRIX_OVERHASUER_UNI[4*4], the basis matrix.
//
{
    double *mbasis;

    mbasis = new double[4*4];

    mbasis[0+0*4] = - 1.0E+00 / 2.0E+00;
    mbasis[0+1*4] =   3.0E+00 / 2.0E+00;
    mbasis[0+2*4] = - 3.0E+00 / 2.0E+00;
    mbasis[0+3*4] =   1.0E+00 / 2.0E+00;

    mbasis[1+0*4] =   2.0E+00 / 2.0E+00;
    mbasis[1+1*4] = - 5.0E+00 / 2.0E+00;
    mbasis[1+2*4] =   4.0E+00 / 2.0E+00;
    mbasis[1+3*4] = - 1.0E+00 / 2.0E+00;

    mbasis[2+0*4] = - 1.0E+00 / 2.0E+00;
    mbasis[2+1*4] =   0.0E+00;
    mbasis[2+2*4] =   1.0E+00 / 2.0E+00;
    mbasis[2+3*4] =   0.0E+00;

    mbasis[3+0*4] =   0.0E+00;
    mbasis[3+1*4] =   2.0E+00 / 2.0E+00;
    mbasis[3+2*4] =   0.0E+00;
    mbasis[3+3*4] =   0.0E+00;

    return mbasis;
}
//******************************************************************************

double *basis_matrix_overhauser_uni_l ( void )

//******************************************************************************
//
//  Purpose: BASIS_MATRIX_OVERHAUSER_UNI_L sets up the left uniform Overhauser spline basis matrix.
//
//  Discussion:
//
//    This basis matrix assumes that the data points P1, P2, and P3
//    are not uniformly spaced in T, and that P1 corresponds to T = 0,
//    and P2 to T = 1.
//
//  Modified:
//
//    14 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double BASIS_MATRIX_OVERHASUER_UNI_L[3*3], the basis matrix.
//
{
    double *mbasis;

    mbasis = new double[3*3];

    mbasis[0+0*3] =   2.0E+00;
    mbasis[0+1*3] = - 4.0E+00;
    mbasis[0+2*3] =   2.0E+00;

    mbasis[1+0*3] = - 3.0E+00;
    mbasis[1+1*3] =   4.0E+00;
    mbasis[1+2*3] = - 1.0E+00;

    mbasis[2+0*3] =   1.0E+00;
    mbasis[2+1*3] =   0.0E+00;
    mbasis[2+2*3] =   0.0E+00;

    return mbasis;
}
//******************************************************************************

double *basis_matrix_overhauser_uni_r ( void )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_OVERHAUSER_UNI_R sets up the right uniform Overhauser spline basis matrix.
//
//  Discussion:
//
//    This basis matrix assumes that the data points P(N-2), P(N-1),
//    and P(N) are uniformly spaced in T, and that P(N-1) corresponds to
//    T = 0, and P(N) to T = 1.
//
//  Modified:
//
//    14 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double BASIS_MATRIX_OVERHASUER_UNI_R[3*3], the basis matrix.
//
{
    double *mbasis;

    mbasis = new double[3*3];

    mbasis[0+0*3] =   2.0E+00;
    mbasis[0+1*3] = - 4.0E+00;
    mbasis[0+2*3] =   2.0E+00;

    mbasis[1+0*3] = - 3.0E+00;
    mbasis[1+1*3] =   4.0E+00;
    mbasis[1+2*3] = - 1.0E+00;

    mbasis[2+0*3] =   1.0E+00;
    mbasis[2+1*3] =   0.0E+00;
    mbasis[2+2*3] =   0.0E+00;

    return mbasis;
}
//******************************************************************************

double basis_matrix_tmp ( int left, int n, double mbasis[], int ndata,
                          double tdata[], double ydata[], double tval )

//******************************************************************************
//
//  Purpose:
//
//    BASIS_MATRIX_TMP computes Q = T * MBASIS * P
//
//  Discussion:
//
//    YDATA is a vector of data values, most frequently the values of some
//    function sampled at uniformly spaced points.  MBASIS is the basis
//    matrix for a particular kind of spline.  T is a vector of the
//    powers of the normalized difference between TVAL and the left
//    endpoint of the interval.
//
//  Modified:
//
//    14 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LEFT, indicats that TVAL is in the interval
//    [ TDATA(LEFT), TDATA(LEFT+1) ], or that this is the "nearest"
//    interval to TVAL.
//    For TVAL < TDATA(1), use LEFT = 1.
//    For TDATA(NDATA) < TVAL, use LEFT = NDATA - 1.
//
//    Input, int N, the order of the basis matrix.
//
//    Input, double MBASIS[N*N], the basis matrix.
//
//    Input, int NDATA, the dimension of the vectors TDATA and YDATA.
//
//    Input, double TDATA[NDATA], the abscissa values.  This routine
//    assumes that the TDATA values are uniformly spaced, with an
//    increment of 1.0.
//
//    Input, double YDATA[NDATA], the data values to be interpolated or
//    approximated.
//
//    Input, double TVAL, the value of T at which the spline is to be
//    evaluated.
//
//    Output, double BASIS_MATRIX_TMP, the value of the spline at TVAL.
//
{
    double arg = 0.0;
    int first = 0;
    int i;
    int j;
    //  double temp;
    double tm;
    double *tvec;
    double yval;
    //
    tvec = new double[n];

    if ( left == 1 ) {
        arg = 0.5E+00 * ( tval - tdata[left-1] );
        first = left;
    } else if ( left < ndata - 1 ) {
        arg = tval - tdata[left-1];
        first = left - 1;
    } else if ( left == ndata - 1 ) {
        arg = 0.5E+00 * ( 1.0E+00 + tval - tdata[left-1] );
        first = left - 1;
    }
    //
    //  TVEC(I) = ARG**(N-I).
    //
    tvec[n-1] = 1.0E+00;
    for ( i = n-2; 0 <= i; i-- ) {
        tvec[i] = arg * tvec[i+1];
    }

    yval = 0.0E+00;
    for ( j = 0; j < n; j++ ) {
        tm = 0.0E+00;
        for ( i = 0; i < n; i++ ) {
            tm = tm + tvec[i] * mbasis[i+j*n];
        }
        yval = yval + tm * ydata[first - 1 + j];
    }

    delete [] tvec;

    return yval;
}
//******************************************************************************

void bc_val ( int n, double t, double xcon[], double ycon[], double *xval,
              double *yval )

//******************************************************************************
//
//  Purpose:
//
//    BC_VAL evaluates a parameterized Bezier curve.
//
//  Discussion:
//
//    BC_VAL(T) is the value of a vector function of the form
//
//      BC_VAL(T) = ( X(T), Y(T) )
//
//    where
//
//      X(T) = Sum ( 0 <= I <= N ) XCON(I) * BERN(I,N)(T)
//      Y(T) = Sum ( 0 <= I <= N ) YCON(I) * BERN(I,N)(T)
//
//    BERN(I,N)(T) is the I-th Bernstein polynomial of order N
//    defined on the interval [0,1],
//
//    XCON(0:N) and YCON(0:N) are the coordinates of N+1 "control points".
//
//  Modified:
//
//    12 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kahaner, Moler, and Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989.
//
//  Parameters:
//
//    Input, int N, the order of the Bezier curve, which
//    must be at least 0.
//
//    Input, double T, the point at which the Bezier curve should
//    be evaluated.  The best results are obtained within the interval
//    [0,1] but T may be anywhere.
//
//    Input, double XCON[0:N], YCON[0:N], the X and Y coordinates
//    of the control points.  The Bezier curve will pass through
//    the points ( XCON(0), YCON(0) ) and ( XCON(N), YCON(N) ), but
//    generally NOT through the other control points.
//
//    Output, double *XVAL, *YVAL, the X and Y coordinates of the point
//    on the Bezier curve corresponding to the given T value.
//
{
    double *bval;
    int i;
    //
    bval = bp01 ( n, t );

    *xval = 0.0E+00;
    for ( i = 0; i <= n; i++ ) {
        *xval = *xval + xcon[i] * bval[i];
    }

    *yval = 0.0E+00;
    for ( i = 0; i <= n; i++ ) {
        *yval = *yval + ycon[i] * bval[i];
    }

    delete [] bval;

    return;
}
//******************************************************************************

double bez_val ( int n, double x, double a, double b, double y[] )

//******************************************************************************
//
//  Purpose:
//
//    BEZ_VAL evaluates a Bezier function at a point.
//
//  Discussion:
//
//    The Bezier function has the form:
//
//      BEZ(X) = Sum ( 0 <= I <= N ) Y(I) * BERN(N,I)( (X-A)/(B-A) )
//
//    BERN(N,I)(X) is the I-th Bernstein polynomial of order N
//    defined on the interval [0,1],
//
//    Y(0:N) is a set of coefficients,
//
//    and if, for I = 0 to N, we define the N+1 points
//
//      X(I) = ( (N-I)*A + I*B) / N,
//
//    equally spaced in [A,B], the pairs ( X(I), Y(I) ) can be regarded as
//    "control points".
//
//  Modified:
//
//    12 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kahaner, Moler, and Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989.
//
//  Parameters:
//
//    Input, int N, the order of the Bezier function, which
//    must be at least 0.
//
//    Input, double X, the point at which the Bezier function should
//    be evaluated.  The best results are obtained within the interval
//    [A,B] but X may be anywhere.
//
//    Input, double A, B, the interval over which the Bezier function
//    has been defined.  This is the interval in which the control
//    points have been set up.  Note BEZ(A) = Y(0) and BEZ(B) = Y(N),
//    although BEZ will not, in general pass through the other
//    control points.  A and B must not be equal.
//
//    Input, double Y[0:N], a set of data defining the Y coordinates
//    of the control points.
//
//    Output, double BEZ_VAL, the value of the Bezier function at X.
//
{
    double *bval;
    int i;
    double value;
    double x01;
    //
    if ( b - a == 0.0E+00 ) {
        cout << "\n";
        cout << "BEZ_VAL - Fatal error!\n";
        cout << "  Null interval, A = B = " << a << "\n";
        exit ( 1 );
    }
    //
    //  X01 lies in [0,1], in the same relative position as X in [A,B].
    //
    x01 = ( x - a ) / ( b - a );

    bval = bp01 ( n, x01 );

    value = 0.0E+00;
    for ( i = 0; i <= n; i++ ) {
        value = value + y[i] * bval[i];
    }

    delete [] bval;

    return value;
}
//******************************************************************************

double bp_approx ( int n, double a, double b, double ydata[], double xval )

//******************************************************************************
//
//  Purpose:
//
//    BP_APPROX evaluates the Bernstein polynomial for F(X) on [A,B].
//
//  Formula:
//
//    BERN(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
//
//    where
//
//      X(I) = ( ( N - I ) * A + I * B ) / N
//      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
//
//  Discussion:
//
//    The Bernstein polynomial BERN(F) for F(X) is an approximant, not an
//    interpolant; in other words, its value is not guaranteed to equal
//    that of F at any particular point.  However, for a fixed interval
//    [A,B], if we let N increase, the Bernstein polynomial converges
//    uniformly to F everywhere in [A,B], provided only that F is continuous.
//    Even if F is not continuous, but is bounded, the polynomial converges
//    pointwise to F(X) at all points of continuity.  On the other hand,
//    the convergence is quite slow compared to other interpolation
//    and approximation schemes.
//
//  Modified:
//
//    12 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kahaner, Moler, and Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989.
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein polynomial to be used.
//
//    Input, double A, B, the endpoints of the interval on which the
//    approximant is based.  A and B should not be equal.
//
//    Input, double YDATA[0:N], the data values at N+1 equally spaced points
//    in [A,B].  If N = 0, then the evaluation point should be 0.5 * ( A + B).
//    Otherwise, evaluation point I should be ( (N-I)*A + I*B ) / N ).
//
//    Input, double XVAL, the point at which the Bernstein polynomial
//    approximant is to be evaluated.  XVAL does not have to lie in the
//    interval [A,B].
//
//    Output, double BP_APPROX, the value of the Bernstein polynomial approximant
//    for F, based in [A,B], evaluated at XVAL.
//
{
    double *bvec;
    int i;
    double yval;
    //
    //  Evaluate the Bernstein basis polynomials at XVAL.
    //
    bvec = bpab ( n, a, b, xval );
    //
    //  Now compute the sum of YDATA(I) * BVEC(I).
    //
    yval = 0.0E+00;

    for ( i = 0; i <= n; i++ ) {
        yval = yval + ydata[i] * bvec[i];
    }
    delete [] bvec;

    return yval;
}
//******************************************************************************

double *bp01 ( int n, double x )

//******************************************************************************
//
//  Purpose:
//
//    BP01 evaluates the Bernstein basis polynomials for [0,1] at a point.
//
//  Discussion:
//
//    For any N greater than or equal to 0, there is a set of N+1 Bernstein
//    basis polynomials, each of degree N, which form a basis for
//    all polynomials of degree N on [0,1].
//
//  Formula:
//
//    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)**(N-I) * X**I
//
//    N is the degree;
//
//    0 <= I <= N indicates which of the N+1 basis polynomials
//    of degree N to choose;
//
//    X is a point in [0,1] at which to evaluate the basis polynomial.
//
//  First values:
//
//    B(0,0,X) = 1
//
//    B(1,0,X) =      1-X
//    B(1,1,X) =                X
//
//    B(2,0,X) =     (1-X)**2
//    B(2,1,X) = 2 * (1-X)    * X
//    B(2,2,X) =                X**2
//
//    B(3,0,X) =     (1-X)**3
//    B(3,1,X) = 3 * (1-X)**2 * X
//    B(3,2,X) = 3 * (1-X)    * X**2
//    B(3,3,X) =                X**3
//
//    B(4,0,X) =     (1-X)**4
//    B(4,1,X) = 4 * (1-X)**3 * X
//    B(4,2,X) = 6 * (1-X)**2 * X**2
//    B(4,3,X) = 4 * (1-X)    * X**3
//    B(4,4,X) =                X**4
//
//  Special values:
//
//    B(N,I,1/2) = C(N,K) / 2**N
//
//  Modified:
//
//    12 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kahaner, Moler, and Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989.
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein basis polynomials.
//
//    Input, double X, the evaluation point.
//
//    Output, double BP01[0:N], the values of the N+1 Bernstein basis
//    polynomials at X.
//
{
    double *bern;
    int i;
    int j;
    //
    bern = new double[n+1];

    if ( n == 0 ) {
        bern[0] = 1.0E+00;
    } else if ( 0 < n ) {
        bern[0] = 1.0E+00 - x;
        bern[1] = x;

        for ( i = 2; i <= n; i++ ) {
            bern[i] = x * bern[i-1];
            for ( j = i-1; 1 <= j; j-- ) {
                bern[j] = x * bern[j-1] + ( 1.0E+00 - x ) * bern[j];
            }
            bern[0] = ( 1.0E+00 - x ) * bern[0];
        }

    }

    return bern;
}
//******************************************************************************

double *bpab ( int n, double a, double b, double x )

//******************************************************************************
//
//  Purpose:
//
//    BPAB evaluates the Bernstein basis polynomials for [A,B] at a point.
//
//  Formula:
//
//    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)**(N-I) * (X-A)**I / (B-A)**N
//
//  First values:
//
//    B(0,0,X) =   1
//
//    B(1,0,X) = (      B-X                ) / (B-A)
//    B(1,1,X) = (                 X-A     ) / (B-A)
//
//    B(2,0,X) = (     (B-X)**2            ) / (B-A)**2
//    B(2,1,X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)**2
//    B(2,2,X) = (                (X-A)**2 ) / (B-A)**2
//
//    B(3,0,X) = (     (B-X)**3            ) / (B-A)**3
//    B(3,1,X) = ( 3 * (B-X)**2 * (X-A)    ) / (B-A)**3
//    B(3,2,X) = ( 3 * (B-X)    * (X-A)**2 ) / (B-A)**3
//    B(3,3,X) = (                (X-A)**3 ) / (B-A)**3
//
//    B(4,0,X) = (     (B-X)**4            ) / (B-A)**4
//    B(4,1,X) = ( 4 * (B-X)**3 * (X-A)    ) / (B-A)**4
//    B(4,2,X) = ( 6 * (B-X)**2 * (X-A)**2 ) / (B-A)**4
//    B(4,3,X) = ( 4 * (B-X)    * (X-A)**3 ) / (B-A)**4
//    B(4,4,X) = (                (X-A)**4 ) / (B-A)**4
//
//  Modified:
//
//    12 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kahaner, Moler, and Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989.
//
//  Parameters:
//
//    Input, integer N, the degree of the Bernstein basis polynomials.
//    For any N greater than or equal to 0, there is a set of N+1
//    Bernstein basis polynomials, each of degree N, which form a basis
//    for polynomials on [A,B].
//
//    Input, double A, B, the endpoints of the interval on which the
//    polynomials are to be based.  A and B should not be equal.
//
//    Input, double X, the point at which the polynomials are to be
//    evaluated.  X need not lie in the interval [A,B].
//
//    Output, double BERN[0:N], the values of the N+1 Bernstein basis
//    polynomials at X.
//
{
    double *bern;
    int i;
    int j;
    //
    if ( b == a ) {
        cout << "\n";
        cout << "BPAB - Fatal error!\n";
        cout << "  A = B = " << a << "\n";
        exit ( 1 );
    }

    bern = new double[n+1];

    if ( n == 0 ) {
        bern[0] = 1.0E+00;
    } else if ( 0 < n ) {
        bern[0] = ( b - x ) / ( b - a );
        bern[1] = ( x - a ) / ( b - a );

        for ( i = 2; i <= n; i++ ) {
            bern[i] = ( x - a ) * bern[i-1] / ( b - a );
            for ( j = i-1; 1 <= j; j-- ) {
                bern[j] = ( ( b - x ) * bern[j] + ( x - a ) * bern[j-1] ) / ( b - a );
            }
            bern[0] = ( b - x ) * bern[0] / ( b - a );
        }
    }

    return bern;
}
//*********************************************************************

double d_max ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    D_MAX returns the maximum of two double precision values.
//
//  Modified:
//
//    10 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double D_MAX, the maximum of X and Y.
//
{
    if ( y < x ) {
        return x;
    } else {
        return y;
    }
}
//*********************************************************************

double d_min ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    D_MIN returns the minimum of two double precision values.
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double D_MIN, the minimum of X and Y.
//
{
    if ( y < x ) {
        return y;
    } else {
        return x;
    }
}
//****************************************************************************

double d_random ( double rlo, double rhi, int *seed )

//****************************************************************************
//
//  Purpose:
//
//    D_RANDOM returns a random double precision value.
//
//  Modified:
//
//    04 February 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RLO, RHI, the minimum and maximum values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double D_RANDOM, the randomly chosen value.
//
{
    double t;

    t = d_uniform_01 ( seed );

    return ( 1.0E+00 - t ) * rlo + t * rhi;

}
//******************************************************************************

double d_uniform_01 ( int *seed )

//******************************************************************************
//
//  Purpose:
//
//    D_UNIFORM_01 is a portable pseudorandom number generator.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Modified:
//
//    11 August 2004
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double D_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
    int k;
    double r;

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 ) {
        *seed = *seed + 2147483647;
    }
    //
    //  Although SEED can be represented exactly as a 32 bit integer,
    //  it generally cannot be represented exactly as a 32 bit real number!
    //
    r = ( double ) ( *seed ) * 4.656612875E-10;

    return r;
}
//******************************************************************************

double *d3_mxv ( int n, double a[], double x[] )

//******************************************************************************
//
//  Purpose:
//
//    D3_MXV multiplies a D3 matrix times a vector.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Modified:
//
//    15 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, double A[3*N], the D3 matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double D3_MXV[N], the product A * x.
//
{
    double *b;
    int i;

    b = new double[n];

    for ( i = 0; i < n; i++ ) {
        b[i] =        a[1+i*3] * x[i];
    }
    for ( i = 0; i < n-1; i++ ) {
        b[i] = b[i] + a[0+(i+1)*3] * x[i+1];
    }
    for ( i = 1; i < n; i++ ) {
        b[i] = b[i] + a[2+(i-1)*3] * x[i-1];
    }

    return b;
}
//**********************************************************************

double *d3_np_fs ( int n, double a[], double b[] )

//**********************************************************************
//
//  Purpose:
//
//    D3_NP_FS factors and solves a D3 system.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    This algorithm requires that each diagonal entry be nonzero.
//    It does not use pivoting, and so can fail on systems that
//    are actually nonsingular.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Modified:
//
//    15 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, double A[3*N].
//    On input, the nonzero diagonals of the linear system.
//    On output, the data in these vectors has been overwritten
//    by factorization information.
//
//    Input, double B[N], the right hand side.
//
//    Output, double D3_NP_FS[N], the solution of the linear system.
//    This is NULL if there was an error because one of the diagonal
//    entries was zero.
//
{
    int i;
    double *x;
    double xmult;
    //
    //  Check.
    //
    for ( i = 0; i < n; i++ ) {
        if ( a[1+i*3] == 0.0E+00 ) {
            return NULL;
        }
    }
    x = new double[n];

    for ( i = 0; i < n; i++ ) {
        x[i] = b[i];
    }

    for ( i = 1; i < n; i++ ) {
        xmult = a[2+(i-1)*3] / a[1+(i-1)*3];
        a[1+i*3] = a[1+i*3] - xmult * a[0+i*3];
        x[i] = x[i] - xmult * x[i-1];
    }

    x[n-1] = x[n-1] / a[1+(n-1)*3];
    for ( i = n-2; 0 <= i; i-- ) {
        x[i] = ( x[i] - a[0+(i+1)*3] * x[i+1] ) / a[1+i*3];
    }

    return x;
}
//******************************************************************************

void d3_print ( int n, double a[], char *title )

//******************************************************************************
//
//  Purpose:
//
//    D3_PRINT prints a D3 matrix.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the D3 matrix.
//
//    Input, char *TITLE, a title to print.
//
{
    if ( 0 < s_len_trim ( title ) ) {
        cout << "\n";
        cout << title << "\n";
    }

    cout << "\n";

    d3_print_some ( n, a, 1, 1, n, n );

    return;
}
//******************************************************************************

void d3_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi )

//******************************************************************************
//
//  Purpose:
//
//    D3_PRINT_SOME prints some of a D3 matrix.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Modified:
//
//    05 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the D3 matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column, to be printed.
//
{
# define INCX 5

    int i;
    int i2hi;
    int i2lo;
    int inc;
    int j;
    int j2;
    int j2hi;
    int j2lo;
    //
    //  Print the columns of the matrix, in strips of 5.
    //
    for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX ) {
        j2hi = j2lo + INCX - 1;
        j2hi = i_min ( j2hi, n );
        j2hi = i_min ( j2hi, jhi );

        inc = j2hi + 1 - j2lo;

        cout << "\n";
        cout << "  Col: ";
        for ( j = j2lo; j <= j2hi; j++ ) {
            j2 = j + 1 - j2lo;
            cout << setw(7) << j << "       ";
        }

        cout << "\n";
        cout << "  Row\n";
        cout << "  ---\n";
        //
        //  Determine the range of the rows in this strip.
        //
        i2lo = i_max ( ilo, 1 );
        i2lo = i_max ( i2lo, j2lo - 1 );

        i2hi = i_min ( ihi, n );
        i2hi = i_min ( i2hi, j2hi + 1 );

        for ( i = i2lo; i <= i2hi; i++ ) {
            //
            //  Print out (up to) 5 entries in row I, that lie in the current strip.
            //
            cout << setw(6) << i << "  ";

            for ( j2 = 1; j2 <= inc; j2++ ) {
                j = j2lo - 1 + j2;

                if ( 1 < i-j || 1 < j-i ) {
                    cout << "              ";
                } else if ( j == i+1 ) {
                    cout << setw(12) << a[0+(j-1)*3] << "  ";
                } else if ( j == i ) {
                    cout << setw(12) << a[1+(j-1)*3] << "  ";
                } else if ( j == i-1 ) {
                    cout << setw(12) << a[2+(j-1)*3] << "  ";
                }

            }
            cout << "\n";
        }
    }

    cout << "\n";

    return;
# undef INCX
}
//******************************************************************************

double *d3_random ( int n, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    D3_RANDOM randomizes a D3 matrix.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double D3_RANDOM[3*N], the D3 matrix.
//
{
    double *a;
    int i;
    double *u;
    double *v;
    double *w;

    a = new double[3*n];

    u = dvec_random ( n-1, 0.0E+00, 1.0E+00, seed );
    v = dvec_random ( n,   0.0E+00, 1.0E+00, seed );
    w = dvec_random ( n-1, 0.0E+00, 1.0E+00, seed );

    a[0+0*3] = 0.0E+00;
    for ( i = 1; i < n; i++ ) {
        a[0+i*3] = u[i-1];
    }
    for ( i = 0; i < n; i++ ) {
        a[1+i*3] = v[i];
    }
    for ( i = 0; i < n-1; i++ ) {
        a[2+i*3] = w[i];
    }
    a[2+(n-1)*3] = 0.0E+00;

    delete [] u;
    delete [] v;
    delete [] w;

    return a;
}
//**********************************************************************

void data_to_dif ( int ntab, double xtab[], double ytab[], double diftab[] )

//**********************************************************************
//
//  Purpose:
//
//    DATA_TO_DIF sets up a divided difference table from raw data.
//
//  Discussion:
//
//    Space can be saved by using a single array for both the DIFTAB and
//    YTAB dummy parameters.  In that case, the difference table will
//    overwrite the Y data without interfering with the computation.
//
//  Modified:
//
//    04 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NTAB, the number of pairs of points
//    (XTAB[I],YTAB[I]) which are to be used as data.
//
//    Input, double XTAB[NTAB], the X values at which data was taken.
//    These values must be distinct.
//
//    Input, double YTAB[NTAB], the corresponding Y values.
//
//    Output, double DIFTAB[NTAB], the divided difference coefficients
//    corresponding to the input (XTAB,YTAB).
//
{
    int i;
    int j;
    //
    //  Copy the data values into DIFTAB.
    //
    for ( i = 0; i < ntab; i++ ) {
        diftab[i] = ytab[i];
    }
    //
    //  Make sure the abscissas are distinct.
    //
    for ( i = 0; i < ntab; i++ ) {
        for ( j = i+1; j < ntab; j++ ) {
            if ( xtab[i] - xtab[j] == 0.0E+00 ) {
                cout << "\n";
                cout << "DATA_TO_DIF - Fatal error!\n";
                cout << "  Two entries of XTAB are equal!\n";
                cout << "  XTAB[%d] = " << xtab[i] << "\n";
                cout << "  XTAB[%d] = " << xtab[j] << "\n";
                exit ( 1 );
            }
        }
    }
    //
    //  Compute the divided differences.
    //
    for ( i = 1; i <= ntab-1; i++ ) {
        for ( j = ntab-1; i <= j; j-- ) {
            diftab[j] = ( diftab[j] - diftab[j-1] ) / ( xtab[j] - xtab[j-i] );
        }
    }

    return;
}
//**********************************************************************

double dif_val ( int ntab, double xtab[], double diftab[], double xval )

//**********************************************************************
//
//  Purpose:
//
//    DIF_VAL evaluates a divided difference polynomial at a point.
//
//  Discussion:
//
//    DATA_TO_DIF must be called first to set up the divided difference table.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer NTAB, the number of divided difference
//    coefficients in DIFTAB, and the number of points XTAB.
//
//    Input, double XTAB[NTAB], the X values upon which the
//    divided difference polynomial is based.
//
//    Input, double DIFTAB[NTAB], the divided difference table.
//
//    Input, double XVAL, a value of X at which the polynomial
//    is to be evaluated.
//
//    Output, double DIF_VAL, the value of the polynomial at XVAL.
//
{
    int i;
    double value;

    value = diftab[ntab-1];
    for ( i = 2; i <= ntab; i++ ) {
        value = diftab[ntab-i] + ( xval - xtab[ntab-i] ) * value;
    }

    return value;
}
//********************************************************************

void dvec_bracket ( int n, double x[], double xval, int *left,
                    int *right )

//********************************************************************
//
//  Purpose:
//
//    DVEC_BRACKET searches a sorted array for successive brackets of a value.
//
//  Discussion:
//
//    If the values in the vector are thought of as defining intervals
//    on the real line, then this routine searches for the interval
//    nearest to or containing the given value.
//
//    It is always true that RIGHT = LEFT+1.
//
//    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
//      XVAL   < X[0] < X[1];
//    If X(1) <= XVAL < X[N-1], then
//      X[LEFT-1] <= XVAL < X[RIGHT-1];
//    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
//      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
//
//    For consistency, this routine computes indices RIGHT and LEFT
//    that are 1-based, although it would be more natural in C and
//    C++ to use 0-based values.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input, double X[N], an array that has been sorted into ascending order.
//
//    Input, double XVAL, a value to be bracketed.
//
//    Output, int *LEFT, *RIGHT, the results of the search.
//
{
    int i;

    for ( i = 2; i <= n - 1; i++ ) {
        if ( xval < x[i-1] ) {
            *left = i - 1;
            *right = i;
            return;
        }

    }

    *left = n - 1;
    *right = n;

    return;
}
//********************************************************************

void dvec_bracket3 ( int n, double t[], double tval, int *left )

//********************************************************************
//
//  Purpose:
//
//    DVEC_BRACKET3 finds the interval containing or nearest a given value.
//
//  Discussion:
//
//    The routine always returns the index LEFT of the sorted array
//    T with the property that either
//    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
//    *  T < T[LEFT] = T[0], or
//    *  T > T[LEFT+1] = T[N-1].
//
//    The routine is useful for interpolation problems, where
//    the abscissa must be located within an interval of data
//    abscissas for interpolation, or the "nearest" interval
//    to the (extreme) abscissa must be found so that extrapolation
//    can be carried out.
//
//    For consistency with other versions of this routine, the
//    value of LEFT is assumed to be a 1-based index.  This is
//    contrary to the typical C and C++ convention of 0-based indices.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of the input array.
//
//    Input, double T[N], an array that has been sorted into ascending order.
//
//    Input, double TVAL, a value to be bracketed by entries of T.
//
//    Input/output, int *LEFT.
//
//    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
//    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
//    is searched first, followed by the appropriate interval to the left
//    or right.  After that, a binary search is used.
//
//    On output, LEFT is set so that the interval [ T[LEFT-1], T[LEFT] ]
//    is the closest to TVAL; it either contains TVAL, or else TVAL
//    lies outside the interval [ T[0], T[N-1] ].
//
{
    int high;
    int low;
    int mid;
    //
    //  Check the input data.
    //
    if ( n < 2 ) {
        cout << "\n";
        cout << "DVEC_BRACKET3 - Fatal error!\n";
        cout << "  N must be at least 2.\n";
        exit ( 1 );
    }
    //
    //  If *LEFT is not between 1 and N-1, set it to the middle value.
    //
    if ( *left < 1 || n - 1 < *left ) {
        *left = ( n + 1 ) / 2;
    }

    //
    //  CASE 1: TVAL < T[*LEFT]:
    //  Search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-1.
    //
    if ( tval < t[*left] ) {

        if ( *left == 1 ) {
            return;
        } else if ( *left == 2 ) {
            *left = 1;
            return;
        } else if ( t[*left-2] <= tval ) {
            *left = *left - 1;
            return;
        } else if ( tval <= t[1] ) {
            *left = 1;
            return;
        }
        //
        //  ...Binary search for TVAL in (T[I-1],T[I]), for I = 2 to *LEFT-2.
        //
        low = 2;
        high = *left - 2;

        for (;;) {

            if ( low == high ) {
                *left = low;
                return;
            }

            mid = ( low + high + 1 ) / 2;

            if ( t[mid-1] <= tval ) {
                low = mid;
            } else {
                high = mid - 1;
            }

        }
    }
    //
    //  CASE 2: T[*LEFT] < TVAL:
    //  Search for TVAL in (T[I-1],T[I]) for intervals I = *LEFT+1 to N-1.
    //
    else if ( t[*left] < tval ) {

        if ( *left == n - 1 ) {
            return;
        } else if ( *left == n - 2 ) {
            *left = *left + 1;
            return;
        } else if ( tval <= t[*left+1] ) {
            *left = *left + 1;
            return;
        } else if ( t[n-2] <= tval ) {
            *left = n - 1;
            return;
        }
        //
        //  ...Binary search for TVAL in (T[I-1],T[I]) for intervals I = *LEFT+2 to N-2.
        //
        low = *left + 2;
        high = n - 2;

        for ( ; ; ) {

            if ( low == high ) {
                *left = low;
                return;
            }

            mid = ( low + high + 1 ) / 2;

            if ( t[mid-1] <= tval ) {
                low = mid;
            } else {
                high = mid - 1;
            }
        }
    }
    //
    //  CASE 3: T[*LEFT-1] <= TVAL <= T[*LEFT]:
    //  T is just where the user said it might be.
    //
    else {}

    return;
}
//******************************************************************************

double *dvec_even ( int n, double alo, double ahi )

//******************************************************************************
//
//  Purpose:
//
//    DVEC_EVEN returns N real values, evenly spaced between ALO and AHI.
//
//  Modified:
//
//    17 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, double ALO, AHI, the low and high values.
//
//    Output, double DVEC_EVEN[N], N evenly spaced values.
//    Normally, A(1) = ALO and A(N) = AHI.
//    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
//
{
    double *a;
    int i;

    a = new double[n];

    if ( n == 1 ) {
        a[0] = 0.5E+00 * ( alo + ahi );
    } else {
        for ( i = 1; i <= n; i++ ) {
            a[i-1] = ( ( double ) ( n - i     ) * alo
                       + ( double ) (     i - 1 ) * ahi )
                     / ( double ) ( n     - 1 );
        }
    }

    return a;
}
//********************************************************************

double *dvec_indicator ( int n )

//********************************************************************
//
//  Purpose:
//
//    DVEC_INDICATOR sets a real vector to the indicator vector.
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, double DVEC_INDICATOR[N], the array to be initialized.
//
{
    int i;
    double *a;

    a = new double[n];

    for ( i = 0; i < n; i++ ) {
        a[i] = ( double ) ( i + 1 );
    }

    return a;
}
//**********************************************************************

void dvec_order_type ( int n, double x[], int *order )

//**********************************************************************
//
//  Purpose:
//
//    DVEC_ORDER_TYPE determines if an array is (non)strictly ascending/descending.
//
//  Modified:
//
//    14 September 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the array.
//
//    Input, double X[N], the array to be checked.
//
//    Output, int *ORDER, order indicator:
//    -1, no discernable order;
//    0, all entries are equal;
//    1, ascending order;
//    2, strictly ascending order;
//    3, descending order;
//    4, strictly descending order.
//
{
    int i;
    //
    //  Search for the first value not equal to X(0).
    //
    i = 0;

    for (;;) {

        i = i + 1;
        if ( n-1 < i ) {
            *order = 0;
            return;
        }

        if ( x[0] < x[i] ) {
            if ( i == 1 ) {
                *order = 2;
                break;
            } else {
                *order = 1;
                break;
            }
        } else if ( x[i] < x[0] ) {
            if ( i == 1 ) {
                *order = 4;
                break;
            } else {
                *order = 3;
                break;
            }
        }
    }
    //
    //  Now we have a "direction".  Examine subsequent entries.
    //
    for (;;) {
        i = i + 1;
        if ( n - 1 < i ) {
            return;
        }

        if ( *order == 1 ) {
            if ( x[i] < x[i-1] ) {
                *order = -1;
                return;
            }
        } else if ( *order == 2 ) {
            if ( x[i] < x[i-1] ) {
                *order = -1;
                return;
            } else if ( x[i] == x[i-1] ) {
                *order = 1;
            }
        } else if ( *order == 3 ) {
            if ( x[i-1] < x[i] ) {
                *order = -1;
                return;
            }
        } else if ( *order == 4 ) {
            if ( x[i-1] < x[i] ) {
                *order = -1;
                return;
            } else if ( x[i] == x[i-1] ) {
                *order = 3;
            }

        }

    }

}
//********************************************************************

void dvec_print ( int n, double a[], char *title )

//********************************************************************
//
//  Purpose:
//
//    DVEC_PRINT prints a real vector.
//
//  Modified:
//
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, char *TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
    int i;

    if ( s_len_trim ( title ) != 0 ) {
        cout << "\n";
        cout << title << "\n";
    }

    cout << "\n";
    for ( i = 0; i <= n-1; i++ ) {
        cout << setw(6)  << i + 1 << "  "
             << setw(14) << a[i]  << "\n";
    }

    return;
}
//**********************************************************************

double *dvec_random ( int n, double alo, double ahi, int *seed )

//**********************************************************************
//
//  Purpose:
//
//    DVEC_RANDOM randomizes a real vector.
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double ALO, AHI, the range allowed for the entries.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double DVEC_RANDOM[N], the vector of randomly chosen integers.
//
{
    double *a;
    int i;

    a = new double[n];

    for ( i = 0; i < n; i++ ) {
        a[i] = d_random ( alo, ahi, seed );
    }

    return a;
}
//********************************************************************

void dvec_sort_bubble_a ( int n, double a[] )

//********************************************************************
//
//  Purpose:
//
//    DVEC_SORT_BUBBLE_A ascending sorts a real vector using bubble sort.
//
//  Modified:
//
//    09 April 1999
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input/output, double A[N].
//    On input, an unsorted array of floats.
//    On output, A has been sorted.
//
{
    int i;
    int j;
    double temp;

    for ( i = 0; i < n-1; i++ ) {
        for ( j = i+1; j < n; j++ ) {
            if ( a[j] < a[i] ) {
                temp = a[i];
                a[i] = a[j];
                a[j] = temp;
            }
        }
    }
}
//****************************************************************************

int i_max ( int i1, int i2 )

//****************************************************************************
//
//  Purpose:
//
//    I_MAX returns the maximum of two integers.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I_MAX, the larger of I1 and I2.
//
//
{
    if ( i2 < i1 ) {
        return i1;
    } else {
        return i2;
    }

}
//****************************************************************************

int i_min ( int i1, int i2 )

//****************************************************************************
//
//  Purpose:
//
//    I_MIN returns the smaller of two integers.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I_MIN, the smaller of I1 and I2.
//
//
{
    if ( i1 < i2 ) {
        return i1;
    } else {
        return i2;
    }

}
//******************************************************************************

void least_set ( int ntab, double xtab[], double ytab[], int ndeg,
                 double ptab[], double b[], double c[], double d[], double *eps, int *ierror )

//******************************************************************************
//
//  Purpose:
//
//    LEAST_SET constructs the least squares polynomial approximation to data.
//
//  Discussion:
//
//    The least squares polynomial is not returned directly as a simple
//    polynomial.  Instead, it is represented in terms of a set of
//    orthogonal polynomials appopriate for the given data.  This makes
//    the computation more accurate, but means that the user can not
//    easily evaluate the computed polynomial.  Instead, the routine
//    LEAST_EVAL should be used to evaluate the least squares polynomial
//    at any point.  (However, the value of the least squares polynomial
//    at each of the data points is returned as part of this computation.)
//
//
//    A discrete unweighted inner product is used, so that
//
//      ( F(X), G(X) ) = sum ( 1 <= I <= NTAB ) F(XTAB(I)) * G(XTAB(I)).
//
//    The least squares polynomial is determined using a set of
//    orthogonal polynomials PHI.  These polynomials can be defined
//    recursively by:
//
//      PHI(0)(X) = 1
//      PHI(1)(X) = X - B(1)
//      PHI(I)(X) = ( X - B(I) ) * PHI(I-1)(X) - D(I) * PHI(I-2)(X)
//
//    The array B(1:NDEG) contains the values
//
//      B(I) = ( X*PHI(I-1), PHI(I-1) ) / ( PHI(I-1), PHI(I-1) )
//
//    The array D(2:NDEG) contains the values
//
//      D(I) = ( PHI(I-1), PHI(I-1) ) / ( PHI(I-2), PHI(I-2) )
//
//    Using this basis, the least squares polynomial can be represented as
//
//      P(X)(I) = sum ( 0 <= I <= NDEG ) C(I) * PHI(I)(X)
//
//    The array C(0:NDEG) contains the values
//
//      C(I) = ( YTAB(I), PHI(I) ) / ( PHI(I), PHI(I) )
//
//  Modified:
//
//    15 May 2004
//
//  Reference:
//
//    Gisela Engeln-Muellges and Frank Uhlig,
//    Numerical Algorithms with C, pages 191-193.
//    Springer, 1996.
//
//  Parameters:
//
//    Input, int NTAB, the number of data points.
//
//    Input, double XTAB[NTAB], the X data.  The values in XTAB
//    should be distinct, and in increasing order.
//
//    Input, double YTAB[NTAB], the Y data values corresponding
//    to the X data in XTAB.
//
//    Input, int NDEG, the degree of the polynomial which the
//    program is to use.  NDEG must be at least 1, and less than or
//    equal to NTAB-1.
//
//    Output, double PTAB[NTAB], the value of the least squares polynomial
//    at the points XTAB(1:NTAB).
//
//    Output, double B[1:NDEG], C[0:NDEG], D[2:NDEG], arrays containing
//    data about the polynomial.
//
//    Output, double *EPS, the root-mean-square discrepancy of the
//    polynomial fit.
//
//    Output, int *IERROR, error flag.
//    zero, no error occurred;
//    nonzero, an error occurred, and the polynomial could not be computed.
//
{
    int B_OFFSET = -1;
    int D_OFFSET = -2;
    int i;
    int i0l1;
    int i1l1;
    int it;
    int k;
    int mdeg;
    double rn0;
    double rn1;
    double s;
    double sum2;
    double y_sum;
    double *ztab;
    //
    *ierror = 0;
    ztab = new double[2*ntab];
    //
    //  Check NDEG.
    //
    if ( ndeg < 1 ) {
        *ierror = 1;
        cout << "\n";
        cout << "LEAST_SET - Fatal error!\n";
        cout << "  NDEG < 1.\n";
        exit ( 1 );
    }

    if ( ntab <= ndeg ) {
        *ierror = 1;
        cout << "\n";
        cout << "LEAST_SET - Fatal error!\n";
        cout << "  NTAB <= NDEG.\n";
        exit ( 1 );
    }
    //
    //  Check that the abscissas are strictly increasing.
    //
    for ( i = 1; i <= ntab-1; i++ ) {
        if ( xtab[i] <= xtab[i-1] ) {
            *ierror = 1;
            cout << "\n";
            cout << "LEAST_SET - Fatal error!\n";
            cout << "  XTAB must be strictly increasing, but\n";
            cout << "  XTAB(" << i-1 << ") = " << xtab[i-1] << "\n";
            cout << "  XTAB(" << i   << ") = " << xtab[i]   << "\n";
            exit ( 1 );
        }
    }

    i0l1 = 0;
    i1l1 = ntab;
    //
    //  The polynomial is of degree at least zero.
    //
    y_sum = 0.0E+00;
    for ( i = 0; i < ntab; i++ ) {
        y_sum = y_sum + ytab[i];
    }

    rn0 = ntab;
    c[0] = y_sum / ( double ) ( ntab );

    for ( i = 0; i < ntab; i++ ) {
        ptab[i] = y_sum / ( double ) ( ntab );
    }

    if ( ndeg == 0 ) {
        *eps = 0.0E+00;
        for ( i = 0; i < ntab; i++ ) {
            *eps = *eps + pow ( ( y_sum / ( double ) ( ntab ) - ytab[i] ), 2 );
        }

        *eps = sqrt ( *eps / ( double ) ( ntab ) );
        delete [] ztab;
        return;
    }
    //
    //  The polynomial is of degree at least 1.
    //
    ztab[0] = 0.0E+00;
    for ( i = 0; i < ntab; i++ ) {
        ztab[0] = ztab[0] + xtab[i];
    }

    b[1+B_OFFSET] = ztab[0] / ( double ) ( ntab );

    s = 0.0E+00;
    sum2 = 0.0E+00;
    for ( i = 0; i < ntab; i++ ) {
        ztab[i1l1+i] = xtab[i] - b[1+B_OFFSET];
        s = s + ztab[i1l1+i] * ztab[i1l1+i];
        sum2 = sum2 + ztab[i1l1+i] * ( ytab[i] - ptab[i] );
    }

    rn1 = s;
    c[1] = sum2 / s;

    for ( i = 0; i < ntab; i++ ) {
        ptab[i] = ptab[i] + c[1] * ztab[i1l1+i];
    }


    if ( ndeg == 1 ) {
        *eps = 0.0E+00;
        for ( i = 0; i < ntab; i++ ) {
            *eps = *eps + pow ( ( ptab[i] - ytab[i] ), 2 );
        }

        *eps = sqrt ( *eps / ( double ) ( ntab ) );
        delete [] ztab;
        return;
    }

    for ( i = 0; i < ntab; i++ ) {
        ztab[i] = 1.0E+00;
    }

    mdeg = 2;
    k = 2;

    for ( ; ; ) {
        d[k+D_OFFSET] = rn1 / rn0;

        sum2 = 0.0E+00;
        for ( i = 0; i < ntab; i++ ) {
            sum2 = sum2 + xtab[i] * ztab[i1l1+i] * ztab[i1l1+i];
        }

        b[k+B_OFFSET] = sum2 / rn1;

        s = 0.0E+00;
        sum2 = 0.0E+00;

        for ( i = 0; i < ntab; i++ ) {
            ztab[i0l1+i] = ( xtab[i] - b[k+B_OFFSET] ) * ztab[i1l1+i]
                           - d[k+D_OFFSET] * ztab[i0l1+i];
            s = s + ztab[i0l1+i] * ztab[i0l1+i];
            sum2 = sum2 + ztab[i0l1+i] * ( ytab[i] - ptab[i] );
        }

        rn0 = rn1;
        rn1 = s;

        c[k] = sum2 / rn1;

        it = i0l1;
        i0l1 = i1l1;
        i1l1 = it;

        for ( i = 0; i < ntab; i++ ) {
            ptab[i] = ptab[i] + c[k] * ztab[i1l1+i];
        }

        if ( ndeg <= mdeg ) {
            break;
        }

        mdeg = mdeg + 1;
        k = k + 1;

    }
    //
    //  Compute the RMS error.
    //
    *eps = 0.0E+00;
    for ( i = 0; i < ntab; i++ ) {
        *eps = *eps + pow ( ( ptab[i] - ytab[i] ), 2 );
    }

    *eps = sqrt ( *eps / ( double ) ( ntab ) );
    delete [] ztab;

    return;
}
//******************************************************************************

double least_val ( double x, int ndeg, double b[], double c[], double d[] )

//******************************************************************************
//
//  Purpose:
//
//    LEAST_VAL evaluates a least squares polynomial defined by LEAST_SET.
//
//  Modified:
//
//    17 December 2004
//
//  Reference:
//
//    Gisela Engeln-Muellges and Frank Uhlig,
//    Numerical Algorithms with C, pages 191-193.
//    Springer, 1996.
//
//  Parameters:
//
//    Input, double X, the point at which the polynomial is to be evaluated.
//
//    Input, int NDEG, the degree of the polynomial fit used.
//    This is the value of NDEG as returned from LEAST_SET.
//
//    Input, double B[1:NDEG], C[0:NDEG], D[2:NDEG], arrays defined by
//    LEAST_SET, and needed to evaluate the polynomial.
//
//    Output, double LEAST_VAL, the value of the polynomial at X.
//
{
    int B_OFFSET = -1;
    int D_OFFSET = -2;
    int k;
    double sk = 0.0;
    double skp1;
    double skp2;
    double value;
    //
    if ( ndeg <= 0 ) {
        value = c[0];
    } else if ( ndeg == 1 ) {
        value = c[0] + c[1] * ( x - b[1+B_OFFSET] );
    } else {
        skp2 = c[ndeg];
        skp1 = c[ndeg-1] + ( x - b[ndeg+B_OFFSET] ) * skp2;

        for ( k = ndeg-2; 0 <= k; k-- ) {
            sk = c[k] + ( x - b[k+1+B_OFFSET] ) * skp1 - d[k+2+D_OFFSET] * skp2;
            skp2 = skp1;
            skp1 = sk;
        }
        value = sk;
    }

    return value;
}
//**********************************************************************

void parabola_val2 ( int ndim, int ndata, double tdata[], double ydata[],
                     int left, double tval, double *yval )

//**********************************************************************
//
//  Purpose:
//
//    PARABOLA_VAL2 evaluates a parabolic function through 3 points in a table.
//
//  Discussion:
//
//    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
//    It constructs the parabolic interpolant through the data in
//    3 consecutive entries of a table and evaluates this interpolant
//    at a given abscissa value.
//
//  Modified:
//
//    11 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer NDIM, the dimension of a single data point.
//    NDIM must be at least 1.
//
//    Input, int NDATA, the number of data points.
//    NDATA must be at least 3.
//
//    Input, double TDATA[NDATA], the abscissas of the data points.  The
//    values in TDATA must be in strictly ascending order.
//
//    Input, double YDATA[NDIM*NDATA], the data points corresponding to
//    the abscissas.
//
//    Input, int LEFT, the location of the first of the three
//    consecutive data points through which the parabolic interpolant
//    must pass.  0 <= LEFT <= NDATA - 3.
//
//    Input, double TVAL, the value of T at which the parabolic interpolant
//    is to be evaluated.  Normally, TDATA[0] <= TVAL <= T[NDATA-1], and
//    the data will be interpolated.  For TVAL outside this range,
//    extrapolation will be used.
//
//    Output, double YVAL[NDIM], the value of the parabolic interpolant
//    at TVAL.
//
{
    double dif1;
    double dif2;
    int i;
    double t1;
    double t2;
    double t3;
    double y1;
    double y2;
    double y3;
    //
    //  Check.
    //
    if ( left < 1 ) {
        cout << "\n";
        cout << "PARABOLA_VAL2 - Fatal error!\n";
        cout << "  LEFT < 0.\n";
        exit ( 1 );
    }

    if ( ndata-2 < left ) {
        cout << "\n";
        cout << "PARABOLA_VAL2 - Fatal error!\n";
        cout << " NDATA-2 < LEFT.\n";
        exit ( 1 );
    }

    if ( ndim < 1 ) {
        cout << "\n";
        cout << "PARABOLA_VAL2 - Fatal error!\n";
        cout << " NDIM < 1.\n";
        exit ( 1 );
    }
    //
    //  Copy out the three abscissas.
    //
    t1 = tdata[left-1];
    t2 = tdata[left];
    t3 = tdata[left+1];

    if ( t2 <= t1 || t3 <= t2 ) {
        cout << "\n" ;
        cout << "PARABOLA_VAL2 - Fatal error!\n";
        cout << "  T2 <= T1 or T3 <= T2.\n";
        cout << "  T1 = " << t1 << "\n";
        cout << "  T2 = " << t2 << "\n";
        cout << "  T3 = " << t3 << "\n";
        exit ( 1 );
    }
    //
    //  Construct and evaluate a parabolic interpolant for the data.
    //
    for ( i = 0; i < ndim; i++ ) {
        y1 = ydata[i+(left-1)*ndim];
        y2 = ydata[i+(left  )*ndim];
        y3 = ydata[i+(left+1)*ndim];

        dif1 = ( y2 - y1 ) / ( t2 - t1 );
        dif2 =
            ( ( y3 - y1 ) / ( t3 - t1 )
              - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 );

        yval[i] = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 );
    }

    return;
}
//******************************************************************************

int s_len_trim ( char* s )

//******************************************************************************
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
    int n;
    char* t;

    n = strlen ( s );
    t = s + strlen ( s ) - 1;

    while ( 0 < n ) {
        if ( *t != ' ' ) {
            return n;
        }
        t--;
        n--;
    }

    return n;
}
//******************************************************************************

double spline_b_val ( int ndata, double tdata[], double ydata[], double tval )

//******************************************************************************
//
//  Purpose:
//
//    SPLINE_B_VAL evaluates a cubic B spline approximant.
//
//  Discussion:
//
//    The cubic B spline will approximate the data, but is not
//    designed to interpolate it.
//
//    In effect, two "phantom" data values are appended to the data,
//    so that the spline will interpolate the first and last data values.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl de Boor,
//    A Practical Guide to Splines,
//    Springer Verlag, 1978.
//
//  Parameters:
//
//    Input, int NDATA, the number of data values.
//
//    Input, double TDATA[NDATA], the abscissas of the data.
//
//    Input, double YDATA[NDATA], the data values.
//
//    Input, double TVAL, a point at which the spline is to be evaluated.
//
//    Output, double SPLINE_B_VAL, the value of the function at TVAL.
//
{
    double bval;
    int left;
    int right;
    double u;
    double yval;
    //
    //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    //
    dvec_bracket ( ndata, tdata, tval, &left, &right );
    //
    //  Evaluate the 5 nonzero B spline basis functions in the interval,
    //  weighted by their corresponding data values.
    //
    u = ( tval - tdata[left-1] ) / ( tdata[right-1] - tdata[left-1] );
    yval = 0.0E+00;
    //
    //  B function associated with node LEFT - 1, (or "phantom node"),
    //  evaluated in its 4th interval.
    //
    bval = ( ( (     - 1.0E+00
                     * u + 3.0E+00 )
               * u - 3.0E+00 )
             * u + 1.0E+00 ) / 6.0E+00;

    if ( 0 < left-1 ) {
        yval = yval + ydata[left-2] * bval;
    } else {
        yval = yval + ( 2.0E+00 * ydata[0] - ydata[1] ) * bval;
    }
    //
    //  B function associated with node LEFT,
    //  evaluated in its third interval.
    //
    bval = ( ( (       3.0E+00
                       * u - 6.0E+00 )
               * u + 0.0E+00 )
             * u + 4.0E+00 ) / 6.0E+00;

    yval = yval + ydata[left-1] * bval;
    //
    //  B function associated with node RIGHT,
    //  evaluated in its second interval.
    //
    bval = ( ( (     - 3.0E+00
                     * u + 3.0E+00 )
               * u + 3.0E+00 )
             * u + 1.0E+00 ) / 6.0E+00;

    yval = yval + ydata[right-1] * bval;
    //
    //  B function associated with node RIGHT+1, (or "phantom node"),
    //  evaluated in its first interval.
    //
    bval = pow ( u, 3 ) / 6.0E+00;

    if ( right+1 <= ndata ) {
        yval = yval + ydata[right] * bval;
    } else {
        yval = yval + ( 2.0E+00 * ydata[ndata-1] - ydata[ndata-2] ) * bval;
    }

    return yval;
}
//******************************************************************************

double spline_beta_val ( double beta1, double beta2, int ndata, double tdata[],
                         double ydata[], double tval )

//******************************************************************************
//
//  Purpose:
//
//    SPLINE_BETA_VAL evaluates a cubic beta spline approximant.
//
//  Discussion:
//
//    The cubic beta spline will approximate the data, but is not
//    designed to interpolate it.
//
//    If BETA1 = 1 and BETA2 = 0, the cubic beta spline will be the
//    same as the cubic B spline approximant.
//
//    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
//    a linear spline.
//
//    In effect, two "phantom" data values are appended to the data,
//    so that the spline will interpolate the first and last data values.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double BETA1, the skew or bias parameter.
//    BETA1 = 1 for no skew or bias.
//
//    Input, double BETA2, the tension parameter.
//    BETA2 = 0 for no tension.
//
//    Input, int NDATA, the number of data values.
//
//    Input, double TDATA[NDATA], the abscissas of the data.
//
//    Input, double YDATA[NDATA], the data values.
//
//    Input, double TVAL, a point at which the spline is to be evaluated.
//
//    Output, double SPLINE_BETA_VAL, the value of the function at TVAL.
//
{
    double a;
    double b;
    double bval;
    double c;
    double d;
    double delta;
    int left;
    int right;
    double u;
    double yval;
    //
    //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    //
    dvec_bracket ( ndata, tdata, tval, &left, &right );
    //
    //  Evaluate the 5 nonzero beta spline basis functions in the interval,
    //  weighted by their corresponding data values.
    //
    u = ( tval - tdata[left-1] ) / ( tdata[right-1] - tdata[left-1] );

    delta = ( ( 2.0E+00
                * beta1 + 4.0E+00 )
              * beta1 + 4.0E+00 )
            * beta1 + 2.0E+00 + beta2;

    yval = 0.0E+00;
    //
    //  Beta function associated with node LEFT - 1, (or "phantom node"),
    //  evaluated in its 4th interval.
    //
    bval = 2.0E+00 * pow ( ( beta1 * ( 1.0E+00 - u ) ), 3 )  / delta;

    if ( 0 < left-1 ) {
        yval = yval + ydata[left-2] * bval;
    } else {
        yval = yval + ( 2.0E+00 * ydata[0] - ydata[1] ) * bval;
    }
    //
    //  Beta function associated with node LEFT,
    //  evaluated in its third interval.
    //
    a = beta2 + ( 4.0E+00 + 4.0E+00 * beta1 ) * beta1;

    b = - 6.0E+00 * beta1 * ( 1.0E+00 - beta1 ) * ( 1.0E+00 + beta1 );

    c = ( (     - 6.0E+00
                * beta1 - 6.0E+00 )
          * beta1 + 0.0E+00 )
        * beta1 - 3.0E+00 * beta2;

    d = ( (     + 2.0E+00
                * beta1 + 2.0E+00 )
          * beta1 + 2.0E+00 )
        * beta1 + 2.0E+00 * beta2;

    bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta;

    yval = yval + ydata[left-1] * bval;
    //
    //  Beta function associated with node RIGHT,
    //  evaluated in its second interval.
    //
    a = 2.0E+00;

    b = 6.0E+00 * beta1;

    c = 3.0E+00 * beta2 + 6.0E+00 * beta1 * beta1;

    d = - 2.0E+00 * ( 1.0E+00 + beta2 + beta1 + beta1 * beta1 );

    bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta;

    yval = yval + ydata[right-1] * bval;
    //
    //  Beta function associated with node RIGHT+1, (or "phantom node"),
    //  evaluated in its first interval.
    //
    bval = 2.0E+00 * pow ( u, 3 ) / delta;

    if ( right+1 <= ndata ) {
        yval = yval + ydata[right] * bval;
    } else {
        yval = yval + ( 2.0E+00 * ydata[ndata-1] - ydata[ndata-2] ) * bval;
    }

    return yval;
}
//******************************************************************************

double spline_constant_val ( int ndata, double tdata[], double ydata[], double tval )

//******************************************************************************
//
//  Purpose:
//
//    SPLINE_CONSTANT_VAL evaluates a piecewise constant spline at a point.
//
//  Discussion:
//
//    NDATA-1 points TDATA define NDATA intervals, with the first
//    and last being semi-infinite.
//
//    The value of the spline anywhere in interval I is YDATA(I).
//
//  Modified:
//
//    10 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDATA, the number of data points defining the spline.
//
//    Input, double TDATA[NDATA-1], the breakpoints.  The values of TDATA should
//    be distinct and increasing.
//
//    Input, double YDATA[NDATA], the values of the spline in the intervals
//    defined by the breakpoints.
//
//    Input, double TVAL, the point at which the spline is to be evaluated.
//
//    Output, double *SPLINE_CONSTANT_VAL, the value of the spline at TVAL.
//
{
    int i;

    for ( i = 0; i < ndata-1; i++ ) {
        if ( tval <= tdata[i] ) {
            return ydata[i];
        }
    }

    return ydata[ndata-1];
}
//**********************************************************************

double *spline_cubic_set ( int n, double t[], double y[], int ibcbeg, double ybcbeg,
                           int ibcend, double ybcend )

//**********************************************************************
//
//  Purpose:
//
//    SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
//
//  Discussion:
//
//    For data interpolation, the user must call SPLINE_SET to determine
//    the second derivative data, passing in the data to be interpolated,
//    and the desired boundary conditions.
//
//    The data to be interpolated, plus the SPLINE_SET output, defines
//    the spline.  The user may then call SPLINE_VAL to evaluate the
//    spline at any point.
//
//    The cubic spline is a piecewise cubic polynomial.  The intervals
//    are determined by the "knots" or abscissas of the data to be
//    interpolated.  The cubic spline has continous first and second
//    derivatives over the entire interval of interpolation.
//
//    For any point T in the interval T(IVAL), T(IVAL+1), the form of
//    the spline is
//
//      SPL(T) = A(IVAL)
//             + B(IVAL) * ( T - T(IVAL) )
//             + C(IVAL) * ( T - T(IVAL) )**2
//             + D(IVAL) * ( T - T(IVAL) )**3
//
//    If we assume that we know the values Y(*) and YPP(*), which represent
//    the values and second derivatives of the spline at each knot, then
//    the coefficients can be computed as:
//
//      A(IVAL) = Y(IVAL)
//      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
//        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
//      C(IVAL) = YPP(IVAL) / 2
//      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
//
//    Since the first derivative of the spline is
//
//      SPL'(T) =     B(IVAL)
//              + 2 * C(IVAL) * ( T - T(IVAL) )
//              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
//
//    the requirement that the first derivative be continuous at interior
//    knot I results in a total of N-2 equations, of the form:
//
//      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
//      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
//
//    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
//
//      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
//      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
//      + YPP(IVAL-1) * H(IVAL-1)
//      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
//      =
//      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
//      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
//
//    or
//
//      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
//      + YPP(IVAL) * H(IVAL)
//      =
//      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
//      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
//
//    Boundary conditions must be applied at the first and last knots.
//    The resulting tridiagonal system can be solved for the YPP values.
//
//  Modified:
//
//    06 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data points.  N must be at least 2.
//    In the special case where N = 2 and IBCBEG = IBCEND = 0, the
//    spline will actually be linear.
//
//    Input, double T[N], the knot values, that is, the points were data is
//    specified.  The knot values should be distinct, and increasing.
//
//    Input, double Y[N], the data values to be interpolated.
//
//    Input, int IBCBEG, left boundary condition flag:
//      0: the cubic spline should be a quadratic over the first interval;
//      1: the first derivative at the left endpoint should be YBCBEG;
//      2: the second derivative at the left endpoint should be YBCBEG.
//
//    Input, double YBCBEG, the values to be used in the boundary
//    conditions if IBCBEG is equal to 1 or 2.
//
//    Input, int IBCEND, right boundary condition flag:
//      0: the cubic spline should be a quadratic over the last interval;
//      1: the first derivative at the right endpoint should be YBCEND;
//      2: the second derivative at the right endpoint should be YBCEND.
//
//    Input, double YBCEND, the values to be used in the boundary
//    conditions if IBCEND is equal to 1 or 2.
//
//    Output, double SPLINE_CUBIC_SET[N], the second derivatives of the cubic spline.
//
{
    double *a;
    double *b;
    int i;
    double *ypp;
    //
    //  Check.
    //
    if ( n <= 1 ) {
        cout << "\n";
        cout << "SPLINE_CUBIC_SET - Fatal error!\n";
        cout << "  The number of data points N must be at least 2.\n";
        cout << "  The input value is " << n << ".\n";
        return NULL;
    }

    for ( i = 0; i < n - 1; i++ ) {
        if ( t[i+1] <= t[i] ) {
            cout << "\n";
            cout << "SPLINE_CUBIC_SET - Fatal error!\n";
            cout << "  The knots must be strictly increasing, but\n";
            cout << "  T(" << i   << ") = " << t[i]   << "\n";
            cout << "  T(" << i+1 << ") = " << t[i+1] << "\n";
            return NULL;
        }
    }
    a = new double[3*n];
    b = new double[n];
    //
    //  Set up the first equation.
    //
    if ( ibcbeg == 0 ) {
        b[0] = 0.0E+00;
        a[1+0*3] = 1.0E+00;
        a[0+1*3] = -1.0E+00;
    } else if ( ibcbeg == 1 ) {
        b[0] = ( y[1] - y[0] ) / ( t[1] - t[0] ) - ybcbeg;
        a[1+0*3] = ( t[1] - t[0] ) / 3.0E+00;
        a[0+1*3] = ( t[1] - t[0] ) / 6.0E+00;
    } else if ( ibcbeg == 2 ) {
        b[0] = ybcbeg;
        a[1+0*3] = 1.0E+00;
        a[0+1*3] = 0.0E+00;
    } else {
        cout << "\n";
        cout << "SPLINE_CUBIC_SET - Fatal error!\n";
        cout << "  IBCBEG must be 0, 1 or 2.\n";
        cout << "  The input value is " << ibcbeg << ".\n";
        delete [] a;
        delete [] b;
        return NULL;
    }
    //
    //  Set up the intermediate equations.
    //
    for ( i = 1; i < n-1; i++ ) {
        b[i] = ( y[i+1] - y[i] ) / ( t[i+1] - t[i] )
               - ( y[i] - y[i-1] ) / ( t[i] - t[i-1] );
        a[2+(i-1)*3] = ( t[i] - t[i-1] ) / 6.0E+00;
        a[1+ i   *3] = ( t[i+1] - t[i-1] ) / 3.0E+00;
        a[0+(i+1)*3] = ( t[i+1] - t[i] ) / 6.0E+00;
    }
    //
    //  Set up the last equation.
    //
    if ( ibcend == 0 ) {
        b[n-1] = 0.0E+00;
        a[2+(n-2)*3] = -1.0E+00;
        a[1+(n-1)*3] = 1.0E+00;
    } else if ( ibcend == 1 ) {
        b[n-1] = ybcend - ( y[n-1] - y[n-2] ) / ( t[n-1] - t[n-2] );
        a[2+(n-2)*3] = ( t[n-1] - t[n-2] ) / 6.0E+00;
        a[1+(n-1)*3] = ( t[n-1] - t[n-2] ) / 3.0E+00;
    } else if ( ibcend == 2 ) {
        b[n-1] = ybcend;
        a[2+(n-2)*3] = 0.0E+00;
        a[1+(n-1)*3] = 1.0E+00;
    } else {
        cout << "\n";
        cout << "SPLINE_CUBIC_SET - Fatal error!\n";
        cout << "  IBCEND must be 0, 1 or 2.\n";
        cout << "  The input value is " << ibcend << ".\n";
        delete [] a;
        delete [] b;
        return NULL;
    }
    //
    //  Solve the linear system.
    //
    if ( n == 2 && ibcbeg == 0 && ibcend == 0 ) {
        ypp = new double[2];

        ypp[0] = 0.0E+00;
        ypp[1] = 0.0E+00;
    } else {
        ypp = d3_np_fs ( n, a, b );

        if ( !ypp ) {
            cout << "\n";
            cout << "SPLINE_CUBIC_SET - Fatal error!\n";
            cout << "  The linear system could not be solved.\n";
            delete [] a;
            delete [] b;
            return NULL;
        }

    }

    delete [] a;
    delete [] b;
    return ypp;
}
//**********************************************************************

double spline_cubic_val ( int n, double t[], double tval, double y[],
                          double ypp[], double *ypval, double *yppval )

//**********************************************************************
//
//  Purpose:
//
//    SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
//
//  Discussion:
//
//    SPLINE_CUBIC_SET must have already been called to define the values of YPP.
//
//    For any point T in the interval T(IVAL), T(IVAL+1), the form of
//    the spline is
//
//      SPL(T) = A
//             + B * ( T - T(IVAL) )
//             + C * ( T - T(IVAL) )**2
//             + D * ( T - T(IVAL) )**3
//
//    Here:
//      A = Y(IVAL)
//      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
//        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
//      C = YPP(IVAL) / 2
//      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
//
//  Modified:
//
//    04 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int n, the number of knots.
//
//    Input, double Y[N], the data values at the knots.
//
//    Input, double T[N], the knot values.
//
//    Input, double TVAL, a point, typically between T[0] and T[N-1], at
//    which the spline is to be evalulated.  If TVAL lies outside
//    this range, extrapolation is used.
//
//    Input, double Y[N], the data values at the knots.
//
//    Input, double YPP[N], the second derivatives of the spline at
//    the knots.
//
//    Output, double *YPVAL, the derivative of the spline at TVAL.
//
//    Output, double *YPPVAL, the second derivative of the spline at TVAL.
//
//    Output, double SPLINE_VAL, the value of the spline at TVAL.
//
{
    double dt;
    double h;
    int i;
    int ival;
    double yval;
    //
    //  Determine the interval [ T(I), T(I+1) ] that contains TVAL.
    //  Values below T[0] or above T[N-1] use extrapolation.
    //
    ival = n - 2;

    for ( i = 0; i < n-1; i++ ) {
        if ( tval < t[i+1] ) {
            ival = i;
            break;
        }
    }
    //
    //  In the interval I, the polynomial is in terms of a normalized
    //  coordinate between 0 and 1.
    //
    dt = tval - t[ival];
    h = t[ival+1] - t[ival];

    yval = y[ival]
           + dt * ( ( y[ival+1] - y[ival] ) / h
                    - ( ypp[ival+1] / 6.0E+00 + ypp[ival] / 3.0E+00 ) * h
                    + dt * ( 0.5E+00 * ypp[ival]
                             + dt * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0E+00 * h ) ) ) );

    *ypval = ( y[ival+1] - y[ival] ) / h
             - ( ypp[ival+1] / 6.0E+00 + ypp[ival] / 3.0E+00 ) * h
             + dt * ( ypp[ival]
                      + dt * ( 0.5E+00 * ( ypp[ival+1] - ypp[ival] ) / h ) );

    *yppval = ypp[ival] + dt * ( ypp[ival+1] - ypp[ival] ) / h;

    return yval;
}
//**********************************************************************

void spline_cubic_val2 ( int n, double t[], double tval, int *left, double y[],
                         double ypp[], double *yval, double *ypval, double *yppval )

//**********************************************************************
//
//  Purpose:
//
//    SPLINE_CUBIC_VAL2 evaluates a piecewise cubic spline at a point.
//
//  Discussion:
//
//    This routine is a modification of SPLINE_CUBIC_VAL; it allows the
//    user to speed up the code by suggesting the appropriate T interval
//    to search first.
//
//    SPLINE_CUBIC_SET must have already been called to define the
//    values of YPP.
//
//    In the LEFT interval, let RIGHT = LEFT+1.  The form of the spline is
//
//    SPL(T) =
//      A
//    + B * ( T - T[LEFT] )
//    + C * ( T - T[LEFT] )**2
//    + D * ( T - T[LEFT] )**3
//
//    Here:
//      A = Y[LEFT]
//      B = ( Y[RIGHT] - Y[LEFT] ) / ( T[RIGHT] - T[LEFT] )
//        - ( YPP[RIGHT] + 2 * YPP[LEFT] ) * ( T[RIGHT] - T[LEFT] ) / 6
//      C = YPP[LEFT] / 2
//      D = ( YPP[RIGHT] - YPP[LEFT] ) / ( 6 * ( T[RIGHT] - T[LEFT] ) )
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of knots.
//
//    Input, double T[N], the knot values.
//
//    Input, double TVAL, a point, typically between T[0] and T[N-1], at
//    which the spline is to be evalulated.  If TVAL lies outside
//    this range, extrapolation is used.
//
//    Input/output, int *LEFT, the suggested T interval to search.
//    LEFT should be between 1 and N-1.  If LEFT is not in this range,
//    then its value will be ignored.  On output, LEFT is set to the
//    actual interval in which TVAL lies.
//
//    Input, double Y[N], the data values at the knots.
//
//    Input, double YPP[N], the second derivatives of the spline at
//    the knots.
//
//    Output, double *YVAL, *YPVAL, *YPPVAL, the value of the spline, and
//    its first two derivatives at TVAL.
//
{
    double dt;
    double h;
    int right;
    //
    //  Determine the interval [T[LEFT], T[RIGHT]] that contains TVAL.
    //
    //  What you want from DVEC_BRACKET3 is that TVAL is to be computed
    //  by the data in interval [T[LEFT-1], T[RIGHT-1]].
    //
    dvec_bracket3 ( n, t, tval, left );
    //
    // In the interval LEFT, the polynomial is in terms of a normalized
    // coordinate  ( DT / H ) between 0 and 1.
    //
    right = *left + 1;

    dt = tval - t[*left-1];
    h = t[right-1] - t[*left-1];

    *yval = y[*left-1]
            + dt * ( ( y[right-1] - y[*left-1] ) / h
                     - ( ypp[right-1] / 6.0E+00 + ypp[*left-1] / 3.0E+00 ) * h
                     + dt * ( 0.5E+00 * ypp[*left-1]
                              + dt * ( ( ypp[right-1] - ypp[*left-1] ) / ( 6.0E+00 * h ) ) ) );

    *ypval = ( y[right-1] - y[*left-1] ) / h
             - ( ypp[right-1] / 6.0E+00 + ypp[*left-1] / 3.0E+00 ) * h
             + dt * ( ypp[*left-1]
                      + dt * ( 0.5E+00 * ( ypp[right-1] - ypp[*left-1] ) / h ) );

    *yppval = ypp[*left-1] + dt * ( ypp[right-1] - ypp[*left-1] ) / h;

    return;
}
//************************************************************************

double *spline_hermite_set ( int ndata, double tdata[], double ydata[],
                             double ypdata[] )

//************************************************************************
//
//  Purpose:
//
//    SPLINE_HERMITE_SET sets up a piecewise cubic Hermite interpolant.
//
//  Discussion:
//
//    Once the array C is computed, then in the interval
//    (TDATA(I), TDATA(I+1)), the interpolating Hermite polynomial
//    is given by
//
//      SVAL(TVAL) =                 C(1,I)
//         + ( TVAL - TDATA(I) ) * ( C(2,I)
//         + ( TVAL - TDATA(I) ) * ( C(3,I)
//         + ( TVAL - TDATA(I) ) *   C(4,I) ) )
//
//  Modified:
//
//    11 February 2004
//
//  Reference:
//
//    Conte and de Boor,
//    Algorithm CALCCF,
//    Elementary Numerical Analysis,
//    1973, page 235.
//
//  Parameters:
//
//    Input, int NDATA, the number of data points.
//    NDATA must be at least 2.
//
//    Input, double TDATA[NDATA], the abscissas of the data points.
//    The entries of TDATA are assumed to be strictly increasing.
//
//    Input, double Y[NDATA], YP[NDATA], the value of the
//    function and its derivative at TDATA(1:NDATA).
//
//    Output, double SPLINE_HERMITE_SET[4*NDATA], the coefficients of
//    the Hermite polynomial.  We will refer to this array as "C".
//    C(1,1:NDATA) = Y(1:NDATA) and C(2,1:NDATA) = YP(1:NDATA).
//    C(3,1:NDATA-1) and C(4,1:NDATA-1) are the quadratic and cubic
//    coefficients.
//
{
    double *c;
    double divdif1;
    double divdif3;
    double dt;
    int i;
    int j;

    c = new double[4*ndata];

    for ( j = 0; j < ndata; j++ ) {
        c[0+j*4] = ydata[j];
    }

    for ( j = 0; j < ndata; j++ ) {
        c[1+j*4] = ypdata[j];
    }

    for ( i = 1; i <= ndata-1; i++ ) {
        dt = tdata[i] - tdata[i-1];
        divdif1 = ( c[0+i*4] - c[0+(i-1)*4] ) / dt;
        divdif3 = c[1+(i-1)*4] + c[1+i*4] - 2.0E+00 * divdif1;
        c[2+(i-1)*4] = ( divdif1 - c[1+(i-1)*4] - divdif3 ) / dt;
        c[3+(i-1)*4] = divdif3 / ( dt * dt );
    }

    c[2+(ndata-1)*4] = 0.0E+00;
    c[3+(ndata-1)*4] = 0.0E+00;

    return c;
}
//************************************************************************

void spline_hermite_val ( int ndata, double tdata[], double c[], double tval,
                          double *sval, double *spval )

//************************************************************************
//
//  Purpose:
//
//    SPLINE_HERMITE_VAL evaluates a piecewise cubic Hermite interpolant.
//
//  Discussion:
//
//    SPLINE_HERMITE_SET must be called first, to set up the
//    spline data from the raw function and derivative data.
//
//    In the interval (TDATA(I), TDATA(I+1)), the interpolating
//    Hermite polynomial is given by
//
//      SVAL(TVAL) =                 C(1,I)
//         + ( TVAL - TDATA(I) ) * ( C(2,I)
//         + ( TVAL - TDATA(I) ) * ( C(3,I)
//         + ( TVAL - TDATA(I) ) *   C(4,I) ) )
//
//    and
//
//      SVAL'(TVAL) =                    C(2,I)
//         + ( TVAL - TDATA(I) ) * ( 2 * C(3,I)
//         + ( TVAL - TDATA(I) ) *   3 * C(4,I) )
//
//  Modified:
//
//    24 February 2004
//
//  Reference:
//
//    Conte and de Boor,
//    Algorithm PCUBIC,
//    Elementary Numerical Analysis,
//    1973, page 234.
//
//  Parameters:
//
//    Input, int NDATA, the number of data points.
//    NDATA is assumed to be at least 2.
//
//    Input, double TDATA[NDATA], the abscissas of the data points.
//    The entries of TDATA are assumed to be strictly increasing.
//
//    Input, double C[4*NDATA], the coefficient data computed by
//    SPLINE_HERMITE_SET.
//
//    Input, double TVAL, the point where the interpolant is to
//    be evaluated.
//
//    Output, double *SVAL, *SPVAL, the value of the interpolant
//    and its derivative at TVAL.
//
{
    double dt;
    int left;
    int right;
    //
    //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains
    //  or is nearest to TVAL.
    //
    dvec_bracket ( ndata, tdata, tval, &left, &right );
    //
    //  Evaluate the cubic polynomial.
    //
    dt = tval - tdata[left-1];

    *sval =        c[0+(left-1)*4]
                   + dt * ( c[1+(left-1)*4]
                            + dt * ( c[2+(left-1)*4]
                                     + dt *   c[3+(left-1)*4] ) );

    *spval =                 c[1+(left-1)*4]
                             + dt * ( 2.0E+00 * c[2+(left-1)*4]
                                      + dt *   3.0E+00 * c[3+(left-1)*4] );

    return;
}
//******************************************************************************

double spline_linear_int ( int ndata, double tdata[], double ydata[], double a, double b )

//******************************************************************************
//
//  Purpose:
//
//    SPLINE_LINEAR_INT evaluates the integral of a piecewise linear spline.
//
//  Modified:
//
//    25 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDATA, the number of data points defining the spline.
//
//    Input, double TDATA[NDATA], YDATA[NDATA], the values of the independent
//    and dependent variables at the data points.  The values of TDATA should
//    be distinct and increasing.
//
//    Input, double A, B, the interval over which the integral is desired.
//
//    Output, double SPLINE_LINEAR_INT, the value of the integral.
//
{
    double a_copy;
    int a_left;
    int a_right;
    double b_copy;
    int b_left;
    int b_right;
    int i_left;
    double int_val;
    double tval;
    double yp;
    double yval;

    int_val = 0.0E+00;

    if ( a == b ) {
        return int_val;
    }

    a_copy = d_min ( a, b );
    b_copy = d_max ( a, b );
    //
    //  Find the interval [ TDATA(A_LEFT), TDATA(A_RIGHT) ] that contains, or is
    //  nearest to, A.
    //
    dvec_bracket ( ndata, tdata, a_copy, &a_left, &a_right );
    //
    //  Find the interval [ TDATA(B_LEFT), TDATA(B_RIGHT) ] that contains, or is
    //  nearest to, B.
    //
    dvec_bracket ( ndata, tdata, b_copy, &b_left, &b_right );
    //
    //  If A and B are in the same interval...
    //
    if ( a_left == b_left ) {
        tval = ( a_copy + b_copy ) / 2.0E+00;

        yp = ( ydata[a_right-1] - ydata[a_left-1] ) /
             ( tdata[a_right-1] - tdata[a_left-1] );

        yval = ydata[a_left-1] + ( tval - tdata[a_left-1] ) * yp;

        int_val = yval * ( b_copy - a_copy );

        return int_val;
    }
    //
    //  Otherwise, integrate from:
    //
    //  A               to TDATA(A_RIGHT),
    //  TDATA(A_RIGHT)  to TDATA(A_RIGHT+1),...
    //  TDATA(B_LEFT-1) to TDATA(B_LEFT),
    //  TDATA(B_LEFT)   to B.
    //
    //  Use the fact that the integral of a linear function is the
    //  value of the function at the midpoint times the width of the interval.
    //
    tval = ( a_copy + tdata[a_right-1] ) / 2.0E+00;

    yp = ( ydata[a_right-1] - ydata[a_left-1] ) /
         ( tdata[a_right-1] - tdata[a_left-1] );

    yval = ydata[a_left-1] + ( tval - tdata[a_left-1] ) * yp;

    int_val = int_val + yval * ( tdata[a_right-1] - a_copy );

    for ( i_left = a_right; i_left <= b_left - 1; i_left++ ) {
        tval = ( tdata[i_left] + tdata[i_left-1] ) / 2.0E+00;

        yp = ( ydata[i_left-1] - ydata[i_left-2] ) /
             ( tdata[i_left-1] - tdata[i_left-2] );

        yval = ydata[i_left-2] + ( tval - tdata[i_left-2] ) * yp;

        int_val = int_val + yval * ( tdata[i_left-1] - tdata[i_left-2] );
    }

    tval = ( tdata[b_left-1] + b_copy ) / 2.0E+0;

    yp = ( ydata[b_right-1] - ydata[b_left-1] ) /
         ( tdata[b_right-1] - tdata[b_left-1] );

    yval = ydata[b_left-1] + ( tval - tdata[b_left-1] ) * yp;

    int_val = int_val + yval * ( b_copy - tdata[b_left-1] );

    if ( b < a ) {
        int_val = -int_val;
    }

    return int_val;
}
//******************************************************************************

void spline_linear_intset ( int n, double int_x[], double int_v[], double data_x[],
                            double data_y[] )

//******************************************************************************
//
//  Purpose:
//
//    SPLINE_LINEAR_INTSET sets a piecewise linear spline with given integral properties.
//
//  Discussion:
//
//    The user has in mind an interval, divided by N+1 points into
//    N intervals.  A linear spline is to be constructed,
//    with breakpoints at the centers of each interval, and extending
//    continuously to the left of the first and right of the last
//    breakpoints.  The constraint on the linear spline is that it is
//    required that it have a given integral value over each interval.
//
//    A tridiagonal linear system of equations is solved for the
//    values of the spline at the breakpoints.
//
//  Modified:
//
//    07 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of intervals.
//
//    Input, double INT_X[N+1], the points that define the intervals.
//    Interval I lies between INT_X(I) and INT_X(I+1).
//
//    Input, double INT_V[N], the desired value of the integral of the
//    linear spline over each interval.
//
//    Output, double DATA_X[N], DATA_Y[N], the values of the independent
//    and dependent variables at the data points.  The values of DATA_X are
//    the interval midpoints.  The values of DATA_Y are determined in such
//    a way that the exact integral of the linear spline over interval I
//    is equal to INT_V(I).
//
{
    double *a;
    double *b;
    double *c;
    int i;

    a = new double[3*n];
    b = new double[n];
    //
    //  Set up the easy stuff.
    //
    for ( i = 1; i <= n; i++ ) {
        data_x[i-1] = 0.5E+00 * ( int_x[i-1] + int_x[i] );
    }
    //
    //  Set up the coefficients of the linear system.
    //
    for ( i = 0; i < n-2; i++ ) {
        a[2+i*3] = 1.0E+00 - ( 0.5E+00 * ( data_x[i+1] + int_x[i+1] )
                               - data_x[i] ) / ( data_x[i+1] - data_x[i] );
    }
    a[2+(n-2)*3] = 0.0E+00;
    a[2+(n-1)*3] = 0.0E+00;

    a[1+0*3] = int_x[1] - int_x[0];

    for ( i = 1; i < n-1; i++ ) {
        a[1+i*3] = 1.0E+00 + ( 0.5E+00 * ( data_x[i] + int_x[i] )
                               - data_x[i-1] ) / ( data_x[i] - data_x[i-1] )
                   - ( 0.5E+00 * ( data_x[i] + int_x[i+1] ) - data_x[i] )
                   / ( data_x[i+1] - data_x[i] );
    }
    a[1+(n-1)*3] = int_x[n] - int_x[n-1];

    a[0+0*3] = 0.0E+00;
    a[0+1*3] = 0.0E+00;
    for ( i = 2; i < n; i++ ) {
        a[0+i*3] = ( 0.5E+00 * ( data_x[i-1] + int_x[i] )
                     - data_x[i-1] ) / ( data_x[i] - data_x[i-1] );
    }
    //
    //  Set up DATA_Y, which begins as the right hand side of the linear system.
    //
    b[0] = int_v[0];
    for ( i = 2; i <= n-1; i++ ) {
        b[i-1] = 2.0E+00 * int_v[i-1] / ( int_x[i] - int_x[i-1] );
    }
    b[n-1] = int_v[n-1];
    //
    //  Solve the linear system.
    //
    c = d3_np_fs ( n, a, b );

    for ( i = 0; i < n; i++ ) {
        data_y[i] = c[i];
    }

    delete [] a;
    delete [] b;
    delete [] c;

    return;
}
//**********************************************************************

void spline_linear_val ( int ndata, double tdata[], double ydata[],
                         double tval, double *yval, double *ypval )

//**********************************************************************
//
//  Purpose:
//
//    SPLINE_LINEAR_VAL evaluates a piecewise linear spline at a point.
//
//  Discussion:
//
//    Because of the extremely simple form of the linear spline,
//    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
//    evaluate the spline at any point.  No processing of the data
//    is required.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDATA, the number of data points defining the spline.
//
//    Input, double TDATA[NDATA], YDATA[NDATA], the values of the independent
//    and dependent variables at the data points.  The values of TDATA should
//    be distinct and increasing.
//
//    Input, double TVAL, the point at which the spline is to be evaluated.
//
//    Output, double *YVAL, *YPVAL, the value of the spline and its first
//    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
//    equal to TDATA(I) for some I.
//
{
    int left;
    int right;
    //
    //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
    //  nearest to, TVAL.
    //
    dvec_bracket ( ndata, tdata, tval, &left, &right );
    //
    //  Now evaluate the piecewise linear function.
    //
    *ypval = ( ydata[right-1] - ydata[left-1] )
             / ( tdata[right-1] - tdata[left-1] );

    *yval = ydata[left-1] +  ( tval - tdata[left-1] ) * (*ypval);

    return;
}
//******************************************************************************

double spline_overhauser_nonuni_val ( int ndata, double tdata[],
                                      double ydata[], double tval )

//******************************************************************************
//
//  Purpose:
//
//    SPLINE_OVERHAUSER_NONUNI_VAL evaluates the nonuniform Overhauser spline.
//
//  Discussion:
//
//    The nonuniformity refers to the fact that the abscissas values
//    need not be uniformly spaced.
//
//  Modified:
//
//    10 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDATA, the number of data points.
//    NDATA must be at least 3.
//
//    Input, double TDATA[NDATA], the abscissas of the data points.
//    The values of TDATA are assumed to be distinct and increasing.
//
//    Input, double YDATA[NDATA], the data values.
//
//    Input, double TVAL, the value where the spline is to
//    be evaluated.
//
//    Output, double SPLINE_OVERHAUSER_NONUNI_VAL, the value of the
//    spline at TVAL.
//
{
    double alpha;
    double beta;
    int left;
    double *mbasis;
    int right;
    double yval;
    //
    //  Check NDATA.
    //
    if ( ndata < 3 ) {
        cout << "\n";
        cout << "SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!\n";
        cout << "  NDATA < 3.\n";
        exit ( 1 );
    }
    //
    //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    //
    dvec_bracket ( ndata, tdata, tval, &left, &right );

    //  printf("nearest interval: %d - %d\n", left, right);

    //
    //  Evaluate the spline in the given interval.
    //
    if ( left == 1 ) {
        alpha = fabs ( tdata[1] - tdata[0] )
                / ( fabs ( tdata[1] - tdata[0] ) + fabs ( tdata[2] - tdata[1] ) );

        mbasis = basis_matrix_overhauser_nul ( alpha );

        yval = basis_matrix_tmp ( left, 3, mbasis, ndata, tdata, ydata, tval );
    } else if ( left < ndata-1 ) {
        alpha = fabs ( tdata[left]   - tdata[left-1]   )
                / ( fabs ( tdata[left]   - tdata[left-1]   )
                    + fabs ( tdata[left+1] - tdata[left]     ) );

        beta = fabs ( tdata[left+1] - tdata[left]     )
               / ( fabs ( tdata[left]   - tdata[left-1]   )
                   + fabs ( tdata[left+1] - tdata[left]     ) );

        mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

        yval = basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, tval );
    } else if ( left == ndata-1 ) {
        beta = fabs ( tdata[ndata-1] - tdata[ndata-2] )
               / ( fabs ( tdata[ndata-2] - tdata[ndata-3] )
                   + fabs ( tdata[ndata-1] - tdata[ndata-2] ) );

        mbasis = basis_matrix_overhauser_nur ( beta );

        yval = basis_matrix_tmp ( left, 3, mbasis, ndata, tdata, ydata, tval );
    } else {
        cout << "\n";
        cout << "SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!\n";
        cout << "  Nonsensical value of LEFT = " << left << "\n";
        cout << "  but 0 < LEFT < NDATA = " << ndata << "\n";
        cout << "  is required.\n";
        exit ( 1 );
    }

    delete [] mbasis;

    return yval;
}
//******************************************************************************

double spline_overhauser_uni_val ( int ndata, double tdata[], double ydata[],
                                   double tval )

//******************************************************************************
//
//  Purpose:
//
//    SPLINE_OVERHAUSER_UNI_VAL evaluates the uniform Overhauser spline.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDATA, the number of data points.
//    NDATA must be at least 3.
//
//    Input, double TDATA[NDATA], the abscissas of the data points.
//    The values of TDATA are assumed to be distinct and increasing.
//    This routine also assumes that the values of TDATA are uniformly
//    spaced; for instance, TDATA(1) = 10, TDATA(2) = 11, TDATA(3) = 12...
//
//    Input, double YDATA[NDATA], the data values.
//
//    Input, double TVAL, the value where the spline is to
//    be evaluated.
//
//    Output, double SPLINE_OVERHAUSER_UNI_VAL, the value of the spline at TVAL.
//
{
    int left;
    double *mbasis = NULL;
    int right;
    double yval = 0.0;
    //
    //  Check NDATA.
    //
    if ( ndata < 3 ) {
        cout << "\n";
        cout << "SPLINE_OVERHAUSER_UNI_VAL - Fatal error!\n";
        cout << "  NDATA < 3.\n";
        exit ( 1 );
    }
    //
    //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    //
    dvec_bracket ( ndata, tdata, tval, &left, &right );
    //
    //  Evaluate the spline in the given interval.
    //
    if ( left == 1 ) {

        mbasis = basis_matrix_overhauser_uni_l ( );

        yval = basis_matrix_tmp ( left, 3, mbasis, ndata, tdata, ydata, tval );
    } else if ( left < ndata-1 ) {
        mbasis = basis_matrix_overhauser_uni ( );

        yval = basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, tval );
    } else if ( left == ndata-1 ) {
        mbasis = basis_matrix_overhauser_uni_r ( );

        yval = basis_matrix_tmp ( left, 3, mbasis, ndata, tdata, ydata, tval );

    }
    if(mbasis != NULL)
        delete [] mbasis;

    return yval;
}
//**********************************************************************

void spline_overhauser_val ( int ndim, int ndata, double tdata[], double ydata[],
                             double tval, double yval[] )

//**********************************************************************
//
//  Purpose:
//
//    SPLINE_OVERHAUSER_VAL evaluates an Overhauser spline.
//
//  Discussion:
//
//    Over the first and last intervals, the Overhauser spline is a
//    quadratic.  In the intermediate intervals, it is a piecewise cubic.
//    The Overhauser spline is also known as the Catmull-Rom spline.
//
//  Modified:
//
//    11 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    H Brewer and D Anderson,
//    Visual Interaction with Overhauser Curves and Surfaces,
//    SIGGRAPH 77, pages 132-137.
//
//    E Catmull and R Rom,
//    A Class of Local Interpolating Splines,
//    in Computer Aided Geometric Design,
//    edited by R Barnhill and R Reisenfeld,
//    Academic Press, 1974, pages 317-326.
//
//    David Rogers and Alan Adams,
//    Mathematical Elements of Computer Graphics,
//    McGraw Hill, 1990, Second Edition, pages 278-289.
//
//  Parameters:
//
//    Input, int NDIM, the dimension of a single data point.
//    NDIM must be at least 1.
//
//    Input, int NDATA, the number of data points.
//    NDATA must be at least 3.
//
//    Input, double TDATA[NDATA], the abscissas of the data points.  The
//    values in TDATA must be in strictly ascending order.
//
//    Input, double YDATA[NDIM*NDATA], the data points corresponding to
//    the abscissas.
//
//    Input, double TVAL, the abscissa value at which the spline
//    is to be evaluated.  Normally, TDATA[0] <= TVAL <= T[NDATA-1], and
//    the data will be interpolated.  For TVAL outside this range,
//    extrapolation will be used.
//
//    Output, double YVAL[NDIM], the value of the spline at TVAL.
//
{
    int i;
    int left;
    int order;
    int right;
    double *yl;
    double *yr;
    //
    //  Check.
    //
    dvec_order_type ( ndata, tdata, &order );

    if ( order != 2 ) {
        cout << "\n";
        cout << "SPLINE_OVERHAUSER_VAL - Fatal error!\n";
        cout << "  The data abscissas are not strictly ascending.\n";
        exit ( 1 );
    }

    if ( ndata < 3 ) {
        cout << "\n";
        cout << "SPLINE_OVERHAUSER_VAL - Fatal error!\n";
        cout << "  NDATA < 3.\n";
        exit ( 1 );
    }
    //
    //  Locate the abscissa interval T[LEFT], T[LEFT+1] nearest to or
    //  containing TVAL.
    //
    dvec_bracket ( ndata, tdata, tval, &left, &right );
    //
    //  Evaluate the "left hand" quadratic defined at
    //  T[LEFT-1], T[LEFT], T[RIGHT].
    //
    yl = new double[ndim];
    yr = new double[ndim];

    if ( 0 < left-1 ) {
        parabola_val2 ( ndim, ndata, tdata, ydata, left-1, tval, yl );
    }
    //
    //  Evaluate the "right hand" quadratic defined at
    //  T[LEFT], T[RIGHT], T[RIGHT+1].
    //
    if ( right+1 <= ndata ) {
        parabola_val2 ( ndim, ndata, tdata, ydata, left, tval, yr );
    }
    //
    //  Blend the quadratics.
    //
    if ( left == 1 ) {
        for ( i = 0; i < ndim; i++ ) {
            yval[i] = yr[i];
        }
    } else if ( right < ndata ) {
        for ( i = 0; i < ndim; i++ ) {
            yval[i] = (
                          ( tdata[right-1] - tval                 ) * yl[i]
                          + (                  tval - tdata[left-1] ) * yr[i] )
                      / ( tdata[right-1]        - tdata[left-1] );
        }
    } else {
        for ( i = 0; i < ndim; i++ ) {
            yval[i] = yl[i];
        }
    }

    delete [] yl;
    delete [] yr;

    return;
}
//******************************************************************************

void spline_quadratic_val ( int ndata, double tdata[], double ydata[],
                            double tval, double *yval, double *ypval )

//******************************************************************************
//
//  Purpose:
//
//    SPLINE_QUADRATIC_VAL evaluates a piecewise quadratic spline at a point.
//
//  Discussion:
//
//    Because of the simple form of a piecewise quadratic spline,
//    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
//    evaluate the spline at any point.  No processing of the data
//    is required.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDATA, the number of data points defining the spline.
//    NDATA should be odd and at least 3.
//
//    Input, double TDATA[NDATA], YDATA[NDATA], the values of the independent
//    and dependent variables at the data points.  The values of TDATA should
//    be distinct and increasing.
//
//    Input, double TVAL, the point at which the spline is to be evaluated.
//
//    Output, double *YVAL, *YPVAL, the value of the spline and its first
//    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
//    equal to TDATA(I) for some I.
//
{
    double dif1;
    double dif2;
    int left;
    int right;
    double t1;
    double t2;
    double t3;
    double y1;
    double y2;
    double y3;

    if ( ndata < 3 ) {
        cout << "\n";
        cout << "SPLINE_QUADRATIC_VAL - Fatal error!\n";
        cout << "  NDATA < 3.\n";
        exit ( 1 );
    }

    if ( ndata % 2 == 0 ) {
        cout << "\n";
        cout << "SPLINE_QUADRATIC_VAL - Fatal error!\n";
        cout << "  NDATA must be odd.\n";
        exit ( 1 );
    }
    //
    //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
    //  nearest to, TVAL.
    //
    dvec_bracket ( ndata, tdata, tval, &left, &right );
    //
    //  Force LEFT to be odd.
    //
    if ( left % 2 == 0 ) {
        left = left - 1;
    }
    //
    //  Copy out the three abscissas.
    //
    t1 = tdata[left-1];
    t2 = tdata[left  ];
    t3 = tdata[left+1];

    if ( t2 <= t1 || t3 <= t2 ) {
        cout << "\n";
        cout << "SPLINE_QUADRATIC_VAL - Fatal error!\n";
        cout << "  T2 <= T1 or T3 <= T2.\n";
        exit ( 1 );
    }
    //
    //  Construct and evaluate a parabolic interpolant for the data
    //  in each dimension.
    //
    y1 = ydata[left-1];
    y2 = ydata[left  ];
    y3 = ydata[left+1];

    dif1 = ( y2 - y1 ) / ( t2 - t1 );

    dif2 = ( ( y3 - y1 ) / ( t3 - t1 )
             - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 );

    *yval = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 );
    *ypval = dif1 + dif2 * ( 2.0E+00 * tval - t1 - t2 );

    return;
}
//**********************************************************************

void timestamp ( void )

//**********************************************************************
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
#define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    size_t len;
    time_t now;

    now = time ( NULL );
    tm = localtime ( &now );

    len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

    cout << time_buffer << "\n";

    return;
#undef TIME_SIZE
}

