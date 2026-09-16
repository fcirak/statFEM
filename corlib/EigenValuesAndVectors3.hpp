// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file EigenValuesAndVectors3.hpp

#ifndef corlib_eigenValuesAndVectors3_h
#define corlib_eigenValuesAndVectors3_h

#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <corlib/verify.hpp>
#include <corlib/linalg.hpp>

//------------------------------------------------------------------------------
namespace corlib{

    struct EigenValuesAndVectors3;

    namespace eigenX = corlib::eigenX;

    namespace detail_{

        //----------------------------------------------------------------------
        /** Symmetric Householder reduction to tridiagonal form.
         * This code has been taken from http://barnesc.blogspot.com/ and 
         * adapted. Its original version claims:
         *
         *    "Eigen decomposition code for symmetric 3x3 matrices, copied 
         *     from the public domain Java Matrix library JAMA."
         *
         * Comments have been kept as in the source.
         */
        static void tred2(double V[3][3], double d[3], double e[3]) 
        {

            //  This is derived from the Algol procedures tred2 by
            //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
            //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutine in EISPACK.

            for (int j = 0; j < 3; j++) {
                d[j] = V[2][j];
            }

            // Householder reduction to tridiagonal form.

            for (int i = 2; i > 0; i--) {

                // Scale to avoid under/overflow.
                
                double scale = 0.0;
                double h = 0.0;
                for (int k = 0; k < i; k++) {
                    scale = scale + fabs(d[k]);
                }
                if (scale == 0.0) {
                    e[i] = d[i-1];
                    for (int j = 0; j < i; j++) {
                        d[j] = V[i-1][j];
                        V[i][j] = 0.0;
                        V[j][i] = 0.0;
                    }
                } 
                else {
                    // Generate Householder vector.
                    for (int k = 0; k < i; k++) {
                        d[k] /= scale;
                        h += d[k] * d[k];
                    }
                    double f = d[i-1];
                    double g = sqrt(h);
                    if (f > 0) {
                        g = -g;
                    }
                    e[i] = scale * g;
                    h = h - f * g;
                    d[i-1] = f - g;
                    for (int j = 0; j < i; j++) {
                        e[j] = 0.0;
                    }

                    // Apply similarity transformation to remaining columns.
                    for (int j = 0; j < i; j++) {
                        f = d[j];
                        V[j][i] = f;
                        g = e[j] + V[j][j] * f;
                        for (int k = j+1; k <= i-1; k++) {
                            g += V[k][j] * d[k];
                            e[k] += V[k][j] * f;
                        }
                        e[j] = g;
                    }
                    f = 0.0;
                    for (int j = 0; j < i; j++) {
                        e[j] /= h;
                        f += e[j] * d[j];
                    }
                    double hh = f / (h + h);
                    for (int j = 0; j < i; j++) {
                        e[j] -= hh * d[j];
                    }
                    for (int j = 0; j < i; j++) {
                        f = d[j];
                        g = e[j];
                        for (int k = j; k <= i-1; k++) {
                            V[k][j] -= (f * e[k] + g * d[k]);
                        }
                        d[j] = V[i-1][j];
                        V[i][j] = 0.0;
                    }
                }
                d[i] = h;
            }

            // Accumulate transformations.
            for (int i = 0; i < 2; i++) {
                V[2][i] = V[i][i];
                V[i][i] = 1.0;
                double h = d[i+1];
                if (h != 0.0) {
                    for (int k = 0; k <= i; k++) {
                        d[k] = V[k][i+1] / h;
                    }
                    for (int j = 0; j <= i; j++) {
                        double g = 0.0;
                        for (int k = 0; k <= i; k++) {
                            g += V[k][i+1] * V[k][j];
                        }
                        for (int k = 0; k <= i; k++) {
                            V[k][j] -= g * d[k];
                        }
                    }
                }
                for (int k = 0; k <= i; k++) {
                    V[k][i+1] = 0.0;
                }
            }
            for (int j = 0; j < 3; j++) {
                d[j] = V[2][j];
                V[2][j] = 0.0;
            }
            V[2][2] = 1.0;
            e[0] = 0.0;
        } 

        //----------------------------------------------------------------------
        // Symmetric tridiagonal QL algorithm.
        static void tql2(double V[3][3], double d[3], double e[3]) 
        {
            //  This is derived from the Algol procedures tql2, by
            //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
            //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutine in EISPACK.
            for (int i = 1; i < 3; i++) {
                e[i-1] = e[i];
            }
            e[2] = 0.0;

            double f = 0.0;
            double tst1 = 0.0;
            double eps = pow(2.0,-52.0);
            for (int l = 0; l < 3; l++) {
                // Find small subdiagonal element
                tst1 = std::max(tst1,fabs(d[l]) + fabs(e[l]));
                int m = l;
                while (m < 3) {
                    if (fabs(e[m]) <= eps*tst1) {
                        break;
                    }
                    m++;
                }

                // If m == l, d[l] is an eigenvalue,
                // otherwise, iterate.
                if (m > l) {
                    int iter = 0;
                    do {
                        iter = iter + 1;  // (Could check iteration count here.)

                        // Compute implicit shift
                        double g = d[l];
                        double p = (d[l+1] - g) / (2.0 * e[l]);
                        double r = std::sqrt( p*p + 1. );
                        if (p < 0) {
                            r = -r;
                        }
                        d[l] = e[l] / (p + r);
                        d[l+1] = e[l] * (p + r);
                        double dl1 = d[l+1];
                        double h = g - d[l];
                        for (int i = l+2; i < 3; i++) {
                            d[i] -= h;
                        }
                        f = f + h;

                        // Implicit QL transformation.
                        p = d[m];
                        double c = 1.0;
                        double c2 = c;
                        double c3 = c;
                        double el1 = e[l+1];
                        double s = 0.0;
                        double s2 = 0.0;
                        for (int i = m-1; i >= l; i--) {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * e[i];
                            h = c * p;
                            r = std::sqrt( p*p + e[i] * e[i] );
                            e[i+1] = s * r;
                            s = e[i] / r;
                            c = p / r;
                            p = c * d[i] - s * g;
                            d[i+1] = h + s * (c * g + s * d[i]);

                            // Accumulate transformation.
                            for (int k = 0; k < 3; k++) {
                                h = V[k][i+1];
                                V[k][i+1] = s * V[k][i] + c * h;
                                V[k][i] = c * V[k][i] - s * h;
                            }
                        }
                        p = -s * s2 * c3 * el1 * e[l] / dl1;
                        e[l] = s * p;
                        d[l] = c * p;

                        // Check for convergence.
                    } while (fabs(e[l]) > eps*tst1);
                }
                d[l] = d[l] + f;
                e[l] = 0.0;
            }
  
            // Sort eigenvalues and corresponding vectors.
            for (int i = 0; i < 2; i++) {
                int k = i;
                double p = d[i];
                for (int j = i+1; j < 3; j++) {
                    if (d[j] < p) {
                        k = j;
                        p = d[j];
                    }
                }
                if (k != i) {
                    d[k] = d[i];
                    d[i] = p;
                    for (int j = 0; j < 3; j++) {
                        p = V[j][i];
                        V[j][i] = V[j][k];
                        V[j][k] = p;
                    }
                }
            }
        }

        //----------------------------------------------------------------------
        /** Simple method to compute the eigenvalues only.
         *  Uses the method given in 
         *  O.K. Smith, 'Eigenvalues of a Symmetric 3x3 Matrix', Communications
         *  of the ACM (4), 1961, p. 168.
         */
        void eig3x3( const eigenX::MatrixSd<3,3> & sig,
                     eigenX::VectorSd<3> & eval )
        {
            // compute first invariant of sigma
            const double Isig =(sig(0,0) + sig(1,1) + sig(2,2) ) / 3.;

            // compute deviator
            const eigenX::MatrixSd<3,3> dev =
                sig - Isig * Eigen::MatrixXd::Identity(3,3);
    
            // compute second invariant of deviator
            const double IIdev = ( dev(0,0)*dev(0,0) + dev(1,1)*dev(1,1) + 
                                   dev(2,2)*dev(2,2) + 2. * dev(0,1)*dev(0,1) + 
                                   2. * dev(0,2)*dev(0,2) + 2. * dev(1,2)*dev(1,2) ) / 6.;

            // compute half the determinant (third invariant) of deviator
            const double IIIdev = 0.5 * corlib::determinant( dev );

            // compute strange angle
            double aux =  IIdev*IIdev*IIdev - IIIdev*IIIdev;
            if ( aux < 0. ) aux = 0.;
            double phi = atan2( sqrt(aux), IIIdev )/ 3.;

            // readjust
            if ( phi < 0. ) phi += M_PI/3.;

            // eigenvalues
            eval(0) = Isig + 2. * sqrt(IIdev) * cos(phi);
            eval(1) = Isig - sqrt(IIdev) * (cos(phi) + sqrt(3.) * sin(phi));
            eval(2) = Isig - sqrt(IIdev) * (cos(phi) - sqrt(3.) * sin(phi));

            return;
        }

    }
}

//------------------------------------------------------------------------------
/** \brief Eigen-analysis of 3x3 symmetric matrices
 *  \details Given a symmetric 3x3 matrix, this object allows for the 
 *  computation of its eigenvalues and -vectors. 
 *  The computation of eigenvalues only uses a very simple technique as cited
 *  above. The computation of eigenvectors (and eigenvalues) is based on the 
 *  successive use of a Householder reduction and a QL-transform. The functions
 *  for the exclusive computation of the eigen-thingees can be easily wrapped
 *  around a functor which, e.g., computes the stress tensor of a an element.
 */
struct corlib::EigenValuesAndVectors3
{
    typedef eigenX::VectorSd<3>   Vec3;
    typedef eigenX::MatrixSd<3,3> Mat3x3;

    //! Compute eigenvalues only
    static Vec3 eigenValues( const Mat3x3 & matrix ) {
        // check for symmetry
        FTL_VERIFY( corlib::EigenValuesAndVectors3::symmetryCheck_( matrix, 1.e-8 ) );

        Vec3 res;
        detail_::eig3x3( matrix, res );
        return res;

        // alternative:
        //Mat3x3 dummy; Vec3   evals;
        //corlib::EigenValuesAndVectors3::callImplementation( matrix, evals, dummy );
        //return evals;
    }
    
    //! Compute eigenvectors
    static Mat3x3 eigenVectors( const Mat3x3 & matrix ) 
    {
        // check for symmetry
        FTL_VERIFY( corlib::EigenValuesAndVectors3::symmetryCheck_( matrix, 1.e-8  ) );

        Mat3x3 evecs;
        Vec3   dummy;
        corlib::EigenValuesAndVectors3::callImplementation_( matrix, dummy, evecs );
        return evecs;
    }

    //! Compute eigenvalues and eigenvectors together
    static void eigenValuesAndVectors( const Mat3x3 & matrix, 
                                       Vec3   & evals, 
                                       Mat3x3 & evecs )
    {
        // check for symmetry
        FTL_VERIFY( corlib::EigenValuesAndVectors3::symmetryCheck_( matrix, 1.e-8 ) );

        corlib::EigenValuesAndVectors3::callImplementation_( matrix, evals, evecs );
        return;
    }

private:
    //! Check of symmetry
    static bool symmetryCheck_( const Mat3x3 & matrix, const double tol )
    {
        const double normMat = matrix.norm( ); // Frobenius norm of matrix
        if ( normMat < tol ) return true; // too small entries
        const double relTol =  normMat * tol;
        return ( matrix - matrix.transpose( ) ).norm( ) < relTol;
    }

    //! Convenience function to call the functions in detail_::
    static void callImplementation_( const Mat3x3 & matrix, 
                                    Vec3   & evals, 
                                    Mat3x3 & evecs )
    {
        // copy matrix to c-array
        double V[3][3], d[3], e[3];
        for ( unsigned i = 0; i < 3; i ++ ) 
            for ( unsigned j = 0; j < 3; j ++ )
                V[i][j] = matrix(i,j);
        
        // compute Householder reduction
        detail_::tred2( V, d, e );
        // do QL decomposition
        detail_::tql2(  V, d, e );

        // copy results to output in ublas-format
        for ( unsigned i = 0; i < 3; i ++ ) {
            evals(i) = d[i];
            for ( unsigned j = 0; j < 3; j ++ ) {
                evecs(i,j) = V[i][j];
            }
        }
        return;
    }

};

#endif
