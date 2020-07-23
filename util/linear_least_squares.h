#ifndef _linear_least_squares_h
#define _linear_least_squares_h

#include "TMatrix.h"
#include "TVectorD.h"

/** This function solves a linear-least squares problem. 
 *
 * Input variables are:
 *
 * residual: A TVectorD of size $n$ that contains the residual values 
 *           (i.e. measurement-prediction) 
 *
 * Jacobian: A TMatrixD, size of n x m with the values of $\left.\partial f/\partial
 *           p\right|_{p_0 = 0}$, where $f(p)$ is the function that predicts the
 *           measurement for a set of given parameters $p$, and the Jacobian is 
 *           evaluated for $p = p_0 = 0$.
 *
 * covariance: A TMatrixD of size n x n which is the covariance matrix of the measurements.
 *
 * Return value:
 *
 * TVectorD of size m, which contains the optimal solution for the parameters $p$.
 *
 * BUGS:
 *
 * The program uses matrix inversion which is very slow with a big matrix.
 *
 * In case a problem appears during matrix inversion, an error message is printed but otherwise 
 * the function continues without e.g. throwing an exception.
 *
 */
TVectorD linear_least_squares(const TVectorD & residual, const TMatrixD & Jacobian, const TMatrixD & covariance);

#endif
