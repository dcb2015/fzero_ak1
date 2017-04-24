// fzero_ak1.cpp
//
// This program seeks a zero of a real function in one dimension.
// The function to be zeroed is defined in the subroutine f
//
// The main root-finding routine is a translation of the FORTRAN routine FZERO.F,
// written by L. F. Shampine (SNLA) and H. A. Watts (SNLA), based upon a method by T. J. Dekker.
// FZERO.F is part of the SLATEC library of programs, and its original code (written in FORTRAN) can be viewed there.
// 
// References:
//
// Shampine, L.F. (SNLA) and H.A.Watts (SNLA)
// "FZERO, A Root-solving Code"   Report SC - TM - 70 - 631
// Sandia Laboratories
// September, 1970
//
// Dekker, T.J.
// "Finding a Zero by Means of Successive Linear Interpolation"
// "Constructive Aspects of the Fundamental Theorem of Algebra"
// edited by B.Dejon and P.Henrici
// Wiley - Interscience
// 1969
//
// To distinguish this version from other possible translations, an '_ak1' suffix has been appended to its name.
//
// The particular equation zeroed in this project (to include an example), is 
//
// f(x) = { sech{ acosh[(1/D)^(1/x)] * 0.7861513777574233 } }^x - 0.805917329564
//
// where D is the Dottie Number (= 0.73908513321516064)
// The value 0.7861513777574233 is the eccenticity of an ellipse currently being examined, and
// the value 0.805917329564 is the value of sin(E) at the quarter-period for this ellipse.
//
// The program requires four inputs:
//
// dumvar.b	the lower bound of the interval over which the zero will be sought. Upon return from fzero_ak1,
//			this variable contains the zero
// dumvar.c	the upper bound of the interval over which the zero will be sought. The sign of the function should change
//			over the interval [b, c]
// dumvar.ae	absolute error
// dumvar.re	relative error
// In the present program, ae and re are assigned the same value: DBL_EPSILON
//
// 23 April 2017
//
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
//

#include <iostream>
#include <cctype>
#include <cmath>
#include <cfloat>

using namespace std;

#define ECC 0.7861513777574233  // The eccentricity for the present example
#define DOTTIE 0.73908513321516064 // The Dottie Number

struct ZTYPE	// mAnom not used for the present equation
{
	double b, c, re, ae;
	short iflag, kount;
	double mAnom;
};

double f(double x);
double sign(double y);
void outCodeandKount(short erFlag, short fxneval);
void fzero_ak1(ZTYPE *zd, double ecc);

double f(double x){	//The function to be zeroed

	static const double INV_DOTTIE = 1.0 / DOTTIE; // Inverse of the Dottie Number
	static const double SET4 = 0.805917329564; // sin(E) at T/4 for ecc = 0.7861...

	double dum1 = 1.0 / x, dum2, k, temp1, temp2;

	k = acosh(pow(INV_DOTTIE, dum1));
	temp1 = ECC * k;
	dum2 = 1.0 / cosh(temp1);
	temp2 = pow(dum2, x);
	return temp2 - SET4;
}

double sign(double y){			//If y < 0, return -1, else +1
	return ((y < 0) ? -1 : 1);
}

void outCodeandKount(short erFlag, short fxneval){ //Output error code and ikount
	cout << "\nThe number of function calls was " << fxneval << ".\n\n";
	switch (erFlag)
	{
	case 1:	cout << "Error Code = 1: The zero is within the requested tolerance, the\n";
		cout << "interval has collapsed to the requested tolerance, the function\n";
		cout << "changes sign over the interval, and the function decreased in\n";
		cout << "magnitude as the interval collapsed.\n";
		break;
	case 2:	cout << "Error Code = 2: A zero has been found, but the interval has not\n";
		cout << "collapsed to the requested tolerance.\n";
		break;
	case 3:	cout << "Error Code = 3: More than MAXIT (=100) function evaluations used.\n";
		break;
	default:	cout << "Invalid Error Code.\n";
	} //End switch
}

void fzero_ak1(ZTYPE *zd, double ecc) {
	/*	This function determines the zero of the function specified in f, above. The method
	employed is an efficient combination of bisection and secant rules.

	Parameters
	b				lower bound of interval in which zero to be found;
					returns as zero to function
	c				upper bound of interval in which zero to be found
	re				relative error
	ae				absolute error
	iflag			status code
	=1	The zero is within the requested tolerance, the interval has collapsed
	to the requested tolerance, the function changes sign over the interval,
	and the function decreased in magnitude as the interval collapsed.
	=2	A zero has been found, but the interval has not collapsed to the
	requested tolerance.
	=3	Possibly near a singular point. The interval has collapsed to the
	requested tolerance and the function changes sign over the interval,
	but the magnitude of the function increased as the interval collapsed.
	=4	The function does not change sign over the specified interval, which
	has collapsed to the requested tolerance. Possibly near a minimum of
	the function or a zero of even multiplicity.
	=5	More than MAXIT function evaluations used.
	kount			number of function calls

	Authors:	Shampine, L.F., SNLA
	Watts, H.A., SNLA */

	static const short MAXIT = 100;
	static const double er = DBL_EPSILON; //Machine epsilon for type double
	double a, acbs, acmb, cmb, fa, fb, fc, fz, p, q, rw, t, tol, z;
	short ic = 0;

	//Initialize

	rw = DBL_EPSILON;  //((zd->re > er) ? zd->re : er);
	z = zd->b + ecc / 2 / 2; //In this case zd->c - zd->b = ecc/2
	fc = fz = f(z); t = zd->b; fb = f(t);
	zd->kount = 2;
	if (sign(fz) == sign(fb))
	{
		t = zd->c; fc = f(t); zd->kount = 3;
		if (sign(fz) != sign(fc))	{ zd->b = z; fb = fz; }
	}
	else zd->c = z;
	a = zd->c; fa = fc; acbs = ecc / 2 / 2; //In this case, fabs(zd->b - zd->c) = ecc/4, since z was bisected in line 1;

	do
	{
		if (fabs(fc) < fabs(fb))		//Interchange if necessary
		{
			a = zd->b; fa = fb; zd->b = zd->c; fb = fc; zd->c = a; fc = fa;
		}
		cmb = (zd->c - zd->b) / 2; acmb = fabs(cmb);
		tol = fabs(zd->b);
		tol++;
		tol *= rw;  //tol = rw*fabs(zd->b) + zd->ae

		//Test-stopping criteria and function count

		if (acmb <= tol) { zd->iflag = 1;	return; }	//Normal completion

		if (fb == 0){ zd->iflag = 2; return; }

		/*Calculate new iterate implicitly as b + p/q, where p is arranged to be
		>= 0. This implicit form is used to prevent overflow.*/

		p = (zd->b - a)*fb; q = fa - fb;
		if (p < 0) { p = -p; q = -q; }

		/*Update a and check for satisfactory reduction in the size of the bracketing
		interval. If not, perform bisection.*/

		a = zd->b; fa = fb; ic++;
		if ((ic >= 4) && (8 * acmb >= acbs))
			zd->b = (zd->c + zd->b) / 2;				//Use bisection
		else
		{
			if (ic >= 4) { ic = 0; acbs = acmb; }
			if (p <= tol*fabs(q))			//Test for too small a change
				zd->b += tol*sign(cmb);
			else							//Root between b and (b + c)/2
			{
				if (p < cmb*q) zd->b += p / q;	//Use secant rule
				else zd->b = (zd->c + zd->b) / 2;
			}
		}	//Use bisection.
		//Have now computed new iterate,b.
		fb = f(zd->b); (zd->kount)++;

		//Decide if next step interpolation or extrapolation.

		if (sign(fb) == sign(fc)) { zd->c = a; fc = fa; }
	} while (zd->kount < MAXIT);				//End while loop
	if (zd->kount >= MAXIT) zd->iflag = 3;	//Too many iterations
	return;
}											//End fzero

int main(){

	char rflag = 0;				//Readiness flag

	cout << "                  fzero_ak1 (23 April 2017)\n";
	cout << "=================================================================\n";
	cout << "This program calculates the power of a sech function which approximates the\n";
	cout << "value of sin(E) at the quarter period for\n";
	cout << "an ellipse of eccentricity e = 0.7861513777574233. \n";
	cout << "At the quarter period, sin(E) = 0.805917329564. \n";
	cout << "\nEverything ready? If yes, Enter y.\n";
	cout << "Otherwise Enter any other key.\n";
	cin >> rflag;

	if (toupper(rflag) == 'Y')
	{
		ZTYPE dumvar;			//Variable for fzero

		cout << "\n=================================================================\n";
		cout.precision(DBL_DIG);

		dumvar.b = 0.125;
		dumvar.c = 0.25;
		dumvar.ae = dumvar.re = DBL_EPSILON;

		fzero_ak1(&dumvar, ECC);
		cout << "\nThe solution p-value is " << dumvar.b << ".\n";
		outCodeandKount(dumvar.iflag, dumvar.kount);
	}
	else cout << "\nNot ready. Try again when ready with information.\n";
	cout << "\nEnter any key to continue.\n";
	cin >> rflag;
	return 0;
}
