#include<stdio.h>
#include<stdlib.h>
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fq.h>
#include<gmp.h>

/* Polynômes de Division:
 * (1) Le l-ième polynôme de division Psi_l a pour racines 
 * les abscisses des points d'ordre l.
 * (2) Psi_1 = 1; Psi_2 = 2*y + a1*x + a3;
 * (3) Pour l>=5:
 * 	Si l=2*k: Psi_2k = \frac{...}{...}
 * 	Si l=2*k+1: ...
 *	https://yx7.cc/docs/tum/thesis_schoof.pdf
 *
 * Pour calculer [l]P
*/


