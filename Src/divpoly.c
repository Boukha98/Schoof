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
 * (3) Pour i>=5:
 * 	Si i=2*m: Psi_i = \frac{...}{...}
 * 	Si i=2*m+1: ...
 *	https://yx7.cc/docs/tum/thesis_schoof.pdf
 *
 * Pour calculer [i]P
*/

void div_poly(fq_poly_t *Psi, elliptic_curve E, 
		unsigned long k, fq_ctx_t ctx){
	fq_t tmp, tmp1;
	slong i,m;
	// On initialise Psi_0...4
	for(i=0; i<3; i++){
	     	fq_set_ui(tmp, i, ctx);
       		fq_poly_set_fq(Psi[i], tmp, ctx);
 	}
	
	/*fq_poly_zero(Psi[0], ctx); //Psi_0 = 0
	fq_poly_one(Psi[1], ctx); //Psi_1 = 1
	fq_poly_set_fq(Psi[2], 2, ctx); //Psi_2 = 2*/
	
	//Psi_3 = 3x^4+6ax^2+12bx-a^2
	fq_set_ui(tmp, 3, ctx);
	fq_poly_set_coeff(Psi[3], 4, tmp, ctx);//3*x^4

	fq_mul_ui(tmp, E.a, 6, ctx);//6*a
	fq_poly_set_coeff(Psi[3], 2, tmp, ctx);//6*a*x^2

	fq_mul_ui(tmp, E.b, 12, ctx);//12*b
	fq_poly_set_coeff(Psi[3], 1, tmp, ctx);//12*b*x

	fq_sqr(tmp, E.a, ctx);
	fq_neg(tmp, tmp, ctx);/*pour a^3)*/fq_mul(tmp1, tmp, E.a, ctx);
	fq_poly_set_coeff(Psi[3], 0, tmp, ctx);

	//Psi_4 = 4y(x^6+5ax^4+20bx^3-5a^2x^2-4abx-8b^2-a^3)
	fq_mul_ui(tmp, tmp, 20, ctx);
	fq_poly_set_coeff(Psi[4], 2, tmp, ctx);// on a deja calculé a^2

	fq_set_ui(tmp, 4, ctx);
	fq_poly_set_coeff(Psi[4], 6, tmp, ctx);

	fq_mul_ui(tmp, E.a, 20, ctx);
	fq_poly_set_coeff(Psi[4], 4, tmp, ctx);
	
	fq_mul_ui(tmp, E.b, 80, ctx);
	fq_poly_set_coeff(Psi[4], 3, tmp, ctx);

	fq_mul(tmp, E.a, E.b, ctx);
	fq_mul_ui(tmp, tmp, 16, ctx);
	fq_neg(tmp, tmp, ctx);
	fq_poly_set_coeff(Psi[4], 1, tmp, ctx);

	fq_sqr(tmp, E.b, ctx);
	fq_mul_ui(tmp, E.b, 8, ctx);
	fq_sub(tmp, tmp1, tmp, ctx);
	fq_mul_ui(tmp, tmp, 4, ctx);
	fq_poly_set_coeff(Psi[4], 0, tmp, ctx);
	
	//Récurrence
	for(i=5; i<k; i++){
		//Cas i pair: i=2*m
		if(!(i&1)){
			m = i>>1;
			fq_poly_sqr(Psi[i], Psi[m+1], ctx);
			fq_poly_mul(Psi[i], Psi[i], Psi[m-2], ctx);
	
			fq_poly_sqr(poly, Psi[m-1], ctx);
			fq_poly_mul(poly, poly, Psi[m+2], ctx);

			fq_poly_sub(poly, poly, Psi[i], ctx);
			fq_set_ui(tmp, 2, ctx);
			fq_inv(tmp, tmp, ctx);
			fq_poly_scalar_mul_fq(Psi[i], Psi[m], tmp, ctx);
			fq_poly_mul(Psi[i], Psi[i], poly, ctx);
		}
		//Cas i impair: i=2*m+1
		else{
			m = i-1;
			m>>=2;
			fq_poly_pow(Psi[i], Psi[m], 3, ctx);
			fq_poly_mul(Psi[i], Psi[i], Psi[m+2], ctx);
			fq_poly_pow(poly, Psi[m+1], 3, ctx);
			fq_poly_mul(poly, poly, Psi[m-1], ctx);
			fq_poly_sub(Psi[i], Psi[i], poly, ctx);
		}
	}
}	
