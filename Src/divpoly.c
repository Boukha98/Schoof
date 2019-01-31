#include<stdio.h>
#include<stdlib.h>
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fq.h>
#include<gmp.h>
#include<flint/fq_poly.h>
//#include"schoof.h"


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
typedef struct{
	fq_t x, y; // Coordonnées du point
}point;

//initialisation
void fq_init_point(point P, fq_ctx_t ctx){
	fq_init(P.x, ctx);
	fq_init(P.y, ctx);
}

//Courbe elliptique y^2 = x^3+a*x+b mod p
typedef struct{
	fmpz_t p;
	fq_t a, b; 
}elliptic_curve;

//initialisation
void fq_init_curve(elliptic_curve E, fq_ctx_t ctx){
	fq_init(E.a, ctx);
	fq_init(E.b, ctx);
}

void div_poly(fq_poly_t *Psi, elliptic_curve E, 
		signed long k, fq_ctx_t ctx, point P){
	fq_t tmp, tmp1;
	fq_poly_t poly;
	fq_poly_init(poly, ctx);
	fq_init(tmp, ctx);fq_init(tmp1, ctx);
	slong i,m;
	//initialisation de P
	//fq_init_point(P, ctx);

	//Polynome y^2 = x^3 + a*x + b
	fq_poly_t y2;
	fq_poly_init(y2, ctx);
	fq_set_ui(tmp, 1, ctx);
	fq_poly_set_coeff(y2, 3, tmp, ctx);//x^3
	
	fq_set(tmp, E.a, ctx);
	fq_poly_set_coeff(y2, 1, tmp, ctx);//a*x
	
	fq_set(tmp, E.b, ctx);
	fq_poly_set_coeff(y2, 0, tmp, ctx);//b

	/*on verifie si P \in E
	fq_sqr(tmp, P.y, ctx);
	fq_poly_evaluate_fq(tmp1, y2, P.x, ctx);
	if(!fq_equal(tmp, tmp1, ctx)){
		fprintf(stderr, "Ce n'est pas un point de la courbe %,%s\n", P.x, P.y); 		exit(EXIT_FAILURE);
	}*/

	// On initialise Psi_0...4
	for(i = 0; i < k + 1; i++){
 		fq_poly_init(Psi[i], ctx);
	}
	/*for(i=0; i<3; i++){
	     	fq_set_ui(tmp, i, ctx);
       		fq_poly_set_fq(Psi[i], tmp, ctx);
 	}
	*/
	fq_poly_zero(Psi[0], ctx); //Psi_0 = 0
	fq_poly_one(Psi[1], ctx); //Psi_1 = 1
	fq_set_ui(tmp, 2, ctx);
	fq_mul(tmp, P.y, tmp, ctx);
	fq_poly_set_fq(Psi[2], tmp, ctx); //Psi_2 = 2 
	
	//Psi_3 = 3x^4+6ax^2+12bx-a^2
	fq_set_ui(tmp1, 3, ctx);
	fq_poly_set_coeff(Psi[3], 4, tmp1, ctx);//3*x^4

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
	fq_poly_scalar_mul_fq(Psi[4], Psi[4], P.y, ctx);

	fq_clear(tmp, ctx);fq_clear(tmp1, ctx);
	//Récurrence
	for(i=5; i<k+1; i++){
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
			fq_mul(tmp, tmp, P.y, ctx);
			fq_poly_scalar_mul_fq(Psi[i], Psi[m], tmp, ctx);
			fq_poly_mul(Psi[i], Psi[i], poly, ctx);
		}
		//Cas i impair: i=2*m+1
		else{
			m = i-1;
			m>>=1;
			fq_poly_sqr(y2, y2, ctx);
			if(m&1){
				fq_poly_pow(Psi[i], Psi[m], 3, ctx);
				fq_poly_mul(Psi[i], Psi[i], Psi[m+2], ctx);
				
				fq_poly_pow(poly, Psi[m+1], 3, ctx);
				fq_poly_mul(poly, poly, Psi[m-1], ctx);
				fq_poly_mul(poly, poly, y2, ctx);
				fq_poly_sub(Psi[i], Psi[i], poly, ctx);
			}
			else{
				fq_poly_pow(Psi[i], Psi[m], 3, ctx);
				fq_poly_mul(Psi[i], Psi[i], Psi[m+2], ctx);
				fq_poly_mul(Psi[i], Psi[i], y2, ctx);

				fq_poly_pow(poly, Psi[m+1], 3, ctx);
				fq_poly_mul(poly, poly, Psi[m-1], ctx);
				fq_poly_sub(Psi[i], Psi[i], poly, ctx);
			}
		}
	}
	fq_clear(tmp1, ctx); fq_clear(tmp, ctx);
	fq_poly_clear(poly, ctx); fq_poly_clear(y2, ctx);
	fq_ctx_clear(ctx);
}


int main(){
	fq_poly_t *Psi; 
	elliptic_curve E;
	signed long k = 13;
	const char* var = "x";
	fmpz_t p; // Premier
	fq_ctx_t ctx; // F_p
	fmpz_set_ui(p, 101);
	fq_ctx_init(ctx, p, 1, var);
	
	point P;
	fq_init_point(P, ctx);
	fq_set_ui(P.x, 1, ctx);
	fq_set_ui(P.y, 0, ctx);// P=(1,0)
	
	fq_init_curve(E, ctx);
	fq_set_ui(E.a, 100, ctx);
	fq_set_ui(E.b, 0, ctx);
	slong i;
	Psi = malloc((k + 1) * sizeof(fq_poly_t));
	div_poly(Psi, E, k, ctx, P);
	for(i = 0; i < k+1; i++){
		flint_printf("Psi[%lu] =",i);fq_poly_print_pretty(Psi[i], var, ctx) ;flint_printf("\n\n");
	}
	fmpz_clear(p);
	fq_ctx_clear(ctx);
	for(i = 0; i < k + 1; i++){
 		fq_poly_clear(Psi[i], ctx);
	}	
	free(Psi);
	return 0;
}
