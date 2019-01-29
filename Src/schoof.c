#include<stdio.h>
#include<stdlib.h>
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fq.h>
#include<gmp.h>
#include<flint/fq_poly.h>

/*
    Entrées:
        1. Une courbe elliptique E = y^2 - x^3 - a*x - b
      	2. Un entier p pour le corps fini F_q avec q=p^d, d >= 1

    Sorite:
        le nombre de points de E sur F_q

    On choisi un ensemble de premiers impairs S qui ne contient pas p
	et tel que N = ∏l > 4*sqrt(q)

*/

typedef struct{}

void prime_set(fmpz * set, fmpz_t p){
	fmpz_t l_i, m, p_sqrt;
	fmpz_init(l_i);fmpz_init(m);fmpz_init(p_sqrt);
	fmpz_sqrt(p_sqrt, p);fmpz_mul_ui(p_sqrt, p_sqrt, 4);
	
}

void point_scalar_mul(point P, fmpz_t n, fq_ctx_t ctx, slong k){
	fq_init_point(P, ctx);
	fq_init_curve(E, ctx);
	fq_poly_t *Psi; fq_poly_t poly;
	slong i;
	for(i = 0; i < k + 1; i++){
 		fq_poly_init(Psi[i], ctx);
	}

	div_poly(Psi, E, k, ctx, P);

	fq_poly_mul(poly, Psi[n-1], Psi[n+1], ctx);


}
void schoof(elliptic_curve E, fq_ctx_t ctx, fmpz_t p, fmpz_t c){
	fmpz_t A, l, C, p_sqrt;
	fq_t tmp;
	fq_poly_t poly;
	fq_init(A, ctx);fq_init(l, ctx);fq_init(C, ctx);
	fmpz_init(p_sqrt);fq_init(t, ctx);
	fq_init_curve(E, ctx);
	fq_set_ui(A, 1, ctx);fq_set_ui(l, 3, ctx);fq_set_ui(C, 0, ctx);
	
	fmpz_sqrt(p_sqrt, p);fmpz_mul_ui(p_sqrt, p_sqrt, 4);// borne 4*sqrt(p)
	
	//Construction du poly de l'end de frobinius Phi1 = x^p - x
	fq_poly_t Phi1;fq_t t;
	fq_set_ui(tmp, 1, ctx);
	fq_poly_set_coeff(Phi, p, tmp, ctx);
	fq_neg(tmp, tmp, ctx);
	fq_poly_set_coeff(Phi, 1, tmp, ctx);
	
	//Poly y2
	fq_poly_init(y2, ctx);
	fq_set_ui(tmp, 1, ctx);
	fq_poly_set_coeff(y2, 3, tmp, ctx);//x^3
	fq_set(tmp, E.a, ctx);
	fq_poly_set_coeff(y2, 1, tmp, ctx);//a*x
	fq_set(tmp, E.b, ctx);
	fq_poly_set_coeff(y2, 0, tmp, ctx);//b
	
	//pgcd(Phi1, y2)
	fq_poly_gcd(poly, Phi1, y2, ctx);
	
	//si pgcd=1 alors t=1, sinon t=0
	if(fq_poly_is_one(poly)){
		fq_set_ui(t, 1, ctx);
	}
	else{fq_set_ui(t, 0, ctx);}
	
	fq_poly_t *Psi; //poly de division
    div_poly(Psi, E, k, ctx, P);//remplissage du tableau 

	

	while(fmpz_cmp(A, p_sqrt) < 0){}
			
	
}	
