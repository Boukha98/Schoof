#include<stdio.h>
#include<stdlib.h>
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fq.h>
#include<gmp.h>

typedef struct{
	fq_t x;
	fq_t y;
}point;

int main(){
	// Déclaration des variables
	point P, Q; // Points P et Q
    	fq_t m, m1, a, tmp1, tmp2;
	fmpz_t p; // Premier
	fmpz_init(p); fmpz_set_ui(p, 5);
	slong d = 2;
	const char* var = "0";
	
	fq_ctx_t ctx;
	fq_ctx_init(ctx, p, d, var);
	
	// Initialisation des variables
	fq_init(P.x, ctx); fq_init(P.y, ctx);
	fq_init(Q.x, ctx); fq_init(Q.y, ctx);
	fq_init(m, ctx); fq_init(m1, ctx);
	fq_init(tmp1, ctx); fq_init(tmp2, ctx);
    fq_init(a, ctx);
	
	// Deux points de la courbe y^2 = x^3 + a*x + b
	fq_set_ui(P.x, 1, ctx); fq_set_ui(P.y, 2, ctx);
	fq_set_ui(Q.x, 2, ctx); fq_set_ui(Q.y, 1, ctx);       
	fq_set_ui(a, 2, ctx);
	
	// P == Q
	// m = (3*(P.x)^2 + a)/2*P.y
	if((P.x == Q.x)&(P.y == Q.y)){
		fq_sqr(tmp1, P.x, ctx);
		fq_mul_ui(tmp1, tmp1, 3, ctx);
		fq_add(tmp1, tmp1, a, ctx);// 
		fq_mul_ui(tmp2, P.y, 2, ctx);
		fq_div(m, tmp1, tmp2, ctx);
	}
	// P \neq Q
	// m = (Q.y - P.y)/(Q.x - P.x)
	else{
		fq_sub(tmp1, Q.y, P.y, ctx);
		fq_sub(tmp2, Q.x, P.x, ctx);
		fq_div(m, tmp1, tmp2, ctx);
	}
	// res.x = m^2 - P.x - Q.x
	fq_sqr(m1, m, ctx);
	fq_add(tmp1, P.x, Q.x, ctx); 
	fq_sub(tmp1, m1, tmp1, ctx);
	

	fq_sub(tmp2, P.x, tmp1, ctx);
	fq_mul(tmp2, tmp2, m, ctx);
	fq_sub(tmp2, tmp2, P.y, ctx);
	
	// P + Q = (tmp1, tmp2)
	fq_print(tmp1,ctx); flint_printf("\n");
	fq_print(tmp2,ctx); flint_printf("\n");
	
	// Libération de la mémoire
	fq_clear(P.x, ctx); fq_clear(P.y, ctx);
	fq_clear(Q.x, ctx); fq_clear(Q.y, ctx);
	fq_clear(tmp1, ctx); fq_clear(tmp2, ctx);
	fq_clear(m, ctx); fq_clear(m1, ctx);
	fq_clear(a, ctx);
	fq_ctx_clear(ctx);
	fmpz_clear(p);
	return 0;
}
	
	
