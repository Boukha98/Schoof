#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<flint/flint.h>
#include<flint/fq.h>
#include<gmp.h>

void Elliptic_inv(fq_t X, fq_t Y, fq_t X1, fq_t Y1, fq_ctx_t ctx){
	fq_set(X1,X,ctx);
	fq_neg(Y1,Y,ctx);
}

void Poly_division(int m, fq_t a, fq_t b, fq_t X, fq_t Y, fq_ctx_t ctx){
	fq_t X2,X4,X6,a2,a3,b2,c1,c2,c3,c4,c5;
	int n=floor(m/2); int i=0;
	fq_t Psi[m+1];
	
	fq_init(X2,ctx);fq_init(X4,ctx); fq_init(X6,ctx);fq_init(a2,ctx); 
	fq_init(a3, ctx);fq_init(b2,ctx);fq_init(c1,ctx);fq_init(c2,ctx);
	fq_init(c3,ctx);fq_init(c4,ctx);fq_init(c5,ctx); 
	for(i=0;i<=m+1;i++){fq_init(Psi[i],ctx);}
	
	fq_sqr(X2,X,ctx);fq_sqr(X4,X2,ctx); fq_mul(X6,X4,X2,ctx);fq_sqr(a2,a,ctx);
	fq_mul(a3,a,a2,ctx);fq_sqr(b2,b,ctx);
	
	fq_set_si(Psi[0],0,ctx);fq_set_si(Psi[1],1,ctx);
	
	fq_mul_si(Psi[2],Y,2,ctx);
	
	fq_mul_si(c1,X4,3,ctx);fq_mul(c2,a,X2,ctx);fq_mul_si(c2,c2,6,ctx);
	fq_mul(c3,b,X,ctx);fq_mul_si(c3,c3,12,ctx);fq_set(c4,a2,ctx);
	fq_add(c1,c1,c2,ctx); fq_sub(c3,c3,c4,ctx);fq_add(Psi[3],c1,c3,ctx);
	
	fq_mul(c1,a,X4,ctx);fq_mul_si(c1,c1,5,ctx);fq_mul(c2,b,X2,ctx);
	fq_mul_si(c2,c2,20,ctx);fq_mul(c3,a2,X2,ctx);fq_mul_si(c3,c3,5,ctx);
	fq_mul(c4,a,b,ctx);fq_mul(c4,c4,X,ctx);fq_mul_si(c4,c4,4,ctx);
	fq_mul_si(c5,b2,8,ctx);	fq_sub(c3,X6,c3,ctx);fq_sub(c1,c1,c4,ctx);
	fq_sub(c2,c2,c5,ctx);fq_sub(c2,c2,a3,ctx);fq_add(c1,c1,c2,ctx);
	fq_add(c1,c1,c3,ctx);fq_mul(c1,c1,Y,ctx);fq_mul_si(Psi[4],c1,4,ctx);
	
	fq_pow_ui(c1,Psi[2],3,ctx);fq_mul(c2,c1,Psi[4],ctx);fq_pow_ui(c3,Psi[3],3,ctx);
	fq_sub(Psi[5],c2,c3,ctx);
	
	if(n>=3){
		for(i=3;i<=n;i++){
			fq_sqr(c1,Psi[n-1],ctx);fq_mul(c2,c1,Psi[n+2],ctx);fq_sqr(c3,Psi[n+1],ctx);
			fq_mul(c4,c3,Psi[n-2],ctx);fq_sub(c2,c2,c4,ctx);fq_mul(c2,c2,Psi[n],ctx);
			fq_mul_si(c4,Y,2,ctx);fq_div(Psi[2*n],c2,c4,ctx);
			
			fq_pow_ui(c1,Psi[n],3,ctx);fq_mul(c2,c1,Psi[n+2],ctx);fq_pow_ui(c3,Psi[n+1],3,ctx);
			fq_mul(c4,Psi[n-1],c3,ctx);fq_sub(Psi[2*n+1],c2,c4,ctx);
		}
	}
	
	fq_print(Psi[m],ctx); flint_printf("\n");

	fq_clear(X2,ctx);fq_clear(X4,ctx);fq_clear(X6,ctx);fq_clear(a2,ctx);
	fq_clear(a3,ctx);fq_clear(b2,ctx);fq_clear(c1,ctx);fq_clear(c2,ctx);
	fq_clear(c3,ctx);fq_clear(c4,ctx);fq_clear(c5,ctx);
	for(i=0;i<=m+1;i++){fq_clear(Psi[i],ctx);}
}

void main(int argc, char** argv){
	fmpz_t p;fmpz_init(p);fmpz_set_ui(p,5);
	signed long d = 1;
	char *var="";
	fq_ctx_t ctx;
	fq_ctx_init (ctx, p, d, var);
	
	fq_t X,Y,a,b;
	
	fq_init(X,ctx);fq_init(Y,ctx);fq_init(a,ctx);fq_init(b,ctx);
	
	fq_set_si(X,1,ctx);fq_set_si(Y,1,ctx);fq_set_si(a,-3,ctx);fq_set_si(b,3,ctx);
	
	Poly_division(8,a,b,X,Y,ctx);

	fq_clear(X,ctx);fq_clear(Y,ctx);fq_clear(a,ctx);fq_clear(b,ctx);
	fq_ctx_clear(ctx);fmpz_clear(p);
}
