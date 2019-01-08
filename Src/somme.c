#include<stdio.h>
#include<stdlib.h>
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fq.h>
#include<gmp.h>

int main(){
	
	fq_t X, Y;
	fq_t X1, X2, Y1;
	fq_t m, m1, a;
	fmpz_t p;	
	fmpz_init(p);fmpz_set_ui(p, 5);
	slong d = 2;
	const char* var = "0";
	
	fq_ctx_t ctx;

	fq_ctx_init(ctx, p, d, var);

	fq_init(X,ctx);fq_init(Y,ctx);fq_init(m,ctx); fq_init(a,ctx);
	fq_init(X1, ctx);fq_init(X2, ctx);fq_init(Y1,ctx);fq_init(m1, ctx);

	fq_set_ui(X,1,ctx);fq_set_ui(Y,1,ctx); fq_set_ui(a,2,ctx);
	fq_set(X1,X,ctx);fq_set(X2,X,ctx);fq_set(Y1,Y,ctx);

	fq_sqr(X, X, ctx);
	fq_mul_ui(X, X, 3, ctx);
	fq_add(X, X, a, ctx);// 
	fq_mul_ui(Y, Y, 2, ctx);
	fq_div(m, X, Y, ctx);
	
	fq_sqr(m1, m, ctx);
	fq_mul_ui(X1, X1, 2, ctx); 
	fq_sub(X1, m1, X1, ctx);
	
	fq_mul_ui(X2, X2, 3, ctx); 
	fq_sub(X2, X2, m1, ctx);
	fq_mul(X2, X2, m, ctx);
	fq_sub(Y1, X2, Y1, ctx);
	
	fq_print(X1,ctx); flint_printf("\n");
	fq_print(Y1,ctx); flint_printf("\n");

	fq_clear(X, ctx);
	fq_clear(Y, ctx);
	fq_clear(X1, ctx);
	fq_clear(Y1, ctx);
	fq_clear(X2, ctx);
	fq_clear(m, ctx);
	fq_clear(m1, ctx);
	fq_clear(a, ctx);
	fq_ctx_clear(ctx);
	fmpz_clear(p);
	return 0;
}
	
	
