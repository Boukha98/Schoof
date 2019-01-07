#include<stdio.h>
#include<stdlib.h>
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fq.h>
#include<gmp.h>

int main(){

	fq_t X,Y;
	fmpz_t p;	
	fmpz_init(p);fmpz_set_ui(p, 5);
	slong d = 2;
	const char* var = "0";
	
	fq_ctx_t ctx;

	fq_ctx_init(ctx, p, d, var);

	fq_init(X,ctx);fq_init(Y,ctx);

	fq_set_ui(X,1,ctx);fq_set_ui(Y,1,ctx);

	fq_neg(Y,Y,ctx);

	fq_print(X,ctx); flint_printf("\n");
	fq_print(Y,ctx); flint_printf("\n");

	fq_clear(X, ctx);
	fq_clear(Y, ctx);
	fq_ctx_clear(ctx);
	fmpz_clear(p);
	return 0;
}
