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

void prime_set(fmpz * set, fmpz_t p){
	fmpz_t l_i, m, p_sqrt;
	fmpz_init(l_i);fmpz_init(m);fmpz_init(p_sqrt);
	fmpz_sqrt(p_sqrt, p);fmpz_mul_ui(p_sqrt, p_sqrt, 4);
	

void schoof(elliptic_curve E, fq_ctx_t ctx, fmpz_t p){
	fmpz_t A, l, C, p_sqrt;
	fq_init(A, ctx);fq_init(l, ctx);fq_init(C, ctx);fmpz_init(p_sqrt);
	fq_init_curve(E, ctx);
	fq_set_ui(A, 1, ctx);fq_set_ui(l, 3, ctx);fq_set_ui(C, 0, ctx);
	
	fmpz_sqrt(p_sqrt, p);fmpz_mul_ui(p_sqrt, p_sqrt, 4);// borne 4*sqrt(p)
	
	while(fmpz_cmp(A, p_sqrt) < 0){
		
	
	
