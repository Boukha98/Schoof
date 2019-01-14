#include<stdio.h>
#include<stdlib.h>
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fq.h>
#include<gmp.h>

//Point
typedef struct{
	fq_t x, y; // Coordonnées du point
}point;

//Courbe elliptique
typedef struct{
	fmpz_t p;
	fq_t a, b; // y^2 = x^3+a*x+b mod p
}elliptic_curve;

//Polynome
typedef struct{
	int deg; // Degré du polynome
	int *coef; // Tableau des coefficients
}poly;

// delta = -(4*a^3 + 27*b^2) mod p
void delta_ec(elliptic_curve E){
	fmpz_t delta;
	fmpz_powm_ui(E.b, E.b, 2, E.p);
	fmpz_mul_ui(E.b, E.b, 27);
	fmpz_mod(E.b, E.b, E.p);
	
	fmpz_powm_ui(E.a, E.a, 3, E.p);
	fmpz_mul_ui(E.a, E.a, 4);
	fmpz_mod(E.a, E.a, E.p);
	
	fmpz_add(delta, E.a, E.b);
	fmpz_neg(delta, delta);
	fmpz_mod(delta, delta, E.p);
}


