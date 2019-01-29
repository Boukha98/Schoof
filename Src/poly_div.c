#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fq.h>
#include<gmp.h>
#include<flint/fq_poly.h>

<<<<<<< HEAD
void poly_div(int m, fq_t a, fq_t b, fq_t Y, fq_ctx_t ctx){
=======
void poly_div(int m, fq_t a, fq_t b, fq_t X, fq_t Y, fq_ctx_t ctx){
>>>>>>> 353c501e729e51e1be2166d4d029c8bc87b972fc
	int n=m>>1; int i=0; int k=5;
	if(m>5){k=m+1;}
	fq_poly_t Psi[k];
	fq_t tmp,tmp1;
	fq_poly_t tmp2,tmp3;
	
	//Initialisation des variables
	for(i=0;i<=k;i++){fq_poly_init(Psi[i],ctx);}
	fq_init(tmp,ctx);fq_init(tmp1,ctx);fq_poly_init(tmp2,ctx);
	fq_poly_init(tmp3,ctx);
	
	//Psi_0 = 0
	fq_poly_zero(Psi[0], ctx); 
	
	//Psi_1 = 1
	fq_poly_one(Psi[1], ctx); 
	
	//Psi_2 = 2 
	fq_set_ui(tmp, 2, ctx);
	fq_mul(tmp, Y, tmp, ctx);
	fq_poly_set_fq(Psi[2], tmp, ctx);
	
	//Psi_3 = 3x^4+6ax^2+12bx-a^2
	fq_set_si(tmp, 3, ctx);fq_poly_set_coeff(Psi[3], 4, tmp, ctx);//3*x^4
	fq_mul_ui(tmp, a, 6, ctx);fq_poly_set_coeff(Psi[3], 2, tmp, ctx);//6*a*x^2
	fq_mul_ui(tmp, b, 12, ctx);fq_poly_set_coeff(Psi[3], 1, tmp, ctx);//12*b*x
	fq_sqr(tmp, a, ctx);fq_neg(tmp, tmp, ctx);
	fq_poly_set_coeff(Psi[3], 0, tmp, ctx);//-a^2

	//Psi_4 = 4y(x^6+5ax^4+20bx^3-5a^2x^2-4abx-8b^2-a^3)
	fq_mul_si(tmp, tmp, 5, ctx);fq_neg(tmp, tmp, ctx);
	fq_poly_set_coeff(Psi[4], 2, tmp, ctx);//-5a^2x^2
	fq_set_si(tmp,1,ctx);fq_poly_set_coeff(Psi[4], 6, tmp, ctx);//x^6
	fq_mul_si(tmp,a,5,ctx);fq_poly_set_coeff(Psi[4], 4, tmp, ctx);//5ax^4
	fq_mul_si(tmp, b, 20, ctx);fq_poly_set_coeff(Psi[4], 3, tmp, ctx);//20bx^3
	fq_mul_si(tmp, a, -4, ctx);fq_mul(tmp, tmp, b, ctx);
	fq_poly_set_coeff(Psi[4], 1, tmp, ctx);//-4abx
	fq_sqr(tmp,b,ctx);fq_mul_si(tmp,tmp,-8,ctx);fq_pow_ui(tmp1,a,3,ctx);
	fq_sub(tmp,tmp,tmp1,ctx);fq_poly_set_coeff(Psi[4], 0, tmp, ctx);//-8b^2-a^3
	
	/* Psi5 = Psi4*(Psi2^3)-(Psi3^3)*Psi1 */
	fq_poly_pow(tmp2,Psi[2],3,ctx);fq_poly_mul(tmp2,tmp2,Psi[4],ctx);
	fq_poly_pow(tmp3,Psi[3],3,ctx);fq_poly_sub(Psi[5],tmp2,tmp3,ctx);
	
	/* Psim, si n>=3 */
	if(n>=3){
		for(i=3;i<=n;i++){
			if(fq_is_zero(Y,ctx)){fq_poly_zero(Psi[2*i],ctx);}
			else{
				fq_poly_sqr(tmp2,Psi[i-1],ctx);fq_poly_mul(tmp2,tmp2,Psi[i+2],ctx);
				fq_poly_sqr(tmp3,Psi[i+1],ctx);fq_poly_mul(tmp3,tmp3,Psi[i-2],ctx);
				fq_poly_sub(tmp2,tmp2,tmp3,ctx);fq_poly_mul(tmp2,tmp2,Psi[i],ctx);
				fq_mul_si(tmp,Y,2,ctx);fq_inv(tmp,tmp,ctx);
				fq_poly_scalar_mul_fq(Psi[2*i],tmp2,tmp,ctx);
			}
					
			fq_poly_pow(tmp2,Psi[i],3,ctx);fq_poly_mul(tmp2,tmp2,Psi[i+2],ctx);
			fq_poly_pow(tmp3,Psi[i+1],3,ctx);fq_poly_mul(tmp3,tmp3,Psi[i-1],ctx);
			fq_poly_sub(Psi[2*i+1],tmp2,tmp3,ctx);
		}
	}
	
	fq_poly_print_pretty(Psi[m],"x",ctx); flint_printf("\n");
	
	fq_clear(tmp,ctx);fq_clear(tmp1,ctx);fq_poly_clear(tmp2,ctx);
	fq_poly_clear(tmp3,ctx);for(i=0;i<=k;i++){fq_poly_clear(Psi[i],ctx);}
}

/*
nP
<<<<<<< HEAD
=======

>>>>>>> 353c501e729e51e1be2166d4d029c8bc87b972fc
fq_t X,Y,a,b,Xn,Yn,P1,P2,P3,P4,c1,c2,c3;
	
	fq_init(X,ctx);fq_init(Y,ctx);fq_init(a,ctx);fq_init(b,ctx);fq_init(Xn,ctx);
	fq_init(Yn,ctx);fq_init(P1,ctx);fq_init(P2,ctx);fq_init(P3,ctx);fq_init(P4,ctx);
	fq_init(c1,ctx);fq_init(c2,ctx);fq_init(c3,ctx);
	
	fq_set_si(X,1,ctx);fq_set_si(Y,1,ctx);fq_set_si(a,-3,ctx);fq_set_si(b,3,ctx);
	Poly_division(t-1,a,b,X,Y,P1,ctx);Poly_division(t,a,b,X,Y,P2,ctx);
	Poly_division(t+1,a,b,X,Y,P3,ctx);Poly_division(2*t,a,b,X,Y,P4,ctx);
	
	
	fq_mul(c1,P1,P3,ctx);fq_sqr(c2,P2,ctx);fq_div(c1,c1,c2,ctx);
	fq_sub(Xn,X,c1,ctx);
	
	fq_pow_ui(c3,P2,4,ctx);fq_mul_si(c3,c3,2,ctx);fq_div(Yn,P4,c3,ctx);
	
	fq_print(Xn,ctx); flint_printf("\n");
	fq_print(Yn,ctx); flint_printf("\n");
<<<<<<< HEAD
=======

>>>>>>> 353c501e729e51e1be2166d4d029c8bc87b972fc
	fq_clear(X,ctx);fq_clear(Y,ctx);fq_clear(a,ctx);fq_clear(b,ctx);
	fq_clear(Xn,ctx);fq_clear(Yn,ctx);fq_clear(P1,ctx);fq_clear(P2,ctx);
	fq_clear(P3,ctx);fq_clear(P4,ctx);fq_clear(c1,ctx);fq_clear(c2,ctx);
	fq_clear(c3,ctx);
<<<<<<< HEAD
=======

>>>>>>> 353c501e729e51e1be2166d4d029c8bc87b972fc
*/

typedef struct{
	fq_t x, y; // Coordonnées du point
}point;

//initialisation
void fq_init_point(point P, fq_ctx_t ctx){
	fq_init(P.x, ctx);
	fq_init(P.y, ctx);
}

void fq_clear_point(point P, fq_ctx_t ctx){
	fq_clear(P.x, ctx);
	fq_clear(P.y, ctx);
}

//Courbe elliptique y^2 = x^3+a*x+b mod p
typedef struct{
	fq_t a, b; 
}elliptic_curve;

//initialisation
void fq_init_curve(elliptic_curve E, fq_ctx_t ctx){
	fq_init(E.a, ctx);
	fq_init(E.b, ctx);
}

void fq_clear_curve(elliptic_curve E, fq_ctx_t ctx){
	fq_clear(E.a, ctx);
	fq_clear(E.b, ctx);
}

int main(){
	fmpz_t p;fmpz_init(p);fmpz_set_ui(p,5);
	signed long d = 1;
	char *var="x";
	fq_ctx_t ctx;
	fq_ctx_init (ctx, p, d, var);
	
<<<<<<< HEAD
	fq_t a,b,Y;fq_init(a,ctx);fq_init(b,ctx);fq_init(Y,ctx);
	
	fq_set_si(a,1,ctx);fq_set_si(b,0,ctx);
	fq_set_si(Y,0,ctx);
	
	poly_div(3,a,b,Y,ctx);
	
	fq_clear(a,ctx);fq_clear(b,ctx);fq_clear(Y,ctx);
=======
	fq_t a,b,X,Y;fq_init(a,ctx);fq_init(b,ctx);fq_init(X,ctx);fq_init(Y,ctx);
	
	fq_set_si(a,-3,ctx);fq_set_si(b,3,ctx);fq_set_si(X,1,ctx);
	fq_set_si(Y,1,ctx);
	
	poly_div(3,a,b,X,Y,ctx);
	
	fq_clear(a,ctx);fq_clear(b,ctx);fq_clear(X,ctx);fq_clear(Y,ctx);
>>>>>>> 353c501e729e51e1be2166d4d029c8bc87b972fc
	fq_ctx_clear(ctx);fmpz_clear(p);
	return 0;
}
