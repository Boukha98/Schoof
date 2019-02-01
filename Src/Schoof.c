#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fq.h>
#include<gmp.h>
#include<flint/fq_poly.h>

void fnext_prime(fmpz_t l){
	fmpz_add_ui(l, l, 2); // res = l+2
	while(!fmpz_is_prime(l)){ // 0 ou -1
		fmpz_add_ui(l,l,2); // res = res+2
	}
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

void Psi(fmpz_t l,fq_poly_t Psi,fq_ctx_t ctx){
	
}

void schoof(elliptic_curve E, fq_ctx_t ctx, fmpz_t q){
	ulong e;
	fq_t tmp1;
	fmpz_t l, a, q_sqrt,tmp,t;
	fq_poly_t poly; fq_poly_init(poly, ctx);fq_poly_t poly1; 
	fq_poly_init(poly1, ctx);fq_poly_t poly2; 
	fq_poly_init(poly2, ctx);
	fq_poly_t Psi1,Psi2;
	fq_poly_init(Psi1, ctx);fq_poly_init(Psi2, ctx);
	fq_init(tmp1,ctx);
	fmpz_init(l);fmpz_init(a);fmpz_init(q_sqrt);
	fmpz_init(tmp);fmpz_init(t);
	fq_init_curve(E, ctx);
	fmpz_set_ui(l, 2);
	
	fmpz_sqrt(q_sqrt, q);fmpz_mul_ui(q_sqrt, q_sqrt, 4);// borne 4*sqrt(p)
	
	//Construction du poly de l'end de frobenius Phi1 = x^q - x
	fq_poly_t Phi1;
	fq_poly_init(Phi1, ctx);
	fq_set_ui(tmp1, 1, ctx);
	fq_poly_set_coeff(Phi1, fmpz_get_si(q), tmp1, ctx);
	fq_neg(tmp1, tmp1, ctx);
	fq_poly_set_coeff(Phi1, 1, tmp1, ctx);
	
	//Construction de Phi2 = x^q² - x
	fq_poly_t Phi2;
	fq_poly_init(Phi2, ctx);
	fq_set_ui(tmp1, 1, ctx);
	fq_poly_set_coeff(Phi2, fmpz_get_si(q)*fmpz_get_si(q), tmp1, ctx);
	fq_neg(tmp1, tmp1, ctx);
	fq_poly_set_coeff(Phi2, 1, tmp1, ctx);
	
	//Poly y2
	fq_poly_t y2;
	fq_poly_init(y2, ctx);
	fq_set_ui(tmp1, 1, ctx);
	fq_poly_set_coeff(y2, 3, tmp1, ctx);//x^3
	fq_poly_set_coeff(y2, 1, E.a, ctx);//a*x
	fq_poly_set_coeff(y2, 0, E.b, ctx);//b
	
	//pgcd(Phi1, y2)
	fq_poly_gcd(poly, Phi1, y2, ctx);
	
	fmpz_set_ui(l,2);fmpz_set_ui(a,2);
	fmpz_sqrt(q_sqrt,q);fmpz_mul_ui(q_sqrt,q_sqrt,4);
	
	//Pour l=2
	if(fq_poly_is_one(poly,ctx)){
		fmpz_set_ui(t, 1);
	}
	else{fmpz_set_ui(t, 0);}
	
	while(fmpz_cmp(a,q_sqrt)<=0){
		fnext_prime(l);
		
    		//q_bar = q mod l
		fmpz_t q_bar; fmpz_init(q_bar);
		fmpz_mod(q_bar, q, l);
		fmpz_divexact_ui(tmp, l, 2); //tmp = l/2
		if(fmpz_cmp(q_bar, tmp)>0){
			fmpz_sub(q_bar, q_bar, l);//q_bar = q_bar - l
		}
		
		//Cas q_bar impair
	    	if(fmpz_is_odd(q_bar)){
			Psi(q_bar,Psi1,ctx);fq_poly_sqr(poly, Psi1, ctx);
		    	fq_poly_mul(poly, poly, Phi2, ctx);fmpz_sub_ui(tmp,q_bar,1);
			Psi(tmp,Psi1,ctx);fmpz_add_ui(tmp,q_bar,1);
			Psi(tmp,Psi2,ctx);
		    	fq_poly_mul(poly1, Psi1, Psi2, ctx);
		    	fq_poly_mul(poly1, poly1, y2, ctx);
		    	fq_poly_add(poly, poly, poly1, ctx);
	   	}

	    	//Cas q_bar pair
	    	else{
		    	Psi(q_bar,Psi1,ctx);fq_poly_sqr(poly, Psi1, ctx);
		    	fq_poly_mul(poly, poly, Phi2, ctx);
		    	fq_poly_mul(poly, poly, y2, ctx);
			fmpz_sub_ui(tmp,q_bar,1);
			Psi(tmp,Psi1,ctx);fmpz_add_ui(tmp,q_bar,1);
			Psi(tmp,Psi2,ctx);
		    	fq_poly_mul(poly1, Psi1, Psi2, ctx);
		    	fq_poly_add(poly, poly, poly1, ctx);
		}
		
		Psi(l,Psi1,ctx);
		fq_poly_gcd(poly, poly, Psi1, ctx);
		
		//Il existe un point de l-torsion P tq Phi^2(P)=+-[q_bar]P
		if(!fq_poly_is_one(poly,ctx)){
			//Cas 1: Phi^2(P)=-[q_bar]P
			if(fmpz_jacobi(q, l) == -1){//q n'est pas un résidu quadratique
				fmpz_set_ui(tmp, 0); // tmp = 0
				fmpz_CRT(t, t, a, tmp, l, 1); // TRC: t=0 mais on fait le calcul
			}

			//Cas 2: Phi^2(P)=[q_bar]P
			else{
				fmpz_t w;fmpz_init(w);
				fmpz_sqrtmod(w, q, l);

				//Cas w impair
				if(fmpz_is_odd(w)){
					Psi(w,Psi1,ctx);
			    		fq_poly_sqr(poly, Psi1, ctx);
			    		fq_poly_mul(poly, poly, Phi1, ctx);
					fmpz_sub_ui(tmp,w,1);
					Psi(tmp,Psi1,ctx);fmpz_add_ui(tmp,w,1);
					Psi(tmp,Psi2,ctx);
			    		fq_poly_mul(poly1, Psi1, Psi2, ctx);
			    		fq_poly_mul(poly1, poly1, y2, ctx);
			    		fq_poly_add(poly, poly, poly1, ctx);
			    	}

			    	//Cas w pair
			    	else{
			    		Psi(w,Psi1,ctx);fq_poly_sqr(poly, Psi1, ctx);
			    		fq_poly_mul(poly, poly, Phi1, ctx);
			    		fq_poly_mul(poly, poly, y2, ctx);
					fmpz_sub_ui(tmp,w,1);
					Psi(tmp,Psi1,ctx);fmpz_add_ui(tmp,w,1);
					Psi(tmp,Psi2,ctx);
			    		fq_poly_mul(poly1, Psi1, Psi2, ctx);
			    		fq_poly_add(poly, poly, poly1, ctx);
			    	}

			    	//Calcul du pgcd
				Psi(l,Psi1,ctx);
				fq_poly_gcd(poly, poly, Psi1, ctx);
				
				//Cas pgcd=1 cad w et -w ne sont pas solution: t=0
				if(fq_poly_is_one(poly, ctx)){
					fmpz_set_ui(t, 0);
				}

				//Cas pgcd \neq 1
				else{
					//Cas w impair
					if(fmpz_is_odd(w)){
				    		fmpz_add_ui(tmp, q, 3);
			    			fmpz_divexact_si(tmp, tmp, 2);
			    		}
			    		//Cas w pair
			    		else{
			    			fmpz_sub_ui(tmp, q, 1);
			    			fmpz_divexact_si(tmp, tmp, 2);
			    		}
			    		e=fmpz_get_ui(tmp);
			    		fq_poly_pow(poly, y2, e, ctx);
					Psi(w,Psi1,ctx);
			    		fq_poly_pow(poly1, Psi1, 3, ctx);
			    		fq_poly_mul(poly, poly, poly1, ctx);
					fq_set_ui(tmp1,4,ctx);
			    		fq_poly_scalar_mul_fq(poly, poly, tmp1, ctx);
					fmpz_add_ui(tmp,w,2);Psi(tmp,Psi1,ctx);
			    		fq_poly_pow(poly1, Psi1, 2, ctx);
					fmpz_sub_ui(tmp,w,1);Psi(tmp,Psi1,ctx);
			    		fq_poly_mul(poly1, poly1, Psi1, ctx);	
			    		fq_poly_sub(poly, poly, poly1, ctx);
					fmpz_sub_ui(tmp,w,2);Psi(tmp,Psi1,ctx);
			    		fq_poly_pow(poly1, Psi1, 2, ctx);
					fmpz_add_ui(tmp,w,1);Psi(tmp,Psi1,ctx);
			    		fq_poly_mul(poly1, poly1, Psi1, ctx);
			    		fq_poly_add(poly, poly, poly1, ctx);
			    				
			    		//Calcul du pgcd
					Psi(l,Psi1,ctx);
				    	fq_poly_gcd(poly, poly, Psi1, ctx);
				    		
			    		fmpz_mul_ui(tmp, w, 2);//tmp=2*w

			    		//pgcd = 1 cad t n'est pas 2w[l]: t = -2w[l]
			    		if(fq_poly_is_one(poly, ctx)){
				    		fmpz_negmod(tmp, tmp, l);
				    	}

			      		fmpz_CRT(t, t, a, tmp, l, 1);//TRC
				    	}
				    }
			    }
		}
		//Il n'existe pas point de l-torsion tq Phi^2(P)=+-[q_bar]P
		else{
			fq_poly_t beta2, alpha;
			fq_poly_init(beta2, ctx);fq_poly_init(alpha, ctx);

			//Calcul de alpha
			fmpz_sub_ui(tmp,q_bar,1);Psi(tmp,Psi1,ctx);
			fq_poly_sqr(poly, Psi1, ctx);
			fmpz_add_ui(tmp,q_bar,2);Psi(tmp,Psi1,ctx);
			fq_poly_mul(poly, poly, Psi1, ctx);

			fmpz_add_ui(tmp,q_bar,1);Psi(tmp,Psi1,ctx);
			fq_poly_sqr(poly1, Psi1, ctx);
			fmpz_sub_ui(tmp,q_bar,1);Psi(tmp,Psi1,ctx);
			fq_poly_mul(poly1, poly, Psi1, ctx);

			fq_poly_sub(poly, poly, poly1, ctx);

			Psi(q_bar,Psi1,ctx);
			fq_poly_pow(poly1, Psi1, 3, ctx);

			fmpz_pow_ui(tmp, q, 2);
			fmpz_add_ui(tmp, tmp, 1);
			fmpz_divexact_ui(tmp, tmp, 2);

			fq_poly_pow(alpha, y2, fmpz_get_si(tmp), ctx);
			fq_poly_mul(alpha, alpha, poly1, ctx);
			fq_set_si(tmp1, 4, ctx);
			fq_poly_scalar_mul_fq(alpha, alpha, tmp1, ctx);

			fq_poly_sub(alpha, poly, alpha, ctx);

			//Construction de X^q² + X^q + X = Phi3
			fq_poly_t Phi3;
			fq_poly_init(Phi3, ctx);
			fq_set_ui(tmp1, 1, ctx);
			fq_poly_set_coeff(Phi3, fmpz_get_si(q)*fmpz_get_si(q), tmp1, ctx);
			fq_poly_set_coeff(Phi3, fmpz_get_si(q), tmp1, ctx);
			fq_poly_set_coeff(Phi3, 1, tmp1, ctx);

			//Calcul de beta^2
			fmpz_sub_ui(tmp,q_bar,1);Psi(tmp,Psi1,ctx);
			fmpz_add_ui(tmp,q_bar,1);Psi(tmp,Psi2,ctx);
			fq_poly_mul(poly, Psi1, Psi2, ctx);
			Psi(q_bar,Psi1,ctx);
			fq_poly_sqr(poly1, Psi1, ctx);
			fq_poly_mul(poly1, poly1, Phi2, ctx);
			fq_poly_add(poly, poly, poly1, ctx);
			fq_poly_sqr(poly, poly, ctx);

			Psi(q_bar,Psi1,ctx);
			fq_poly_sqr(poly1, Psi1, ctx);
			fq_poly_mul(poly1, poly1, y2, ctx);
			fq_poly_mul(poly, poly, poly1, ctx);
			fq_set_ui(tmp1, 16, ctx);
			fq_poly_scalar_mul_fq(beta2, poly, tmp1, ctx);//beta^2

			//Calcul de gamma1 et teta1
			fq_poly_t gamma1,teta1;
			fq_poly_init(gamma1,ctx);fq_poly_init(teta1,ctx);
			fmpz_sub_ui(tmp,q_bar,1);Psi(tmp,Psi1,ctx);
			fmpz_add_ui(tmp,q_bar,1);Psi(tmp,Psi2,ctx);
			fq_poly_mul(poly,Psi1,Psi2,ctx);
			Psi(q_bar,Psi1,ctx);
			fq_poly_mul(poly1,Psi1,Phi3,ctx);
			fq_poly_sub(poly,poly,poly1,ctx);
			fq_poly_mul(poly,poly,beta2,ctx);
			Psi(q_bar,Psi1,ctx);
			fq_poly_sqr(poly1,Psi1,ctx);
			fq_poly_sqr(poly2,alpha,ctx);
			fq_poly_mul(poly1,poly1,poly2,ctx);
			fq_poly_add(gamma1,poly,poly1,ctx);//gamma1

			Psi(q_bar,Psi1,ctx);
			fq_poly_sqr(poly,Psi1,ctx);
			fq_poly_mul(teta1,poly,beta2,ctx);//teta1			
			
			//Calcul de gamma2 et teta2
			fq_poly_t gamma2,teta2;
			fq_poly_init(gamma2,ctx);fq_poly_init(teta2,ctx);	
			fq_set_si(tmp1,2,ctx);fq_poly_set_coeff(poly,fmpz_get_si(q)*fmpz_get_si(q),tmp1,ctx);
			fq_set_si(tmp1,1,ctx);fq_poly_set_coeff(poly,1,tmp1,ctx);
			Psi(q_bar,Psi1,ctx);
			fq_poly_sqr(poly1,Psi1,ctx);fq_poly_mul(poly,poly,poly1,ctx);
			fmpz_sub_ui(tmp,q_bar,1);Psi(tmp,Psi1,ctx);
			fmpz_add_ui(tmp,q_bar,1);Psi(tmp,Psi2,ctx);
			fq_poly_mul(poly2,Psi1,Psi2,ctx);fq_poly_sub(poly,poly,poly2,ctx);
			fq_poly_mul(poly,poly,alpha,ctx);

			fq_poly_neg(poly1,Phi2,ctx);Psi(q_bar,Psi1,ctx);
			fq_poly_sqr(poly2,Psi1,ctx);
			fq_poly_mul(poly1,poly1,poly2,ctx);
			fmpz_sub_ui(tmp,q_bar,1);Psi(tmp,Psi1,ctx);
			fmpz_add_ui(tmp,q_bar,1);Psi(tmp,Psi2,ctx);
			fq_poly_mul(poly2,Psi1,Psi2,ctx);
			fq_poly_sub(poly1,poly1,poly2,ctx);
			Psi(q_bar,Psi1,ctx);
			fq_poly_pow(poly2,Psi1,3,ctx);
			fq_poly_mul(poly1,poly1,poly2,ctx);fq_set_ui(tmp1,4,ctx);
			fq_poly_scalar_mul_fq(teta2,poly1,tmp1,ctx);//teta2

			fmpz_mul(tmp,q,q);fmpz_add_ui(tmp,tmp,1);
			fmpz_divexact_ui(tmp,tmp,2);e=fmpz_get_ui(tmp);
			fq_poly_pow(poly1,y2,e,ctx);fq_poly_mul(poly1,poly1,teta2,ctx);
			fq_poly_sub(poly,poly,poly1,ctx);fmpz_sub_ui(tmp,q,1);
			fmpz_divexact_ui(tmp,tmp,2);e=fmpz_get_ui(tmp);
			fq_poly_pow(poly1,y2,e,ctx);fq_set_ui(tmp1,4,ctx);
			fq_poly_scalar_mul_fq(poly1,poly1,tmp1,ctx);fq_poly_mul(gamma2,poly,poly1,ctx);//gamma2
			
			fmpz_t j; fmpz_init(j);fmpz_one(j);
			while(j<l){
				fq_poly_t absc,ordo;
				fq_poly_init(absc,ctx);fq_poly_init(ordo,ctx);

				//Verification abcisse
				fmpz_mul_ui(tmp,q,2);
				e=fmpz_get_ui(tmp);
				Psi(j,Psi1,ctx);
				fq_poly_pow(poly,Psi1,e,ctx);
				fq_poly_mul(poly,poly,gamma1,ctx);
				
				e=fmpz_get_ui(q);
				fmpz_sub_ui(tmp,j,1);Psi(tmp,Psi1,ctx);				
				fq_poly_pow(poly1,Psi1,e,ctx);
				fmpz_add_ui(tmp,j,1);Psi(tmp,Psi1,ctx);
				fq_poly_pow(poly2,Psi1,e,ctx);
				fq_poly_mul(poly1,poly1,poly2,ctx);
				fq_poly_mul(poly1,poly1,teta1,ctx);

				fq_poly_add(absc,poly,poly1,ctx);
				
				//Verification ordonnée
				e=fmpz_get_ui(q);e*=3;
				Psi(j,Psi1,ctx);
				fq_poly_pow(poly,Psi1,e,ctx);fq_poly_mul(poly,poly,gamma2,ctx);

				fmpz_sub_ui(tmp,j,1);Psi(tmp,Psi1,ctx);
				fq_poly_sqr(poly1,Psi1,ctx);
				fmpz_add_ui(tmp,j,2);Psi(tmp,Psi1,ctx);
				fq_poly_mul(poly1,poly1,Psi1,ctx);
				fmpz_add_ui(tmp,j,1);Psi(tmp,Psi1,ctx);
				fq_poly_sqr(poly2,Psi1,ctx);
				fmpz_sub_ui(tmp,j,2);Psi(tmp,Psi1,ctx);
				fq_poly_mul(poly2,poly2,Psi1,ctx);
				fq_poly_sub(poly1,poly1,poly2,ctx);e=fmpz_get_ui(q);
				fq_poly_pow(poly1,poly1,e,ctx);fq_poly_mul(poly1,poly1,teta2,ctx);

				fq_poly_sub(ordo,poly,poly1,ctx);

				fmpz_add_ui(j,j,1);
			}
		fmpz_mul(a,a,l);
	}

}

int main(){
	
}
