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

void schoof(elliptic_curve E, fq_ctx_t ctx, fmpz_t q, fmpz_t c){
	fmpz_t l, a, q_sqrt,tmp, tmp1,t;
	fq_poly_t poly; fq_poly_init(poly, ctx);fq_poly_t poly1; 
	fq_poly_init(poly1, ctx);
	fmpz_init(l);fmpz_init(a);fmpz_init(q_sqrt);
	fmpz_init(tmp);fmpz_init(tmp1);fmpz_init(t);
	fq_init(t, ctx);fq_init(tmp, ctx);
	fq_init_curve(E, ctx);
	fq_set_ui(l, 2, ctx);
	
	fmpz_sqrt(q_sqrt, q);fmpz_mul_ui(q_sqrt, q_sqrt, 4);// borne 4*sqrt(p)
	
	//Construction du poly de l'end de frobinius Phi1 = x^q - x
	fq_poly_t Phi1;
	fq_poly_init(Phi1, ctx);
	fq_set_ui(tmp, 1, ctx);
	fq_poly_set_coeff(Phi, fmpz_get_si(q), tmp, ctx);
	fq_neg(tmp, tmp, ctx);
	fq_poly_set_coeff(Phi, 1, tmp, ctx);
	
	//Construction de Phi2 = x^q² - x
	fq_poly_t Phi2;
	fq_poly_init(Phi2, ctx);
	fq_set_ui(tmp, 1, ctx);
	fq_poly_set_coeff(Phi, fmpz_get_si(q)*fmpz_get_si(q), tmp, ctx);
	fq_neg(tmp, tmp, ctx);
	fq_poly_set_coeff(Phi, 1, tmp, ctx);
	
	//Poly y2
	fq_poly_t(y2, ctx);
	fq_poly_init(y2, ctx);
	fq_set_ui(tmp, 1, ctx);
	fq_poly_set_coeff(y2, 3, tmp, ctx);//x^3
	fq_poly_set_coeff(y2, 1, E.a,f ctx);//a*x
	fq_poly_set_coeff(y2, 0, E.b, ctx);//b
	
	//pgcd(Phi1, y2)
	fq_poly_gcd(poly, Phi1, y2, ctx);
	
	fmpz_set_ui(l,2);fmpz_set_ui(a,2);
	fmpz_sqrt(q_sqrt,q);fmpz_mul_ui(q_sqrt,q_sqrt,4);
	
	//Pour l=2
	if(fq_poly_is_one(poly)){
		fq_set_ui(t, 1, ctx);
	}
	else{fq_set_ui(t, 0, ctx);}
	
	while(fmpz_cmp(a,q_sqrt)<=0){
		fnext_prime(l);
    		//q_bar = q mod l
		fmpz_t q_bar; fmpz_init(q_bar);
		fmpz_mod(q_bar, q, l);
		fmpz_div_si(tmp, l, 2); //tmp = l/2
		if(fmpz_cmp(q_bar, tmp)>=0){
			fmpz_sub(q_bar, q_bar, l);//q_bar = q_bar - l
		}
		
		//Cas q_bar impair
	    	if(q_bar&1){
		    	fq_poly_sqr(poly, Psi[q_bar], ctx);
		    	fq_poly_mul(poly, poly, Phi2, ctx);
		    	fq_poly_mul(poly1, Psi[q_bar-1], Psi[q_bar+1], ctx);
		    	fq_poly_mul(poly1, poly1, y2, ctx);
		    	fq_poly_add(poly, poly, poly1, ctx);
	   	}

	    	//Cas q_bar pair
	    	else{
		    	fq_poly_sqr(poly, Psi[q_bar], ctx);
		    	fq_poly_mul(poly, poly, Phi2, ctx);
		    	fq_poly_mul(poly, poly, y2, ctx);
		    	fq_poly_mul(poly1, Psi[q_bar-1], Psi[q_bar+1], ctx);
		    	fq_poly_add(poly, poly, poly1, ctx);
		}
		
		fq_poly_gcd(poly, poly, Psi[l], ctx);
		
		//Il existe un point de l-torsion P tq Phi^2(P)=+-[q_bar]P
		if(!(fq_poly_is_one(poly,ctx)){
			//Cas 1: Phi^2(P)=-[q_bar]P
			if(fmpz_jacobi(q, l) == -1){//q n'est pas un résidu quadratique
				fmpz_set_si(tmp, 0); // tmp = 0
				fmpz_CRT(t, t, a, tmp, l, 1); // TRC: t=0 mais on fait le calcul
			}

			//Cas 2: Phi^2(P)=[q_bar]P
			else{
				fmpz_t w;fmpz_init(w);
				fmpz_sqrtmod(w, q, l);
				//Cas w impair
				if(w&1){
			    		fq_poly_sqr(poly, Psi[w], ctx);
			    		fq_poly_mul(poly, poly, Phi1, ctx);
			    		fq_poly_mul(poly1, Psi[w-1], Psi[w+1], ctx);
			    		fq_poly_mul(poly1, poly1, y2, ctx);
			    		fq_poly_add(poly, poly, poly1, ctx);
			    	}
			    	//Cas w pair
			    	else{
			    		fq_poly_sqr(poly, Psi[w], ctx);
			    		fq_poly_mul(poly, poly, Phi1, ctx);
			    		fq_poly_mul(poly, poly, y2, ctx);
			    		fq_poly_mul(poly1, Psi[w-1], Psi[w+1], ctx);
			    		fq_poly_add(poly, poly, poly1, ctx);
			    	}
			    	//Calcul du pgcd
				fq_poly_gcd(poly, poly, Psi[l], ctx);
				
				//Cas pgcd=1 cad w et -w ne sont pas solution: t=0
				if(fq_poly_is_one(poly, ctx)){
					fq_set_ui(t, 0, ctx);
				}

				//Cas pgcd \neq 1
				else{
					//Cas w impair
					if(w&1){
				    		fmpz_add_ui(tmp, q, 3);
			    			fmpz_divexact_si(tmp, tmp, 2);
			    		}
			    		//Cas w pair
			    		else{
			    			fmpz_sub_ui(tmp, q, 1);
			    			fmpz_divexact_si(tmp, tmp, 2);
			    		}
			    		tmp=fmpz_get_ui(tmp);
			    		fq_poly_pow(poly, y2, tmp, ctx);
			    		fq_poly_pow(poly1, Psi[w], 3, ctx);
			    		fq_poly_mul(poly, poly, poly1, ctx);
			    		fq_poly_scalar_mul_fq(poly, poly, 4, ctx);
			    		fq_poly_pow(poly1, Psi[w+2], 2, ctx);
			    		fq_poly_mul(poly1, poly1, Psi[w-1], ctx);	
			    		fq_poly_sub(poly, poly, poly1, ctx);
			    		fq_poly_pow(poly1, Psi[w-2], 2, ctx);
			    		fq_poly_mul(poly1, poly1, Psi[w+1], ctx);
			    		fq_poly_add(poly, poly, poly1, ctx);
			    				
			    		//Calcul du pgcd
				    	fq_poly_gcd(poly, poly, Psi[l], ctx);
				    		
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
			fq_poly_sqr(poly, Psi[q_bar-1], ctx);
			fq_poly_mul(poly, poly, Psi[q_bar+2], ctx);

			fq_poly_sqr(poly1, Psi[q_bar+1], ctx);
			fq_poly_mul(poly1, poly, Psi[q_bar-1], ctx);

			fq_poly_sub(poly, poly, poly1, ctx);

			fq_poly_pow(poly1, Psi[q_bar], 3, ctx);

			fmpz_pow_si(tmp, q_bar, 2);
			fmpz_add_si(tmp, tmp, 1);
			fmpz_divexact_si(tmp, tmp, 2);

			fq_poly_pow(alpha, y2, fmpz_get_si(tmp), ctx);
			fq_poly_mul(alpha, alpha, poly1, ctx);
			fq_set_si(tmp1, 4, ctx);
			fq_poly_scalar_mul(alpha, alpha, tmp1, ctx);

			fq_poly_sub(alpha, poly, alpha, ctx);

			//Construction de X^q² + X^q + X = Phi3
			fq_poly_t Phi3;
			fq_poly_init(Phi3, ctx);
			fq_set_ui(tmp, 1, ctx);
			fq_poly_set_coeff(Phi3, fmpz_get_si(q)*fmpz_get_si(q), tmp, ctx);
			fq_poly_set_coeff(Phi3, fmpz_get_si(q), tmp, ctx);
			fq_poly_set_coeff(Phi3, 1, tmp, ctx);

			//Calcul de beta^2
			fq_poly_mul(poly, Psi[q_bar-1], Psi[q_bar+1], ctx);
			fq_poly_sqr(poly1, Psi[q_bar], ctx);
			fq_poly_mul(poly1, poly1, Phi2, ctx);
			fq_poly_add(poly, poly, poly1, ctx);
			fq_poly_sqr(poly, poly, ctx);

			fq_poly_sqr(poly1, Psi[q_bar], ctx);
			fq_poly_mul(poly1, poly1, y2, ctx);
			fq_poly_mul(poly, poly, poly1, ctx);
			fq_set_ui(tmp1, 16, ctx);
			fq_poly_scalar_mul_fq(beta2, poly, tmp1, ctx);//beta^2

			//Calcul de gamma1 et teta1
			fq_poly_t gamma1,teta1;
			fq_poly_init(gamma1,ctx);fq_poly_init(teta1,ctx);
			fq_poly_mul(poly,Psi[q_bar-1],Psi[q_bar+1],ctx);
			fq_poly_mul(poly1,Psi[q_bar],Phi3,ctx);
			fq_poly_sub(poly,poly,poly1,ctx);
			fq_poly_mul(poly,poly,beta2,ctx);
			fq_poly_sqr(poly1,Psi[q_bar],ctx);
			fq_poly_sqr(poly2,alpha,ctx);
			fq_poly_mul(poly1,poly1,poly2,ctx);
			fq_poly_add(gamma1,poly,poly1,ctx);//gamma1

			fq_poly_sqr(poly,Psi[q_bar],ctx);
			fq_poly_mul(teta1,poly,beta2,ctx);//teta1			
			
			//Calcul de gamma2 et teta2
			fq_poly_t gamma2,teta2;
			fq_poly_init(gamma2,ctx);fq_poly_init(teta2,ctx);	
			fq_set_si(tmp,2,ctx);fq_poly_set_coeff(poly,fmpz_get_si(q)*fmpz_get_si(q),tmp,ctx);
			fq_set_si(tmp,1,ctx);fq_poly_set_coeff(poly,1,tmp,ctx);
			fq_poly_sqr(poly1,Psi[q_bar],ctx);fq_poly_mul(poly,poly,poly1,ctx);
			fq_poly_mul(poly2,Psi[q_bar-1],Psi[q_bar+1],ctx);fq_poly_sub(poly,poly,poly2,ctx);
			fq_poly_mul(poly,poly,alpha,ctx);

			fq_poly_neg(poly1,Phi2,ctx);fq_poly_sqr(poly2,Psi[q_bar],ctx);
			fq_poly_mul(poly1,poly1,poly2,ctx);fq_poly_mul(poly2,Psi[q_bar-1],Psi[q_bar+1],ctx);
			fq_poly_sub(poly1,poly1,poly2,ctx);fq_poly_pow(poly2,Psi[q_bar],3,ctx);
			fq_poly_mul(poly1,poly1,poly2,ctx);fq_set_ui(tmp,4,ctx);
			fq_poly_scalar_mul(teta2,poly1,tmp,ctx);fmpz_sqr(tmp,q);//teta2
			fmpz_add_ui(tmp,tmp,1);fmpz_divexact_ui(tmp,tmp,2);
			tmp=fmpz_get_ui(tmp);fq_poly_pow(poly1,y2,tmp,ctx);
			fq_poly_mul(poly1,poly1,teta2,ctx);

			fq_poly_sub(poly,poly,poly1,ctx);fmpz_sub_ui(tmp,q,1);
			fmpz_divexact_ui(tmp,tmp,2);tmp=fmpz_get_ui(tmp);
			fq_poly_pow(poly1,y2,tmp,ctx);fq_set_ui(tmp,4,ctx);
			fq_poly_scalar_mul(poly1,poly1,tmp,ctx);fq_poly_mul(gamma2,poly,poly1,ctx);//gamma2
			
			for(j=0,j<l,j++){
				fq_poly_t absc,ordo;
				fq_poly_init(absc,ctx);fq_poly_init(ordo,ctx);
				//Verification abcisse
				fmpz_mul_ui(tmp,q,2,ctx);
				tmp=fmp_get_ui(tmp);
				fq_poly_pow(poly,Psi[j],tmp,ctx);
				fq_poly_mul(poly,poly,gamma1,ctx);
				
				tmp=fmp_get_ui(q);				
				fq_poly_pow(poly1,Psi[j-1],tmp,ctx);
				fq_poly_pow(poly2,Psi[j+1],tmp,ctx);
				fq_poly_mul(poly1,poly1,poly2,ctx);
				fq_poly_mul(poly1,poly1,teta1,ctx);

				fq_poly_add(absc,poly,poly1,ctx);
				
				//Verification ordonnée
				tmp=fmpz_get_ui(q);tmp*=3;
				fq_poly_pow(poly,Psi[j],tmp,ctx);fq_poly_mul(poly,poly,gamma2,ctx);

				fq_poly_sqr(poly1,Psi[j-1],ctx);fq_poly_mul(poly1,poly1,Psi[j+2],ctx);
				fq_poly_sqr(poly2,Psi[j+1],ctx);fq_poly_mul(poly2,poly2,Psi[j-2],ctx);
				fq_poly_sub(poly1,poly1,poly2,ctx);tmp=fmpz_get_ui(q);
				fq_poly_pow(poly1,poly1,tmp,ctx);fq_poly_mul(poly1,poly1,teta2,ctx);

				fq_poly_sub(ordo,poly,poly1,ctx);
			}
		fmpz_mul(a,a,l);
	}

}

int main(){
	fmpz_t p;fmpz_init(p);fmpz_set_ui(p,101);
	signed long d = 1;char *var="x";fq_ctx_t ctx;
	fq_ctx_init(ctx, p, d, var);
	
	printf("hello\n");
	point P;fq_init_point(P,ctx);elliptic_curve E;fq_init_curve(E,ctx);
	fq_set_ui(P.x, 1, ctx);fq_set_ui(P.y, 0, ctx);

	printf("hello\n");
	fq_ctx_clear(ctx);
	fmpz_clear(p);
	return 0;
}
