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

void fnext_prime(fmpz_t res, fmpz_t l){
	fmpz_add_ui(res, l, 2); // res = l+2
	while(!fmpz_is_prime(res)){ // 0 ou -1
		fmpz_add_ui(res, res, 2); // res = res+2
	}
}

void produit_l(fmpz_t a, fmpz_t lmax){
	fmpz_set_si(a, 6); fmpz_t(l);fmpz_init(l);fmpz_set_si(l, 3);
	while(!fmpz_cmp(a, q_sqrt)){
		fnext_prime(l, l);
		if(fmpz_equal(l, p)){
			fnext_prime(l, l);
		}
		fmpz_mul(a, a, l);
	}
	fmpz_set(lmax, l);
}	
	

void schoof(elliptic_curve E, fq_ctx_t ctx, fmpz_t q, fmpz_t c){
	fmpz_t l, q_sqrt;
	fq_t tmp, tmp1, t;
	fq_poly_t poly; fq_poly_init(poly, ctx);fq_poly_t poly1; 
	fq_poly_init(poly1, ctx);
	fq_init(l, ctx);fmpz_init(q_sqrt);
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
	
	//si pgcd=1 alors t=1, sinon t=0
	if(fq_poly_is_one(poly)){
		fq_set_ui(t, 1, ctx);
	}
	else{fq_set_ui(t, 0, ctx);}
	
	point P;
	fq_poly_t *Psi; //poly de division
    	div_poly(Psi, E, lmax+2, ctx, P);//remplissage du tableau 
    
   	//q_bar = q mod l
   	fmpz_t q_bar; fmpz_init(p_bar);
   	fmpz_mod(q_bar, q, l);
   	fmpz_div_si(tmp, l, 2); //tmp = l/2
   	if(fmpz_cmp(q_bar, tmp)>=0){
   		fmpz_sub(q_bar, q_bar, l);//q_bar = q_bar - l
   	}
   
   	/* Phi^2(P) = +-[q_bar]P ssi:
   	
   	 Si k pair: 
   	 	(x^q^2 - x)*(Psi_k)^2(x)*(x^3 + ax + b) + Psi_(k-1)(x)*Psi_(k+1)(x) = 0
   	 Si k impair: 
   	 	(x^q^2 - x)*(Psi_k)^2(x) + Psi_(k-1)(x)*Psi_(k+1)(x)*(x^3 + ax + b) = 0     
    	*/
  
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
    /* 1) Si pgcd(poly, Psi[l]) \neq 1 alors:
     il existe un point P non trivial de E[l] avec (Phi_l)^2*P = ±q_bar*P
    	-Si (Phi_l)^2*P = -q_bar*P Alors t = 0[l]
     	-Si (Phi_l)^2*P = +q_bar*P Alors t^2 = 4*q_bar[l] 
     		--> q_bar RC ou non?		(Jacobi)
   		2) Sinon (pgcd=1): alors t \neq 0[l], on peut appliquer les formules
        	d'addition de deux points distincts[ (Phi_l)^2(x,y) + q(x,y) ]
    */

    fq_poly_gcd(poly, poly, Psi[l], ctx);
    //Cas ou pgcd \neq 1
    if(!fq_poly_is_one(poly, ctx)){
	if(fmpz_jacobi(q, l) == -1){//q n'est pas un résidu quadratique
            fmpz_set_si(tmp, 0); // tmp = 0
          fmpz_CRT(t, t, a, tmp, l, 1); // TRC: t=0 mais on fait le calcul
        }
        //On cherche w tel que q = w²[l]
      	else{
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
            //Cas pgcd=1: t=0
            if(fq_poly_is_one(poly, ctx)){
           	fq_set_ui(t, 0, ctx);
            }
            //Cas pgcd \neq 1
            else{
            	//Cas w impair
            	if(w&1){
           		fmpz_add_ui(tmp, q, 3);
    			fmpz_divexact_si(tmp, tmp, 2);
    
    			fq_poly_gcd(poly, poly, Psi[l], ctx);
    	//Cas ou pgcd \neq 1
    	if(!fq_poly_is_one(poly, ctx)){
		    if(fmpz_jacobi(q, l) == -1){//q n'est pas un résidu quadratique
                fmpz_set_si(tmp, 0); // tmp = 0
                fmpz_CRT(t, t, a, tmp, l, 1); // TRC: t=0 mais on fait le calcul
            }
            //On cherche w tel que q = w²[l]
            else{
                fmpz_sqrtmod(w, q, l);
                //Cas w impair
                if(w&1){
                    fq_poly_sqr(poly, Psi[w], ctx);
    				fq_poly_mul(poly, poly, Phi1, ctx);
    				fq_poly_mul(poly1, Psi[w-1], Psi[w+1], ctx);
    				fq_poly_mul(poly1, poly1, y2, ctx);
    				fq_poly_add(poly, poly, poly1, ctx);
>>>>>>> 76473c80101d3089cc4ceb2176ba6ee2c0a72dab
    			}
    			//Cas w pair
    			else{
    				fmpz_sub_ui(tmp, q, 1);
    				fmpz_divexact_si(tmp, tmp, 2);
    			}
    			tmp = fmpz_get_ui(tmp);
    			fq_poly_pow(poly, y2, tmp, ctx);
    			fq_poly_pow(poly1, Psi[w], 3, ctx);
    			fq_poly_mul(poly, poly, poly1, ctx);
    			fq_set_ui(tmp1, 4, ctx);
    			fq_poly_scalar_mul_fq(poly, poly, 4, ctx);
    			
    			fq_poly_pow(poly1, Psi[w+2], 2, ctx);
    			fq_poly_mul(poly1, poly1, Psi[w-1], ctx);
    			
    			fq_poly_sub(poly, poly, poly1, ctx);
    			
    			fq_poly_pow(poly1, Psi[w-2], 2, ctx);
    			fq_poly_mul(poly1, poly1, Psi[w+1], ctx);
    			
    			fq_poly_add(poly, poly, poly1, ctx);
    			
    			//Calcul du pgcd
<<<<<<< HEAD
       			fq_poly_gcd(poly, poly, Psi[l], ctx);
            		
    			fmpz_mul_ui(tmp, w, 2);//tmp=2*w
    			//pgcd = 1: t = -2w[l]
    			//pgcd \neq 1: t = 2w[l]  
    			if(fq_poly_is_one(poly, ctx)){
            		fmpz_negmod(tmp, tmp, l);
=======
                fq_poly_gcd(poly, poly, Psi[l], ctx);
                //Cas pgcd=1: t=0
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
    				tmp = fmpz_get_ui(tmp);
    				fq_poly_pow(poly, y2, tmp, ctx);
    				fq_poly_pow(poly1, Psi[w], 3, ctx);
    				fq_poly_mul(poly, poly, poly1, ctx);
    				fq_set_ui(tmp1, 4, ctx);
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
    				//pgcd = 1: t = -2w[l]
    				//pgcd \neq 1: t = 2w[l]  
    				if(fq_poly_is_one(poly, ctx)){
                        fmpz_negmod(tmp, tmp, l);
            		}
      				fmpz_CRT(t, t, a, tmp, l, 1);//TRC
>>>>>>> 76473c80101d3089cc4ceb2176ba6ee2c0a72dab
            	}
      			fmpz_CRT(t, t, a, tmp, l, 1);//TRC
            }
        }
    }
    //Cas ou Phi^2(P) \neq +-[q_bar]P
    else{
    	//Construction de X^q² + X^q + X = Phi3
<<<<<<< HEAD
    	fq_poly_t Phi3;
=======
        fq_poly_t Phi3;
>>>>>>> 76473c80101d3089cc4ceb2176ba6ee2c0a72dab
		fq_poly_init(Phi3, ctx);
		fq_set_ui(tmp, 1, ctx);
		fq_poly_set_coeff(Phi3, fmpz_get_si(q)*fmpz_get_si(q), tmp, ctx);
		fq_poly_set_coeff(Phi3, fmpz_get_si(q), tmp, ctx);
		fq_poly_set_coeff(Phi3, 1, tmp, ctx);
		
		fq_poly_t beta2, alpha;
		fq_poly_init(beta2, ctx);fq_poly_init(alpha, ctx);
		
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
		fq_poly_scalar_mul_fq(beta2, poly, tmp1, ctx);
			
		//Calcul de alpha
		fq_poly_sqr(poly, Psi[q_bar-1], ctx);
		fq_mul_mul(poly, poly, Psi[q_bar+2], ctx);
		
		fq_poly_sqr(poly1, Psi[q_bar+1], ctx);
		fq_mul_mul(poly1, poly, Psi[q_bar-1], ctx);
		
		fq_poly_sub(poly, poly, poly1, ctx);
		
		fq_set_si(tmp1, 3, ctx);
		fq_poly_pow(poly1, Psi[q_bar], tmp1, ctx);
		
		fmpz_pow_si(tmp, q_bar, 2);
		fmpz_add_si(tmp, tmp, 1);
		fmpz_divexact_si(tmp, tmp, 2);
		
		fq_poly_pow(alpha, y2, fmpz_get_si(tmp), ctx);
		fq_poly_mul(alpha, alpha, poly1, ctx);
		fq_set_si(tmp1, 4, ctx);
		fq_poly_scalar_mul(alpha, alpha, 4, ctx);
			
		fq_poly_sub(alpha, poly, alpha, ctx);
			
		for(j = 1; j <= l-1; j++){
			fq_poly_mul(poly, Psi
			
		}
			
    	
    	
    	
	
}
