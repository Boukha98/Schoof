#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>

//int  fmpz_is_prime(const  fmpz_t l) où l>2

/* Si l est premier, retourne 1. 
 * Si l est composé, retourne 0. 
 * Sinon retourne -1.
*/

void fnext_prime(fmpz_t res, fmpz_t l){
	fmpz_add_ui(res, l, 2); // res = l+2
	while(!fmpz_is_prime(res)){ // 0 ou -1
		fmpz_add_ui(res, res, 2); // res = res+2
	}
}

