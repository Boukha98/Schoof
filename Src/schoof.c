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

