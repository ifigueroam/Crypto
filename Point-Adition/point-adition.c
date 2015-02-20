/*
 * =====================================================================================
 *
 *       Filename:  point-adition.c
 *
 *    Description:  Algoritmo Point Adition in Jacobian - Affine Coordinates
 *		    Para curvas elipticas de tipo y^2=x^3-3x+b	
 *                  
 * 		    Fórmula P+Q punto de adicion simple
 * Donde:
 *
 *  alpha = Z1^3*Y2 - Z2^3*Y1 
 *  beta  = Z1^2*X2 - Z2^2*X1
 *  X3    = alpha^2 - beta^3 - 2*X2^2*X2
 *  Y3    = alpha*(X2^2*X1*beta^2 -X3) - z2^3*Y1*beta^3
 *  Z3    = z1*z2*beta
 *
 * ser transformadas como:
 *    X3 = X3/(Z3^2) = X3*(Z3^2)^{-1}
 *    Y3 = X3/(Z3^3) = Y3*(Z3^3)^{-1}
 *    Z3 = Z3/Z3 = 1
 *
 * Obs: las formulas de MAGMA son entregadas en coordenadas
 * Affine, luego la salida de las formulas explicitas anteriores deben
 * 		
 *        Version:  1.0
 *        Created:  24/10/13 15:22:21
 *       Revision:  none
 *       Compiler:  gcc
 *       	    
 *       	    El programa se compila como: gcc point-adition.c -o padition -lgmp
 *       	    El programa se ejecuta como: ./padition
 *
 *         Author:  Iván Figueroa,
 *   Organization:  LCC
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include "time.h"
#define MILLION 1000000;
time_586 start, end;

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>


int main(void){
    //declara variables nomales
    double sec; 
		
    //declara variables mpz_t
    mpz_t P; 	// p:= numero primo para la definición del cuerpo primo
    mpz_t X1;	// X1, Y1, Z1 coordenadas del punto P de la curva Eliptica
    mpz_t Y1;
    mpz_t Z1;   // valor de Z1 = 1
    mpz_t X2;
    mpz_t Y2;
    mpz_t Z2;   // valor de Z2 = 1
    mpz_t AUX;  // variable auxuliar
    mpz_t ALFA;
    mpz_t BETA;
    mpz_t POW;

    mpz_t X3;
    mpz_t Y3;
    mpz_t Z3;

    //Inicio de variables.

    mpz_init(P);
    mpz_init(X1);
    mpz_init(Y1);
    mpz_init(Z1);
    mpz_init(X2);
    mpz_init(Y2);
    mpz_init(Z2);
    mpz_init(X3);
    mpz_init(Y3);
    mpz_init(Z3);

    mpz_init(ALFA);
    mpz_init(BETA);
    mpz_init(AUX);
    mpz_init(POW);

    //ALFA para multiplicar por 3
    mpz_set_str(ALFA,"3",10); 

    //numero primo para 192-bits de nivel de seguridad NIST
    mpz_set_str(P,  "6277101735386680763835789423207666416083908700390324961279", 10);
    
    //Coordenadas de un punto P en la curva EC.
    mpz_set_str(X1, "3055563715971644849423804859220265481316585815874058627244", 10); 

    mpz_set_str(Y1, "5301271154980119921208363180474818072856515133833450622956", 10);
    mpz_set_str(Z1, "1", 10);

    //Coordenadas de un punto P en la curva EC.
    mpz_set_str(X2, "4938825309294466132101294379631095895831893683098995938660", 10); 
    
    mpz_set_str(Y2, "2643698054889397759420891149012532976831983742899227222633", 10);
    mpz_set_str(Z2, "1", 10);

    mpz_set_str(X3 , "0" , 10);
    mpz_set_str(Y3 , "1" , 10);
    mpz_set_str(Z3 , "1" , 10);


    gmp_printf("\n P    = %Zd\n", P);
    gmp_printf("\n X1   = %Zd", X1);
    gmp_printf("\n Y1   = %Zd", Y1);
    gmp_printf("\n Z1   = %Zd", Z1);
    gmp_printf("\n X2   = %Zd", X2);
    gmp_printf("\n Y2   = %Zd", Y2);
    gmp_printf("\n Z2   = %Zd\n", Z2);
    gmp_printf("\n POW  = %Zd\n", POW);

    time(start);

    /* CÁLCULOS DE LAS FORMULAS EXPLICITAS */

    /* Cálculo de alpha:= Z1^3*Y2-Y1 */

    mpz_sub (ALFA,Y2,Y1); // Z1 vale 1, optimizo el cálculo
    mpz_mod(ALFA , ALFA , P);
    gmp_printf("\n alfa = %Zd\n", ALFA);

    /* Cálculo de beta:= Z1^2*X2 - Z^2*X1 */

    mpz_sub (BETA , X2 , X1); // Z1 y Z2 valen 1, optimizo el cálculo
    mpz_mod(BETA , BETA , P);
    gmp_printf(" beta = %Zd\n", BETA);
    
    /* Cálculo X3 = alpha^2 - beta^3 - 2*Z2^2*X1*beta^2 */

    mpz_set_str(POW, "2" , 10);
    mpz_powm(AUX,ALFA,POW,P);
    mpz_add(X3, X3, AUX);
    mpz_set_str(POW, "3" , 10);
    mpz_powm(AUX,BETA,POW,P);
    mpz_sub(X3, X3, AUX);
    mpz_mod(X3, X3, P);
    mpz_set_str(AUX, "2", 10);
    mpz_mul(AUX, AUX, BETA);
    mpz_mod(AUX, AUX, P);
    mpz_mul(AUX, AUX, BETA);
    mpz_mod(AUX, AUX, P);
    mpz_mul(AUX, AUX, X1);
    mpz_mod(AUX, AUX, P);
    mpz_sub(X3, X3, AUX);
    mpz_mod(X3, X3, P);
    gmp_printf("\n X3   = %Zd\n", X3);

    /* Cálculo Y3 = alpha*(Z2^2*X1*beta^2 - X3) - Z2^3*Y1*beta^3 */

    mpz_mul(AUX, BETA, BETA);
    mpz_mod(AUX, AUX, P);
    mpz_mul(Y3, AUX, X1);
    mpz_mod(Y3, Y3, P);
    mpz_sub(Y3, Y3, X3);
    mpz_mod(Y3, Y3, P);
    mpz_mul(Y3, Y3, ALFA);
    mpz_mul(AUX, AUX, BETA);
    mpz_mod(AUX, AUX, P);
    mpz_mul(AUX, AUX, Y1);
    mpz_mod(AUX, AUX, P);
    mpz_sub(Y3, Y3, AUX);
    mpz_mod(Y3, Y3, P);
    gmp_printf(" Y3   = %Zd\n", Y3);
    
    /* Calculo Z3 = Z1*beta */
    
    mpz_mul (Z3,Z1,BETA); 
    mpz_mod(Z3,Z3,P);
    gmp_printf(" Z3   = %Zd\n", Z3);

    printf("\nTRANSFORMACION COORDENADAS JACOBIAN TO AFFINE\n");
    
    mpz_mul(AUX, Z3, Z3);
    mpz_mod(AUX, AUX, P);
    mpz_invert(AUX, AUX, P);
    mpz_mul(X3, X3, AUX);
    mpz_mod(X3, X3, P);

    mpz_mul(AUX, Z3, Z3);
    mpz_mod(AUX, AUX, P);
    mpz_mul(AUX, AUX, Z3);
    mpz_mod(AUX, AUX, P);
    mpz_invert(AUX, AUX, P);
    mpz_mul(Y3, Y3, AUX);
    mpz_mod(Y3, Y3, P);

    mpz_div(Z3, Z3, Z3);
    mpz_mod(Z3, Z3, P);


    gmp_printf("\n X3   = %Zd", X3);
    gmp_printf("\n Y3   = %Zd", Y3);
    gmp_printf("\n Z3   = %Zd\n", Z3);

    time(end);
    sec = time_diff(end,start)/MILLION;
    printf("\n Tiempo de Ejecución: %.9lf\n\n",sec);

    /* free memory */

    mpz_clear(P);
    mpz_clear(X1);
    mpz_clear(Y1);
    mpz_clear(Z1);
    mpz_clear(AUX);
    mpz_clear(ALFA);
    mpz_clear(BETA);
    mpz_clear(X2);
    mpz_clear(Y2);
    mpz_clear(X3);
    mpz_clear(Y3);
    mpz_clear(Z3);

    return EXIT_SUCCESS;
}





