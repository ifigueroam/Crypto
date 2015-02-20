/*
 * =====================================================================================
 *
 *       Filename:  point-adition.c
 *
 *    Description:  Algoritmo Point Adition in Jacobian - Modified Coordinates
 *		    Para curvas elipticas de tipo y^2=x^3-3x+b	
 *                  
 * 		    Fórmula P+Q punto de adicion modificada
 * Donde:
 *
 *  U1    = X1*Z2^2 
 *  U2    = X2*Z1^2
 *  S1    = Y1*Z2^3
 *  S2    = Y2*Z1^3
 *  H     = U2 - U1
 *  R     = S2 - S1
 *  X3    = -H^3 - 2U1*H^2 + R
 *  Y3    = -S1*H^3 + R*(U1*H^2 - X3)
 *  Z3    = z1*z2*H
 *  Z4    = a*Z3^4
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
 *       	    El programa se compila como: gcc point-adition.c time.h -o padition -lgmp
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
    mpz_t U1;
    mpz_t U2;
    mpz_t S1;
    mpz_t S2;
    mpz_t H;
    mpz_t R;
    mpz_t POW; // variable potencia
    mpz_t A;   // variable representa a = -3, tambien como auxiliar

    mpz_t X3;
    mpz_t Y3;
    mpz_t Z3;
    mpz_t W;   // representa el valor de a*Z3^4

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

    mpz_init(U1);
    mpz_init(U2);
    mpz_init(S1);
    mpz_init(S2);
    mpz_init(R);
    mpz_init(H);
    mpz_init(AUX);
    mpz_init(POW);
    mpz_init(A);
    mpz_init(W);



    //numero primo para 512-bits de nivel de seguridad NIST
    mpz_set_str(P,  "6277101735386680763835789423207666416083908700390324961279", 10);
    
    //Coordenadas de un punto P en la curva EC.
    mpz_set_str(X1, "3055563715971644849423804859220265481316585815874058627244", 10); 
    mpz_set_str(Y1, "5301271154980119921208363180474818072856515133833450622956", 10);
    mpz_set_str(Z1, "1", 10);

    //Coordenadas de un punto Q en la curva EC.
    mpz_set_str(X2, "2603643047814032369115788269704066881802595496052627674363", 10); 
    mpz_set_str(Y2, "692745300022675652233671640120168188610163334847153663825", 10);
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

    time(start);

    /* CÁLCULOS DE LAS FORMULAS EXPLICITAS */

    /* Cálculo de U1:= X1*Z2^2 */

    mpz_set_str(POW, "2", 10);
    mpz_powm(AUX, Z2, POW, P);
    mpz_mul(U1, X1, AUX);
    mpz_mod(U1 ,U1 ,P);
    gmp_printf("\n U1 = %Zd", U1);

    /* Cálculo de U2:= X2*Z1^2 */

    mpz_powm(AUX, Z1, POW, P);
    mpz_mul(U2, X2, AUX);
    mpz_mod(U2 , U2 , P);
    gmp_printf("\n U2 = %Zd", U2);

    /* Cálculo de S1:= Y1*Z2^3 */

    mpz_set_str(POW, "3", 10);
    mpz_powm(AUX, Z2, POW, P);
    mpz_mul(S1, Y1, AUX);
    mpz_mod(S1 , S1 , P);
    gmp_printf("\n S1 = %Zd", S1);

    /* Cálculo de S2:= Y2*Z1^3 */

    mpz_powm(AUX, Z1, POW, P);
    mpz_mul(S2, Y2, AUX);
    mpz_mod(S2 , S2 , P);
    gmp_printf("\n S2 = %Zd", S2);

    
    /* Cálculo de H = U2 - U1 */

    mpz_sub(H, U2, U1);
    mpz_mod(H, H, P);
    gmp_printf("\n H  = %Zd", H);

    /* Cálculo de R = S2 - S1 */
    mpz_sub(R, S2, S1);
    mpz_mod(R, R, P);
    gmp_printf("\n R  = %Zd", R);

    /* Cálculo de X3 = -H^3 - 2*U1*H^2 + R^2 */

    mpz_neg(AUX, H);
    mpz_powm(AUX, AUX, POW, P);
    mpz_add(X3, X3, AUX);
    mpz_set_str(POW, "2", 10);
    mpz_powm(AUX, H, POW, P);
    mpz_set_str(A, "2", 10);
    mpz_mul(A,U1,A);
    mpz_mul(AUX, AUX, A);
    mpz_mod(AUX, AUX, P);
    mpz_sub(X3, X3, AUX);
    mpz_mod(X3, X3, P);
    mpz_powm(AUX,R,POW,P);
    mpz_add(X3, X3, AUX);
    mpz_mod(X3, X3, P);  
    gmp_printf("\n X3   = %Zd", X3);

    /* Cálculo de Y3 = -S1*H^3 + R*(U1*H^2 - X3) */
    
    mpz_powm(AUX, H, POW, P);
    mpz_mul(AUX, U1, AUX);
    mpz_sub(AUX, AUX, X3);
    mpz_mul(AUX, R, AUX);
    mpz_mod(AUX, AUX, P);
    mpz_set_str(POW, "3", 10);
    mpz_powm(Y3, H, POW, P);
    mpz_neg(S1, S1);
    mpz_mul(Y3, S1, Y3);
    mpz_add(Y3, AUX, Y3);
    mpz_mod(Y3, Y3, P);
    gmp_printf("\n Y3   = %Zd", Y3);
    
    
    /* Calculo Z3 = Z1*Z2*H */
    
    mpz_mul(AUX,Z3,H); // Z1 = 1, Z2 = 1 
    mpz_mod(Z3,AUX,P);
    gmp_printf("\n Z3   = %Zd", Z3);

    /* Cálculo de W = a*Z3^4 */

    mpz_set_str(POW, "4", 10);
    mpz_powm(AUX, Z3, POW, P);
    mpz_set_str(A, "-3", 10); 
    mpz_mul(AUX, A, AUX);
    mpz_mod(AUX, AUX, P);
    mpz_mod(W,AUX, P);
    gmp_printf("\n W    = %Zd\n", W);

    printf("\nTRANSFORMACION COORDENADAS JACOBIAN MODIFIED\n");
    
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
    mpz_clear(U1);
    mpz_clear(U2);
    mpz_clear(S1);
    mpz_clear(S2);
    mpz_clear(H);
    mpz_clear(R);
    mpz_clear(POW);
    mpz_clear(A);
    mpz_clear(X2);
    mpz_clear(Y2);
    mpz_clear(X3);
    mpz_clear(Y3);
    mpz_clear(Z3);

    return EXIT_SUCCESS;
}





