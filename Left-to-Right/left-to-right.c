/*
 * =====================================================================================
 *
 *       Filename:  left-to-right.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/13 16:28:49
 *       Revision:  none
 *       Compiler:  gcc
 * 		Compilation 
 * 	 & Ejecuction:	gcc left-to-right.c -o ltr -lgmp -std=c99 && ./ltr
 *
 *         Author:  Iván Figueroa, 
 *   Organization:  LCC
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "time.h"

#define TIME 1000000;
time_686 start, end;

void bin_to_vector(int *vector_bin, int slots, mpz_t bits, mpz_t binary){

	mpz_t aux,div,abin,cont,cont2,mod;

	mpz_init_set_str(aux,   "1" , 10);
	mpz_init_set_str(div,   "10", 10);
	mpz_init_set_str(cont,  "0" , 10);
	mpz_init_set_str(cont2, "0" , 10);
	mpz_init_set_str(mod,   "1" , 10);
	mpz_init_set_str(abin,  "1" , 10);

	mpz_set (abin, binary);

	/* generacion del modulo respecto al tamaño del binario */
	while(mpz_cmp(bits, cont2) > 0){
		
		mpz_mul(mod,div,mod);
		mpz_add_ui(cont2,cont2,1);
	}
	mpz_div(mod,mod,div);

	int i = 0;

	/* Transformacion: almacenamos el binario en un arreglo */
	do{

		if(mpz_cmp(abin,mod) >= 0){ 
			mpz_div(aux,abin,mod);
			mpz_sub(abin,abin,mod);
			vector_bin[i] = mpz_get_ui(aux);
			mpz_div(mod,mod,div);
		
			mpz_add_ui(cont,cont,1);
			i++;
		}
		else{ 
			mpz_div(aux,abin,mod);
			vector_bin[i] = mpz_get_ui(aux);
			mpz_div(mod,mod,div);

			mpz_add_ui(cont,cont,1);
			i++;
		}

	}while(mpz_cmp(bits, cont) >0);

	/* liberar memoria */

	mpz_clear (aux);
	mpz_clear (div);
	mpz_clear (cont);
	mpz_clear (cont2);
	mpz_clear (mod);
	mpz_clear (abin);
}

void num_bits(mpz_t bits, mpz_t number){

	mpz_t pow, ground, counter, module;

	mpz_init_set_str (counter,"100", 													10); /* numero a comparar */
	mpz_init_set_str (pow, 		"0", 														10); /* exponente */
	mpz_init_set_str (ground, "2", 														10); /* base */
	mpz_init_set_str (module, "1000000000000000000000000000", 10); /* modulo */

	/* Calculamos el tamaño del arreglo en base a el valor binario obtenido */
	do{
			mpz_powm(counter,ground,pow,module);
			mpz_add (pow,bits,pow);
	}while(mpz_cmp(number, counter) >= 0);
	mpz_sub(bits,pow,bits);

	/* liberar memoria de la funcion */
	mpz_clear (pow);
	mpz_clear (ground);
	mpz_clear (counter);
	mpz_clear (module);
}

void radix(mpz_t binary, mpz_t number, mpz_t foot, mpz_t bits){
    
	int zero;
  mpz_t a,x,bin,q,aux2,cero,base;

  mpz_init (a);
  mpz_init (x);
  mpz_init (bin);
  mpz_init (q);
  mpz_init (aux2);
  mpz_init (cero);
  mpz_init (bin);
	mpz_init (base);

  zero = 0;

	mpz_set(a   ,number);
	mpz_set(base,foot);
  mpz_set(x 	, a);

	mpz_set_str(cero , "1"  , 10);
  mpz_set_str(bin  , "0"  , 10);
  mpz_set_str(a  	 , "0"  , 10);

	/* Calculo de RADIX  */
  mpz_div (q, x, base); 
  mpz_mul (aux2, q, base); 
  mpz_sub (a, x, aux2);
  if(mpz_cmp_ui(a, zero)==0)
  	mpz_mul_ui (cero, cero, 10);

  while(mpz_cmp_ui(q, zero)>0){
  	if(mpz_cmp_ui(q, zero)==0)
    	break;
   	mpz_div (q, x, base);
    mpz_set (x, q);
    mpz_div (q, x, base); 
    mpz_mul (aux2, q, base); 
    mpz_sub (aux2, x, aux2);
    mpz_mul_ui (a, a, 10);
    mpz_add (a, a, aux2);

    if(mpz_cmp_ui(a, zero)==0)
    	mpz_mul_ui (cero, cero,10);

  }
	while(mpz_cmp_ui(a, zero)!=0){
		mpz_mod_ui (aux2, a, 10);
    mpz_div_ui (a, a, 10);
    mpz_mul_ui (bin, bin,10);
    mpz_add (bin, bin, aux2);
  }
    
	/* Retorno de Conversion*/
  mpz_mul(binary, bin, cero);

	/* Liberar Memoria */
  mpz_clear (a);
  mpz_clear (x);
  mpz_clear (q);
  mpz_clear (aux2);
  mpz_clear (cero);
  mpz_clear (bin);
	mpz_clear (base);
}

void add_jacobian(mpz_t X2, mpz_t Y2, mpz_t Z2, mpz_t X1, mpz_t Y1, mpz_t Z1, mpz_t P){

	mpz_t ALFA,BETA,AUX,X3,Y3,Z3;

	mpz_init_set_str(ALFA, "3", 10);
	mpz_init_set_str(BETA, "1", 10);
	mpz_init_set_str(AUX,  "3", 10);
	mpz_init (X3);
	mpz_init (Y3);
	mpz_init (Z3);

	/* Cálculo de alpha:= Z1^3*Y2-Y1 */

  mpz_sub (ALFA,Y2,Y1); // Z1 vale 1, optimizo el cálculo
  mpz_mod(ALFA , ALFA , P);

  /* Cálculo de beta:= Z1^2*X2 - Z^2*X1 */

  mpz_sub (BETA , X2 , X1); // Z1 y Z2 valen 1, optimizo el cálculo
  mpz_mod(BETA , BETA , P);
    
  /* Cálculo X3 = alpha^2 - beta^3 - 2*Z2^2*X1*beta^2 */

  mpz_mul(X3, ALFA, ALFA);
  mpz_mod(X3, X3, P);
	mpz_powm(AUX,BETA,AUX,P);
  mpz_sub(X3, X3, AUX);
  mpz_mod(X3, X3, P);
  
	mpz_set_str(AUX, "2", 10);
	mpz_powm(AUX,BETA,AUX,P);
  mpz_mul(AUX, AUX, X1);
	mpz_mul_ui(AUX,AUX,2);
	mpz_mod(AUX, AUX, P);
  mpz_sub(X3, X3, AUX);
  mpz_mod(X3, X3, P);

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
    
 	/* Calculo Z3 = Z1*beta */
    
  mpz_mul(Z3,Z1,BETA); 
  mpz_mod(Z3,Z3,P);

	/* Traspaso de variables */
	
	mpz_add_ui(X2,X3,0);
	mpz_add_ui(Y2,Y3,0);
	mpz_add_ui(Z2,Z3,0);
	
	/* gmp_printf("\nTRANSFORMACIÓN COORDENADAS JACOBIAN TO AFFINE */
     
    mpz_t(Z32);
    mpz_t(Z33);
    mpz_init(Z32);
    mpz_init(Z33);
      
    mpz_mul(Z32, Z2, Z2); //Z3^2
    mpz_mod(Z32, Z32, P);
    mpz_mul(Z33, Z2, Z32); //Z3^3
    mpz_mod(Z33, Z33, P); //reducciendo Z2^3 mod p
    
    mpz_invert(Z32, Z32, P); //calcula en inverso multiplicativo en mod p
    mpz_mul(X2, X2, Z32); //X3/Z3^2  Coordenada Affine 
    mpz_mod(X2, X2, P);   //reducción de la coordenada mod p
   
    mpz_invert(Z33, Z33, P); //inverso multiplicativo de Z3^3 mod p 
    mpz_mul(Y2, Y2, Z33);
    mpz_mod(Y2, Y2, P);   //reducción de la coordenada mod p
    
    mpz_cdiv_q(Z2, Z2, Z2);

	/* Liberar Memoria transformacion Afine */

	mpz_clear(Z32);
	mpz_clear(Z33);
	
	/* Liberar Memoria */
	mpz_clear (AUX);
	mpz_clear (ALFA);
	mpz_clear (BETA);
	mpz_clear (X3);
	mpz_clear (Y3);
	mpz_clear (Z3);
}


void doub_jacobian(mpz_t X2, mpz_t Y2, mpz_t Z2, mpz_t P){

 	mpz_t ALFA, AUX, AUX2, BETA, X3,Y3,Z3;

	mpz_init_set_str(ALFA, "3", 10);
	mpz_init_set_str(AUX,  "1", 10);	
	mpz_init_set_str(AUX2, "1", 10);	
	mpz_init_set_str(BETA, "4", 10);
	mpz_init (X3);
	mpz_init (Y3);
	mpz_init (Z3);
    
	/*Calculo de alpha:= 3*(X1+Z1^2)*(X1-Z1^2)*/
  
  mpz_add(ALFA,X2,Z2); //X1+Z1^2
 	mpz_mul_ui(ALFA, ALFA, 3); //3*(X1+Z1^2)
  mpz_mod(ALFA,ALFA,P);  
  mpz_sub(AUX,X2,Z2);    //X1-Z1^2
  mpz_mul(ALFA,ALFA,AUX);//3*(X1+Z1^2)*(X1-Z1^2)
  mpz_mod(ALFA,ALFA,P);

  /* Calculo de beta:=4*X1*Y1^2 */
      
  mpz_mul(AUX, Y2, Y2);  /*Y1^2*/
  mpz_mul(AUX, X2, AUX); /*X1Y1^2*/  
  mpz_mul_ui(BETA, AUX, 4); /*2X1Y1^2*/  
	mpz_mod(BETA, BETA, P);     
	
	/*Calculo X3:=alpha^2-2beta;*/
    
  mpz_mul(AUX, ALFA, ALFA);  /*alpha^2*/ 
  mpz_add(AUX2, BETA, BETA); /*2beta*/
  mpz_sub(X3, AUX, AUX2);  /*X3*/
  mpz_mod(X3, X3, P);    

  /*  Calculo de Y3:=alpha*(beta-X3)-8Y1^4;*/
     
  mpz_set_str(AUX,  "4", 10);
  mpz_powm(AUX2,Y2,AUX,P);
	mpz_mul_ui(AUX2,AUX2,8);
	mpz_sub(AUX,BETA,X3);
	mpz_mul(AUX,AUX,ALFA);
	mpz_sub(Y3,AUX,AUX2);
	mpz_mod(Y3,Y3,P);
 
  /*Calculo de Z3:=2*Y1*Z1*/
   
  mpz_mul(Z3, Y2, Z2);  /*Y1Z1*/
  mpz_add(Z3, Z3, Z3);  /*2Y1Z1*/
  mpz_mod(Z3, Z3, P);   /*Z3 mod p*/

	/* Traspaso de variables */
	mpz_add_ui(X2,X3,0);
	mpz_add_ui(Y2,Y3,0);
	mpz_add_ui(Z2,Z3,0);
	
	/* TRANSFORMACIÓN COORDENADAS JACOBIAN TO AFFINE */
     
    mpz_t(Z32);
    mpz_t(Z33);
    mpz_init(Z32);
    mpz_init(Z33);
      
    mpz_mul(Z32, Z2, Z2); //Z3^2
    mpz_mod(Z32, Z32, P);
    mpz_mul(Z33, Z2, Z32); //Z3^3
    mpz_mod(Z33, Z33, P); //reducciendo Z2^3 mod p
    
    mpz_invert(Z32, Z32, P); //calcula en inverso multiplicativo en mod p
    mpz_mul(X2, X2, Z32); //X3/Z3^2  Coordenada Affine 
    mpz_mod(X2, X2, P);   //reducción de la coordenada mod p
   
    mpz_invert(Z33, Z33, P); //inverso multiplicativo de Z3^3 mod p 
    mpz_mul(Y2, Y2, Z33);
    mpz_mod(Y2, Y2, P);   //reducción de la coordenada mod p
    
    mpz_cdiv_q(Z2, Z2, Z2);

	/* Liberar Memoria transformacion Afine */

	mpz_clear(Z32);
	mpz_clear(Z33);
	
	/* Liberar Memoria */
	mpz_clear (AUX);
	mpz_clear (AUX2);
	mpz_clear (ALFA);
	mpz_clear (BETA);
	mpz_clear (X3);
	mpz_clear (Y3);
	mpz_clear (Z3);
}


int main (){
  /* declaracion de variables del programa */
	int slots, A=0, D=0;	
	double sec;
  mpz_t number, binary, foot, bits, P, X1, Y1, Z1, X2, Y2, Z2;

  /* inicio de variables de gmp */
  
	mpz_init_set_str (number, "599999999998", 10);
  mpz_init_set_str (foot, "2", 10);
	mpz_init_set_str (bits, "1", 10); /* cantidad de bits de memoria a asignar */
  mpz_init (binary);
	

 	/* p, número generador del cuerpo  */
  mpz_init_set_str(P,  "6277101735386680763835789423207666416083908700390324961279", 10);

  /* Coordenadas de un punto P en EC. */
  mpz_init_set_str(X1, "3055563715971644849423804859220265481316585815874058627244", 10); 
  mpz_init_set_str(Y1, "5301271154980119921208363180474818072856515133833450622956", 10);
  mpz_init_set_str(Z1, "1", 10);

  /* Coordenadas de un punto Q en EC. definidas por Q <- oo */
  mpz_init_set_str(X2, "1", 10); 
  mpz_init_set_str(Y2, "1", 10);
  mpz_init_set_str(Z2, "0", 10);

	/* Imprimir Coordenadas */
	gmp_printf("\n generador del cuerpo: ");
  gmp_printf("\n P  = %Zd\n", P);
	
	gmp_printf("\n Coordenadas del punto P en la Curva Elíptica: \n");
  gmp_printf(" Xp = %Zd\n", X1);
  gmp_printf(" Yp = %Zd\n", Y1);
  gmp_printf(" Zp = %Zd\n", Z1);
  
	gmp_printf("\n Coordenadas del punto Q en la Curva Elíptica: \n");
	gmp_printf(" Xq = %Zd\n", X2);
  gmp_printf(" Yq = %Zd\n", Y2);
  gmp_printf(" Zq = %Zd\n", Z2);
	
	/* LLamada a RADIX */
  radix(binary, number, foot, bits);
	
	/* llamada a num_bits para calcular el tamaño del arreglo */
	num_bits(bits, number);
  
	gmp_printf("\n Number       = %Zd\n",number); 	/* numero decimal empleado */
	gmp_printf(" Radix Number = %Zd\n",binary); 	/* numero binario convertido */
	gmp_printf(" Size_Array   = %Zd\n",bits); 	/* tamaño en bits a ocupar en el arreglo */

	/* Llamada a bin_to_vector */

	slots = mpz_get_ui(bits);
	int vector_bin[slots];
	
	bin_to_vector(vector_bin, slots, bits, binary);

	/*  Ejecución de LEFT_TO_RIGTH */
	time(start);
	
	mpz_add_ui(X2,X1,0);
	mpz_add_ui(Y2,Y1,0);
	mpz_add_ui(Z2,Z1,0);
	

  for (int i=1; i<=slots-1; i++){
 		doub_jacobian(X2,Y2,Z2,P);
		D++;
	 	if(vector_bin[i] == 1){
  		add_jacobian(X2,Y2,Z2,X1,Y1,Z1,P); 
			A++;
		}
	}
	time(end);

	/* Resultados de ECSM */
	gmp_printf("\n Resultado SMEC kP\n");
  gmp_printf(" Xq = %Zd\n",X2); 
  gmp_printf(" Yq = %Zd\n",Y2); 
  gmp_printf(" Zq = %Zd\n",Z2); 

	printf("\n Total de Adiciones:  %d",A);
	printf("\n Total de Doblados:   %d",D);

	sec = time_diff(end,start)/TIME;
	printf("\n Tiempo de Ejecución: %.9lf\n\n",sec);

  /* Liberar Memoria  */ 
  mpz_clear (number);
  mpz_clear (foot);
  mpz_clear (binary);
	mpz_clear (bits);
	mpz_clear (P);
	mpz_clear (X1);
	mpz_clear (X2);
	mpz_clear (Y1);
	mpz_clear (Y2);
	mpz_clear (Z1);
	mpz_clear (Z2);

  return EXIT_SUCCESS;
}
