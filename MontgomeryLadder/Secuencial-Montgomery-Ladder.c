/*
 * =====================================================================================
 *
 *       Filename:  Secuencial-Montgomery-Ladder.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/13 16:28:49
 *       Revision:  none
 *       Compiler:  gcc
 * 		Compilation 
 * 	 & Ejecuction:	gcc Secuencial-Montgomery-Ladder.c -o sml -lgmp -std=c99 && ./sml
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

	//int slots;
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

		if(mpz_cmp(abin,mod) >= 0){ //abin/mod = 1
			mpz_div(aux,abin,mod);
			mpz_sub(abin,abin,mod);
			vector_bin[i] = mpz_get_ui(aux);
			mpz_div(mod,mod,div);
		
			mpz_add_ui(cont,cont,1);
			i++;
		}
		else{ //cambiar condicion en 0 no es igual
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

	mpz_init_set_str (counter,"100", 						10); /* numero a comparar */
	mpz_init_set_str (pow, 		"0", 							10); /* exponente */
	mpz_init_set_str (ground, "2", 							10); /* base */
	mpz_init_set_str (module, "10000000000000", 10); /* modulo */

	/* Calculamos el tamaño del arreglo en base a el valor binario obtenido */
	do{
			mpz_powm(counter,ground,pow,module);
			mpz_add (pow,bits,pow);
	}while(mpz_cmp(number, counter) > 0);
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
  mpz_div (q, x, base); //ui
  mpz_mul (aux2, q, base); //ui
  mpz_sub (a, x, aux2);
  if(mpz_cmp_ui(a, zero)==0)
  	mpz_mul_ui (cero, cero, 10);

  while(mpz_cmp_ui(q, zero)>0){
  	if(mpz_cmp_ui(q, zero)==0)
    	break;
   	mpz_div (q, x, base); //ui
    mpz_set (x, q);
    mpz_div (q, x, base); //ui
    mpz_mul (aux2, q, base); //ui
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
  // mpz_clear (bin);
  mpz_clear (q);
  mpz_clear (aux2);
  mpz_clear (cero);
  mpz_clear (bin);
	mpz_clear (base);
}

void montgomery_ladder(mpz_t Xp, mpz_t Zp, mpz_t Xq, mpz_t Zq, mpz_t P, mpz_t a, mpz_t b, mpz_t xd, int R){

	/* Declaración de variables de gmp */
	mpz_t R0,R1,R2,R3,R4,R5,R6,R7;

	/* Inicio y Asignación de Variables */
	
	mpz_init(R0);
	mpz_init(R1);
	mpz_init(R2);
	mpz_init(R3);
	mpz_init(R4);
	mpz_init(R5);
	mpz_init(R6);
	mpz_init(R7);

	mpz_set(R0,Xp);
	mpz_set(R1,Zp);
	mpz_set(R2,Xq);
	mpz_set(R3,Zq);

	/* Inicio del Algoritmo Montgomery Ladder */

	mpz_mul(R6,R2,R1);
	mpz_mul(R7,R3,R0);
	mpz_add(R4,R7,R6);
	mpz_sub(R5,R7,R6);
	mpz_mul(R5,R5,R5);
	mpz_mul(R7,R1,R3);
	mpz_mul(R1,a,R7);
	mpz_mul(R6,R7,R7);
	mpz_mul(R0,R0,R2);
	mpz_mul(R6,b,R6);
	mpz_add(R0,R0,R1);
	mpz_add(R6,R6,R6);
	mpz_mul(R0,R0,R4);
	mpz_mul(R1,xd,R5);
	mpz_add(R4,R0,R6);
	mpz_add(R4,R4,R4);
	mpz_add(R6,R2,R2);
	mpz_sub(R4,R4,R1);
	mpz_add(R7,R3,R3);
	mpz_mul(R0,R6,R7);
	mpz_mul(R1,R3,R3);
	mpz_mul(R2,R2,R2);
	mpz_mul(R3,a,R1);
	mpz_sub(R6,R2,R3);
	mpz_add(R7,R2,R3);
	mpz_add(R1,R1,R1);
	mpz_mul(R2,b,R1);
	mpz_mul(R7,R7,R0);
	mpz_mul(R1,R2,R1);
	mpz_mul(R0,R0,R2);
	mpz_mul(R6,R6,R6);
	mpz_sub(R6,R6,R0);
	mpz_add(R7,R7,R1);
	
	/* reducir valores al cuerpo */

	mpz_mod(R4,R4,P);
	mpz_mod(R5,R5,P);
	mpz_mod(R6,R6,P);
	mpz_mod(R7,R7,P);
	  
	gmp_printf("\n VALORES con Reduccion\n");
	gmp_printf(" Xp = %Zd\n",Xp);
	gmp_printf(" Zp = %Zd\n",Zp);
	gmp_printf(" Xq = %Zd\n",Xq);
	gmp_printf(" Zq = %Zd\n",Zp);
	
	/* Reasignar variables */
	if (R == 1){
		mpz_add_ui(Xp,R4,0);
		mpz_add_ui(Zp,R5,0);
		mpz_add_ui(Xq,R6,0);
		mpz_add_ui(Zq,R7,0);
	/*  
	gmp_printf("\n VALORES R=1\n");
	gmp_printf(" Xp = %Zd\n",Xp);
	gmp_printf(" Zp = %Zd\n",Zp);
	gmp_printf(" Xq = %Zd\n",Xq);
	gmp_printf(" Zq = %Zd\n",Zp);*/
	}

	if (R == 2){
		mpz_add_ui(Xp,R6,0);
		mpz_add_ui(Zp,R7,0);
		mpz_add_ui(Xq,R4,0);
		mpz_add_ui(Zq,R5,0);
	/*  
	gmp_printf("\n VALORES R=2\n");
	gmp_printf(" Xp = %Zd\n",Xp);
	gmp_printf(" Zp = %Zd\n",Zp);
	gmp_printf(" Xq = %Zd\n",Xq);
	gmp_printf(" Zq = %Zd\n",Zp);*/
	}

	/* Transformacion a Coordenadas Afines */
	mpz_invert(Zp,Zp,P);
	mpz_mul(Xp,Xp,Zp);
	mpz_mod(Xp,Xp,P);
	mpz_cdiv_q(Zp,Zp,Zp);

	mpz_invert(Zq,Zq,P);
	mpz_mul(Xq,Xq,Zq);
	mpz_mod(Xq,Xq,P);
	mpz_cdiv_q(Zq,Zq,Zq);

	gmp_printf("\n VALORES en Coordenadas Afines\n");
	gmp_printf(" Xp = %Zd\n",Xp);
	gmp_printf(" Zp = %Zd\n",Zp);
	gmp_printf(" Xq = %Zd\n",Xq);
	gmp_printf(" Zq = %Zd\n",Zp);
	
	/*  Liberar Memoria */
	
	mpz_clear(R0);
	mpz_clear(R1);
	mpz_clear(R2);
	mpz_clear(R3);
	mpz_clear(R4);
	mpz_clear(R5);
	mpz_clear(R6);
	mpz_clear(R7);
}
int main (){
  /* declaracion de variables del programa */
	int slots, R;	
	double sec;
  mpz_t number, binary, foot, bits, P, Xp, Yp, Zp, Xq, Yq, Zq, a, b, xd;

  /* inicio de variables de gmp */
  
	/* Tamaño maximo permitido: 599999999998 */
	mpz_init_set_str (number, "5", 10);
  mpz_init_set_str (foot, "2", 10);
	mpz_init_set_str (bits, "1", 10); /* cantidad de bits de memoria a asignar */
  mpz_init (binary);
	

 	/* p, número de 192 bits generador del cuerpo definido bajo el estandar NIST */
  mpz_init_set_str(P,  "6277101735386680763835789423207666416083908700390324961279", 10);

	/* a, b parámetros de la curva elíptica  */
  mpz_init_set_str(a, "6277101735386680763835789423207666416083908700390324961276", 10);
  mpz_init_set_str(b, "2455155546008943817740293915197451784769108058161191238065", 10);
  mpz_init_set_str(xd, "3055563715971644849423804859220265481316585815874058627244", 10);  

	/* Coordenadas de un punto P en EC. */
  mpz_init_set_str(Xp, "3055563715971644849423804859220265481316585815874058627244", 10); 
  mpz_init_set_str(Yp, "5301271154980119921208363180474818072856515133833450622956", 10);
  mpz_init_set_str(Zp, "1", 10);

  /* Coordenadas de un punto Q en EC. definidas por Q <- oo */
  mpz_init_set_str(Xq, "2308451064198626516353983124879239968795376585073748052063", 10); 
  mpz_init_set_str(Yq, "1366591923459201440060947799383796565631656810572609006762", 10);
  mpz_init_set_str(Zq, "1", 10);

	/* Imprimir Coordenadas */
	gmp_printf("\n generador del cuerpo: ");
  gmp_printf("\n P  = %Zd\n", P);
	
	gmp_printf("\n Coordenadas del punto P en la Curva Elíptica: \n");
  gmp_printf(" Xp = %Zd\n", Xp);
  gmp_printf(" Yp = %Zd\n", Yp);
  gmp_printf(" Zp = %Zd\n", Zp);
  
	gmp_printf("\n Coordenadas del punto Q en la Curva Elíptica: \n");
	gmp_printf(" Xq = %Zd\n", Xq);
  gmp_printf(" Yq = %Zd\n", Yq);
  gmp_printf(" Zq = %Zd\n", Zq);
	
	gmp_printf(" a = %Zd\n", a);
	gmp_printf(" b = %Zd\n", b);
	gmp_printf(" xd = %Zd\n", xd);

	/* LLamada a RADIX */
  radix(binary, number, foot, bits);
  //gmp_printf("\n binary = %Zd\n",binary);
	
	/* llamada a num_bits para calcular el tamaño del arreglo */
	num_bits(bits, number);
  
	gmp_printf("\n Number       = %Zd\n",number); 	/* numero decimal empleado */
	gmp_printf(" Radix Number = %Zd\n",binary); 	/* numero binario convertido */
	gmp_printf(" Size_Array   = %Zd\n",bits); 	/* tamaño en bits a ocupar en el arreglo */

	/* Llamada a bin_to_vector */

	slots = mpz_get_ui(bits);
	int vector_bin[slots];
	
	bin_to_vector(vector_bin, slots, bits, binary);


	printf("\n\n INICIO MONTGOMERY LADDER\n");

	/*  Ejecución de MONTGOMERY_LADDER */
	time(start);
	
	/* Transformación a Coordenadas Afines */
	mpz_invert(Zp,Zp,P);
	mpz_mul(Xp,Xp,Zp);
	mpz_mod(Xp,Xp,P);
	mpz_cdiv_q(Zp,Zp,Zp);

	mpz_invert(Zq,Zq,P);
	mpz_mul(Xq,Xq,Zq);
	mpz_mod(Xq,Xq,P);
	mpz_cdiv_q(Zq,Zq,Zq);

	gmp_printf("\n Xp= %Zd\n",Xp);
	//gmp_printf(" Yp= %Zd\n",Yp);
	gmp_printf(" Zp= %Zd\n",Zp);
	gmp_printf(" Xq= %Zd\n",Xq);
	//gmp_printf(" Yq= %Zd\n",Yq);
	gmp_printf(" Zq= %Zd\n",Zq);
  
	/*  Inicio  */
	for (int i=1; i<=slots-1; i++){

	 	if(vector_bin[i] == 1){
			R = 1;
  		montgomery_ladder(Xp,Zp,Xq,Zq,P,a,b,xd,R);
		}
		if(vector_bin[i] == 0){
			R = 2;
  		montgomery_ladder(Xp,Zp,Xq,Zq,P,a,b,xd,R); 
		}
	}

	time(end);

	/* Resultados de ECSM */
	gmp_printf("\n Resultado SMEC kP\n");
  gmp_printf(" Xp = %Zd\n",Xp); 
  //gmp_printf(" Yq = %Zd\n",Yp); 
  gmp_printf(" Zp = %Zd\n",Zp); 

	gmp_printf("\n Resultado SMEC k(P+1)\n");
  gmp_printf(" Xq = %Zd\n",Xq); 
  //gmp_printf(" Yq = %Zd\n",Yp); 
  gmp_printf(" Zq = %Zd\n",Zq); 
	
	sec = time_diff(end,start)/TIME;
	printf("\n Tiempo de Ejecución: %.9lf seg.\n\n",sec);

  /* Liberar Memoria  */ 
  mpz_clear (number);
  mpz_clear (foot);
  mpz_clear (binary);
	mpz_clear (bits);
	mpz_clear (P);
	mpz_clear (Xp);
	mpz_clear (Xq);
	mpz_clear (Yp);
	mpz_clear (Yq);
	mpz_clear (Zp);
	mpz_clear (Zq);

  return EXIT_SUCCESS;
}
