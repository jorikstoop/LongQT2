/*
 * Ten Tusscher's human ventricle model Re-organization in C++.
 * Alternans and spiral breakup in a human ventricular tissue model.
 * 
 * Single cell stimulus protocal suggestion: -52 pA/pF * 1ms.
 * Timestep suggestion: dt = 0.2 ms.
 * Under Intellectual Property Protection.
 * 
 * 
 * Author      	: Shugang Zhang
 * Modified by	: Jorik Stoop, Barton Wei
 * Last update	: Dec-2024
 */

#ifndef _TP06_H_
#define _TP06_H_
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <typeinfo>
#include "omp.h"
#include <ctime>
#include "Cell.cc"
#include "config.h"


/*------------------------------------------------------------------------------
FLAG TO CHOOSE BETWEEN EPICARDIAL ENDOCARDIAL AND MIDMYOCARDIAL CELL TYPES
------------------------------------------------------------------------------*/

using namespace std;


class TP06: public Cell
{
private:
	// define all 20 states 
	double sm     ; //(*V).M
	double sh     ; //(*V).H
	double sj     ; //(*V).J
	double sml	  ; // from ORd model
	double shl	  ; // from ORd model
	double sxs    ; //(*V).Xs
	double ss     ; //(*V).S
	double sr     ; //(*V).R
	double sd     ; //(*V).D
	double sf     ; //(*V).F
	double sf2    ; //(*V).F2
	double sfcass ; //(*V).FCass
	double sRR    ; //(*V).RR
	double sOO    ; //(*V).OO
	double svolt  ; //(*V).Volt
	double Cai    ; //(*V).Cai 
	double CaSR   ; //(*V).CaSR
	double CaSS   ; //(*V).CaSS
	double Nai    ; //(*V).Nai
	double Ki     ; //(*V).Ki

	// Markov IKr
    double C3  ;
    double C2  ;
    double C1  ;
    double O   ;
    double I   ;

	// current
	double Istim,IKr,IKs,IK1,Ito,INa,INaL,IbNa,ICaL,IbCa,INaCa,IpCa,IpK,INaK,Irel,Ileak,Iup,Ixfer;

	// derivatives
	double dvdt, dvgap_dt;

	// CaMK related
	double CaMKII, CaMKII_CaMCa4, CaMKIIP_CaMCa4, OX_CaMCa4, OX, CaMKIIP, CaM, CaMCa, CaMCa2, CaMCa3, CaMCa4;

	// others
	double dt;
	double quin_conc;// umol

	CellType ctype;

public:

	/* CONSTRUCTOR FUNCTION*/
	TP06(CellType ct);

	/* INITIALIZATION FUNCTION*/
	void init(CellType ct);

	/* UPDATE FUNCITON */
	void update();

	/* output all states, mainly used for generating initialization values */
	void outputAllStates(FILE *datafile);

	/* readin all states from initialization files */
	void readinAllStates(FILE *datafile);

	/* GET & SET */
	void setQuin(double param); // umol
	void setIstim(double param); // pA/pF
	void setDt(double param); // ms
	void setDVgap_dt(double param); // mV/ms
	void setVolt(double param);
	double getDVgap_dt(double param); // mV/ms
	double getAbsDvdt(); // mV/ms
	double getIstim(); // pA/pF
	double getICaL();
	double getINa();
	double getINaL();
	double getIKr();
	double getV(); // mV
	//ADDED
	double getC1();
	double getC2();
	double getC3();
	double getI();
	double getO();
	CellType getCellType();


};
#endif