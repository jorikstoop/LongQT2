/*
 * Ten Tusscher's human ventricle model in C++.
 * Alternans and spiral breakup in a human ventricular tissue model.
 * 
 * Single cell stimulus protocal suggestion: -52 pA/pF * 1ms.
 * Timestep suggestion: dt = 0.2 ms.
 * Under Intellectual Property Protection.
 * 
 * 
 * Author		: Shugang Zhang
 * Modified by	: Jorik Stoop, Barton Wei
 * Last update	: Nov-2026
 */


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
#include "TP06.h"
#include "config.h"

using namespace std;

// Description can be found in TP06.h
TP06::TP06(CellType ct)
{
	init(ct);
}

void TP06::init(CellType ct)
{
	//Initial values of state variables
	svolt = -86.2;
	Cai = 0.00007;
	CaSR = 1.3;
	CaSS = 0.00007;
	Nai = 7.67;
	Ki = 138.3;
	sm = 0.;
	sh = 0.75;
	sj = 0.75; 
	sml = 0;
	shl = 0;
	sxs = 0.;
	sr = 0.;
	ss = 1.;
	sd = 0.;
	sf = 1.;
	sf2 =1.;
	sfcass = 1.;                                        
	sRR = 1.;  
	sOO = 0.;

	// Markov IKr
	C1  = 0;
	C2  = 0;
	C3  = 1;
	O   = 0;
	I   = 0;

	ctype = ct;
}

void TP06::outputAllStates(FILE *datafile)
{
	fprintf(datafile,"%4.10f\n", svolt   );
    fprintf(datafile,"%4.10f\n", Cai     );
    fprintf(datafile,"%4.10f\n", CaSR    );
    fprintf(datafile,"%4.10f\n", CaSS    );
    fprintf(datafile,"%4.10f\n", Nai     );
    fprintf(datafile,"%4.10f\n", Ki      );
    fprintf(datafile,"%4.10f\n", sm      );
    fprintf(datafile,"%4.10f\n", sh      );
    fprintf(datafile,"%4.10f\n", sml     );// late sodium channel gates from ORd.
    fprintf(datafile,"%4.10f\n", shl     );// late sodium channel gates from ORd.
    fprintf(datafile,"%4.10f\n", sj      );
    fprintf(datafile,"%4.10f\n", sxs     );
    fprintf(datafile,"%4.10f\n", sr      );
    fprintf(datafile,"%4.10f\n", ss      );
    fprintf(datafile,"%4.10f\n", sd      );
    fprintf(datafile,"%4.10f\n", sf      );
    fprintf(datafile,"%4.10f\n", sf2     );
    fprintf(datafile,"%4.10f\n", sfcass  );
    fprintf(datafile,"%4.10f\n", sRR     );
    fprintf(datafile,"%4.10f\n", sOO     );
    
	// IKr Markov
    fprintf(datafile,"%4.10f\n", C3      );
    fprintf(datafile,"%4.10f\n", C2      );
    fprintf(datafile,"%4.10f\n", C1      );
    fprintf(datafile,"%4.10f\n", O       );
    fprintf(datafile,"%4.10f\n", I       );
}

void TP06::readinAllStates(FILE *datafile)
{
	double value;
	fscanf(datafile,"%lf", &value); 
	svolt  = value;
	fscanf(datafile,"%lf", &value); 
	Cai    = value;
	fscanf(datafile,"%lf", &value); 
	CaSR   = value;
	fscanf(datafile,"%lf", &value); 
	CaSS   = value;
	fscanf(datafile,"%lf", &value); 
	Nai    = value;
	fscanf(datafile,"%lf", &value); 
	Ki     = value;
	fscanf(datafile,"%lf", &value); 
	sm     = value;
	fscanf(datafile,"%lf", &value); 
	sh     = value;
	fscanf(datafile,"%lf", &value); // late sodium channel from ORd model
	sml    = value;
	fscanf(datafile,"%lf", &value); // late sodium channel from ORd model
	shl    = value;
	fscanf(datafile,"%lf", &value); 
	sj     = value;
	fscanf(datafile,"%lf", &value); 
	sxs    = value;
	fscanf(datafile,"%lf", &value); 
	sr     = value;
	fscanf(datafile,"%lf", &value); 
	ss     = value;
	fscanf(datafile,"%lf", &value); 
	sd     = value;
	fscanf(datafile,"%lf", &value); 
	sf     = value;
	fscanf(datafile,"%lf", &value); 
	sf2    = value;
	fscanf(datafile,"%lf", &value); 
	sfcass = value;
	fscanf(datafile,"%lf", &value); 
	sRR    = value;
	fscanf(datafile,"%lf", &value); 
	sOO    = value;

	// IKr Markov
	fscanf(datafile,"%lf", &value); 
	C3     = value;
	fscanf(datafile,"%lf", &value); 
	C2     = value;
	fscanf(datafile,"%lf", &value); 
	C1     = value;
	fscanf(datafile,"%lf", &value); 
	O     = value;
	fscanf(datafile,"%lf", &value); 
	I     = value;
}

void TP06::update()
{
	/*-----------------------------------------------------------------------------
	ELECTROPHYSIOLOGICAL PARAMETERS:
	-----------------------------------------------------------------------------*/

	//External concentrations
	double Ko = 5.4;

#ifdef K_supplement
	// ref = https://doi.org/10.1016/j.jacc.2003.07.006
	Ko = 5.4*1.3;
#endif
	double Cao = 2.0;
	double Nao = 140.0;

	//Intracellular volumes
	double Vc = 0.016404;
	double Vsr = 0.001094;
	double Vss = 0.00005468;

	//Calcium buffering dynamics
	double Bufc = 0.2;
	double Kbufc = 0.001;
	double Bufsr = 10.;
	double Kbufsr = 0.3;
	double Bufss = 0.4;
	double Kbufss = 0.00025;

	//Intracellular calcium flux dynamics
	double Vmaxup = 0.006375;
#ifdef EAD
	Vmaxup = 0.6*Vmaxup; //MODIFIED -40%
#endif
	double Kup = 0.00025;
	double Vrel = 0.102;
	double k1_ = 0.15;
	double k2_ = 0.045;
	double k3 = 0.060;
	double k4 = 0.005;
	double EC = 1.5;
	double maxsr = 2.5;
	double minsr = 1.;
	double Vleak = 0.00036;
	double Vxfer = 0.0038;

	//Constants
	double R = 8314.472;
	double F = 96485.3415;
	double T = 310.0;
	double RTONF = (R*T)/F;
	double amp_ikr = 1.0;
	double amp_iks = 1.0;
	double amp_ica = 1.0;
	double amp_ina = 1.0;

	//Cellular capacitance         
	double CAPACITANCE = 0.185;

	//Parameters for currents
	//Parameters for IKr
	double Gkr = 0.0135*pow(Ko,0.59);	// from Clancy Rudy 2001
#ifdef T474I
	Gkr = 0.65*Gkr;	// from Clancy Rudy 2001, 35% reduction compare to WT
#endif
#ifdef EAD
	Gkr = 1.5*Gkr;	// MODIFIED +50%
#endif
	if(ctype == EPI)
		Gkr = 1.2*Gkr;

	//Parameters for Iks
	double pKNa = 0.03;
	double Gks;
	if(ctype == MCELL){
		Gks = 0.098;
	} else{
		// EPI & ENDO
		Gks = 0.392;
	}
#ifdef EAD
	Gks = 0.4*Gks;	// MODIFIED -60% FROM ORIGINAL
#endif
	//Parameters for Ik1
	double GK1 = 5.405;
	
  	//Parameters for Ito
	double Gto;
	if(ctype == ENDO) 
		Gto = 0.073;
	else // EPI & MCELL
		Gto = 0.294;


  	//Parameters for INa
	double GNa = 14.838;
	
  	// INaL from ORd model
	double GNaL = 0.0075;
	if (ctype == EPI)  
		GNaL*=0.6;
	
  	//Parameters for IbNa
	double GbNa = 0.00029;
	
  	//Parameters for INaK
	double KmK = 1.0;
	double KmNa = 40.0;
	double knak = 2.724;
	
  	//Parameters for ICaL
	double GCaL = 0.00003980;
#ifdef EAD
	GCaL = 2.2*GCaL; // MODIFIED
#else
	
  	//Parameters for IbCa
	double GbCa = 0.000592;
	
  	//Parameters for INaCa
	double knaca = 1000;
	double KmNai = 87.5;
	double KmCa = 1.38;
	double ksat = 0.1;
	double n = 0.35;
	
  	//Parameters for IpCa
	double GpCa = 0.1238;
	double KpCa = 0.0005;
	
  	//Parameters for IpK;
	double GpK = 0.0146;

	// current	
	double k1;
	double k2;
	double kCaSR;

	// derivatives
	double dNai,dKi,dCai,dCaSR,dCaSS,dRR,dbadt,dbidt,dbrdt;

 	// reverse potential
	double Ek,Ena,Eks,Eca;

	double CaCSQN;
	double bjsr;
	double cjsr;
	double CaSSBuf;
	double bcss;
	double ccss;
	double CaBuf;
	double bc;
	double cc;
	double Ak1;
	double Bk1;
	double rec_iK1;
	double rec_ipK;
	double rec_iNaK;
	double AM;
	double BM;
	double AH_1;
	double BH_1;
	double AH_2;
	double BH_2;
	double AJ_1;
	double BJ_1;
	double AJ_2;
	double BJ_2;
	double M_INF;
	double H_INF;
	double ML_INF;// late sodium channel gate variables from ORd model
	double HL_INF;// late sodium channel gate variables from ORd model
	double J_INF;
	double TAU_M;
	double TAU_H;
	double TAU_J;
	double TAU_ML;// late sodium channel gate variables from ORd model
	double TAU_HL;// late sodium channel gate variables from ORd model
	double axr1;
	double bxr1;
	double axr2;
	double bxr2;
	double Xr1_INF;
	double Xr2_INF;
	double TAU_Xr1;
	double TAU_Xr2;
	double Axs;
	double Bxs;
	double Xs_INF;
	double TAU_Xs;
	double R_INF;
	double TAU_R;
	double S_INF;
	double TAU_S;
	double Ad;
	double Bd;
	double Cd;
	double Af;
	double Bf;
	double Cf;
	double Af2;
	double Bf2;
	double Cf2;
	double TAU_D;
	double D_INF;
	double TAU_F;
	double F_INF;
	double TAU_F2;
	double F2_INF;
	double TAU_FCaSS;
	double FCaSS_INF;

	
	static double inverseVcF2 = 1/(2*Vc*F);
	static double inverseVcF = 1./(Vc*F);
	static double inversevssF2 = 1/(2*Vss*F);

	 
	static char s[200];
	FILE *FF;

#ifdef ISO_10nM
	// ref = https://doi.org/10.1111/j.1469-7793.2000.t01-2-00391.x
	// and https://doi.org/10.1111/j.1476-5381.2012.02015.x
	amp_ikr = 1.8;
	amp_iks = 2.3;
	amp_ica = 1.6;
#endif

#ifdef Verapamil_300nM
	// ref = https://doi.org/10.1007/PL00004921
	amp_ica = 1.17*exp(-9*BCL/10000);
#endif

    //Needed to compute currents
    Ek = RTONF*(log((Ko/Ki)));
    Ena = RTONF*(log((Nao/Nai)));
    Eks = RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    Eca = 0.5*RTONF*(log((Cao/Cai)));
    Ak1 = 0.1/(1.+exp(0.06*(svolt-Ek-200)));
    Bk1 = (3.*exp(0.0002*(svolt-Ek+100))+exp(0.1*(svolt-Ek-10)))/(1.+exp(-0.5*(svolt-Ek)));
    rec_iK1 = Ak1/(Ak1+Bk1);
    rec_iNaK = (1./(1.+0.1245*exp(-0.1*svolt*F/(R*T))+0.0353*exp(-svolt*F/(R*T))));
    rec_ipK = 1./(1.+exp((25-svolt)/5.98));

    //Compute currents
    INa = amp_ina * GNa*sm*sm*sm*sh*sj*(svolt-Ena);
		INaL = GNaL*(svolt-Ena)*sml*shl;     // Late Sodium from ORd model
    ICaL = amp_ica * GCaL*sd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
      (0.25*exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/(exp(2*(svolt-15)*F/(R*T))-1.);
    Ito = Gto*sr*ss*(svolt-Ek);
    IKr = amp_ikr * Gkr*O*(svolt-Ek);// Markov model
    // IKr = Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek); // HH model
    IKs = amp_iks * Gks*sxs*sxs*(svolt-Eks);
    IK1 = GK1*rec_iK1*(svolt-Ek); // different from Lu
    INaCa = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
      (1./(1+ksat*exp((n-1)*svolt*F/(R*T))))*
      (exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
       exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    INaK = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    IpCa = GpCa*Cai/(KpCa+Cai);
    IpK = GpK*rec_ipK*(svolt-Ek);
    IbNa = GbNa*(svolt-Ena);
    IbCa = GbCa*(svolt-Eca);

    //Determine total current pA/pF
    dvdt = -(IKr +
      IKs   +
      IK1   +
      Ito   +
      INa   +
      INaL  + // from ORd model
      IbNa  +
      ICaL  +
      IbCa  +
      INaK  +
      INaCa +
      IpCa  +
      IpK   +
      Istim);


    //update concentrations    
    kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR))); 
    k1=k1_/kCaSR;
    k2=k2_*kCaSR;
    dRR=k4*(1-sRR)-k2*CaSS*sRR;
    sRR+=dt*dRR;
    sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);

    Irel=Vrel*sOO*(CaSR-CaSS);
    Ileak=Vleak*(CaSR-Cai);
    Iup=Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
    Ixfer=Vxfer*(CaSS-Cai);
    
    CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
    dCaSR=dt*(Iup-Irel-Ileak);
    bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
    cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);
    CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
   
    CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);
    dCaSS=dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
    bcss=Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;
    ccss=Kbufss*(CaSSBuf+dCaSS+CaSS);
    CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;

    CaBuf=Bufc*Cai/(Cai+Kbufc);
    dCai=dt*((-(IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
    bc=Bufc-CaBuf-dCai-Cai+Kbufc;
    cc=Kbufc*(CaBuf+dCai+Cai);
    Cai=(sqrt(bc*bc+4*cc)-bc)/2;
        
    dNai = -(INa+INaL+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE; // INaL from ORd model
    Nai += dt*dNai;
    
    dKi = -(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
    Ki += dt*dKi;

    //compute steady state values and time constants 
    AM=1./(1.+exp((-60.-svolt)/5.));
    BM=0.1/(1.+exp((svolt+35.)/5.))+0.10/(1.+exp((svolt-50.)/200.));
    TAU_M=AM*BM;
    M_INF=1./((1.+exp((-56.86-svolt)/9.03))*(1.+exp((-56.86-svolt)/9.03)));
    if (svolt>=-40.)
    {
		AH_1=0.; 
		BH_1=(0.77/(0.13*(1.+exp(-(svolt+10.66)/11.1))));
		TAU_H= 1.0/(AH_1+BH_1);
	  }
    else
    {
		AH_2=(0.057*exp(-(svolt+80.)/6.8));
		BH_2=(2.7*exp(0.079*svolt)+(3.1e5)*exp(0.3485*svolt));
		TAU_H=1.0/(AH_2+BH_2);
    }
    H_INF=1./((1.+exp((svolt+71.55)/7.43))*(1.+exp((svolt+71.55)/7.43)));
    if(svolt>=-40.)
    {
		AJ_1=0.;      
		BJ_1=(0.6*exp((0.057)*svolt)/(1.+exp(-0.1*(svolt+32.))));
		TAU_J= 1.0/(AJ_1+BJ_1);
    }
    else
    {
		AJ_2=(((-2.5428e4)*exp(0.2444*svolt)-(6.948e-6)*
			exp(-0.04391*svolt))*(svolt+37.78)/
		      (1.+exp(0.311*(svolt+79.23))));    
		BJ_2=(0.02424*exp(-0.01052*svolt)/(1.+exp(-0.1378*(svolt+40.14))));
		TAU_J= 1.0/(AJ_2+BJ_2);
    }
    J_INF=H_INF;


    //************************Markov IKr***************************
    double a_1, b_1, a_2, b_2, a_i, b_i, a, b, u;
    double C1_o, C2_o, C3_o, I_o, O_o;
    double C1_n, C2_n, C3_n, I_n, O_n;

    C1_o = C1;
    C2_o = C2;
    C3_o = C3;
    I_o = I;
    O_o = O;

#ifdef WT
	// markov parameter values from Clancy-Rudy 2001
    a_1 = 2.172; //a_in
    b_1 = 1.077; //b_in
    a_2 = 0.0655*exp(0.05547153*(svolt-36)); //aa
    b_2 = 0.0029357*exp(-0.02158*svolt); //bb
    a_i = 0.439*exp(-0.02352*(svolt+25))*(4.5/Ko);
    b_i = 0.656*exp(0.000942*svolt)*pow((4.5/Ko),0.3);
    a = 0.0555*exp(0.05547153*(svolt-12));
    b = 0.002357*exp(-0.036588*svolt);
    u = (a_i*b_2)/b_i;
#endif

#ifdef T474I
	// markov parameter values from Clancy-Rudy 2001
    a_1 = 2.172; //a_in
    b_1 = 1.077; //b_in
    a_2 = 0.0655*exp(0.05547153*(svolt+25)); //mutated aa
    b_2 = 0.0029357*exp(-0.02158*svolt); //bb
    a_i = 0.439*exp(-0.02352*(svolt+25))*(4.5/Ko); 
    b_i = 0.656*exp(0.000942*svolt)*pow((4.5/Ko),0.3);
    a = 0.0555*exp(0.05547153*(svolt+6)); //mutated
    b = 0.002357*exp(-0.036588*svolt);
    u = (a_i*b_2)/b_i;
#endif

#ifdef R56Q
    // markov parameter values from Clancy-Rudy 2001
	a_1 = 2.172; //a_in
    b_1 = 1.077; //b_in
    a_2 = 0.0655*exp(0.05547153*(svolt-36)); //aa
    b_2 = (6.3)*0.0029357*exp(-0.02158*svolt); //mutated bb
    a_i = 0.439*exp(-0.02352*(svolt+25))*(4.5/Ko); 
    b_i = 0.656*exp(0.000942*svolt)*pow((4.5/Ko),0.3);
    a = 0.0555*exp(0.05547153*(svolt-12)); 
    b = (10.5)*0.002357*exp(-0.036588*svolt); //mutated
    u = (a_i*b_2)/b_i;
#endif

    // outward branches forming a coefficient for each state
    const double co_C1 = 1./(1+dt*a);
    const double co_C2 = 1./(1+dt*(b + a_1));
    const double co_C3 = 1./(1+dt*(b_1 + a_2 + a_2));
    const double co_O = 1./(1+dt*(b_2 + b_i));
    const double co_I = 1./(1+dt*(u + a_i));

    double iter = 0;
    double err_sum = 1;
    while (err_sum > 1E-100 && iter < 100)
    {
        C1_n = co_C1 * (C1_o + dt*(C2*b));
        C2_n = co_C2 * (C2_o + dt*(C1*a + C3*b_1));
        C3_n = co_C3 * (C3_o + dt*(C2*a_1 + O*b_2 + I*u));
        O_n  = co_O  * (O_o  + dt*(C3*a_2 + I*a_i));
        I_n  = co_I  * (I_o  + dt*(C3*a_2 + O*b_i));
        err_sum = fabs(C1 - C1_n)+fabs(C2 - C2_n)+fabs(C3 - C3_n )+fabs(O - O_n)+fabs(I - I_n);
        C1 = C1_n;
        C2 = C2_n;
        C3 = C3_n;
        O  = O_n;
        I  = I_n;

        iter++;
        
    }
    C1 = C1_n;
    C2 = C2_n;
    C3 = C3_n;
    O = O_n;
    I  = I_n;

    Xs_INF=1./(1.+exp((-5.-svolt)/14.));
    Axs=(1400./(sqrt(1.+exp((5.-svolt)/6))));
    Bxs=(1./(1.+exp((svolt-35.)/15.)));
    TAU_Xs=Axs*Bxs+80;
    
	if(ctype == EPI)
	{
	    R_INF=1./(1.+exp((20-svolt)/6.));
	    S_INF=1./(1.+exp((svolt+20)/5.));
	    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
	    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
	}
	else if(ctype == ENDO) 
	{
	    R_INF=1./(1.+exp((20-svolt)/6.));
	    S_INF=1./(1.+exp((svolt+28)/5.));
	    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
	    TAU_S=1000.*exp(-(svolt+67)*(svolt+67)/1000.)+8.;
	}
	else // MCELL
	{
			R_INF=1./(1.+exp((20-svolt)/6.));
	    S_INF=1./(1.+exp((svolt+20)/5.));
	    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
	    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
	}



    D_INF=1./(1.+exp((-8-svolt)/7.5)); // original
    Ad=1.4/(1.+exp((-35-svolt)/13))+0.25;
    Bd=1.4/(1.+exp((svolt+5)/5));
    Cd=1./(1.+exp((50-svolt)/20));
    TAU_D=Ad*Bd+Cd;
    F_INF=1./(1.+exp((svolt+20)/7)); // original
    Af=1102.5*exp(-(svolt+27)*(svolt+27)/225);
    Bf=200./(1+exp((13-svolt)/10.));
    Cf=(180./(1+exp((svolt+30)/10)))+20;
    TAU_F=Af+Bf+Cf; //original
#ifdef EAD
	TAU_F = TAU_F/2.2; //ICaL fast inactivation MODIFIED
#endif
    F2_INF=0.67/(1.+exp((svolt+35)/7))+0.33;
    Af2=600*exp(-(svolt+25)*(svolt+25)/170);
    Bf2=31/(1.+exp((25-svolt)/10));
    Cf2=16/(1.+exp((svolt+30)/10));
    TAU_F2=Af2+Bf2+Cf2; //ICaL slow inactivation
    FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
    TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;
    // Late sodium channel gate variables in ORd model
		ML_INF = 1.0/(1.0+exp((-(svolt+42.85))/5.264));
		TAU_ML = 1.0/(6.765*exp((svolt+11.64)/34.77)+8.552*exp(-(svolt+77.42)/5.955));
		HL_INF = 1.0/(1.0+exp((svolt+87.61)/7.488));
		TAU_HL = 200.0;


		//Update gates
		sm = M_INF-(M_INF-sm)*exp(-dt/TAU_M);
		sh = H_INF-(H_INF-sh)*exp(-dt/TAU_H);
		sj = J_INF-(J_INF-sj)*exp(-dt/TAU_J);
		// Late sodium channel gates in ORd model
		sml = ML_INF-(ML_INF-sml)*exp(-dt/TAU_ML);
		shl = HL_INF-(HL_INF-shl)*exp(-dt/TAU_HL);
		sxs = Xs_INF-(Xs_INF-sxs)*exp(-dt/TAU_Xs);
		ss = S_INF-(S_INF-ss)*exp(-dt/TAU_S);
		sr = R_INF-(R_INF-sr)*exp(-dt/TAU_R);
		sd = D_INF-(D_INF-sd)*exp(-dt/TAU_D); 
		sf = F_INF-(F_INF-sf)*exp(-dt/TAU_F); 
		sf2 = F2_INF-(F2_INF-sf2)*exp(-dt/TAU_F2); 
		sfcass = FCaSS_INF-(FCaSS_INF-sfcass)*exp(-dt/TAU_FCaSS);

    //update voltage
    svolt = svolt + dt*(dvdt + dvgap_dt);
}

void TP06::setQuin(double param) {quin_conc = param;} // umol
void TP06::setIstim(double param) {Istim = param;} // pA/pF
void TP06::setVolt(double value) {svolt = value;}
void TP06::setDt(double param) {dt = param;} // ms
void TP06::setDVgap_dt(double param) {dvgap_dt = param;}// mV/ms
double TP06::getDVgap_dt(double param) {return dvgap_dt;} // mV/ms
double TP06::getAbsDvdt() {return fabs(dvdt);} // mV/ms
double TP06::getIstim() {return Istim;} // pA/pF
double TP06::getICaL() {return ICaL;} // pA/pF
double TP06::getINa() {return INa;} // pA/pF
double TP06::getINaL() {return INaL;} // pA/pF
double TP06::getIKr() {return IKr;} // pA/pF
double TP06::getV() {return svolt;} // mV
//ADDED GET MARKOV STATE VALUES
double TP06::getC1() {return C1;}
double TP06::getC2() {return C2;}
double TP06::getC3() {return C3;}
double TP06::getI() {return I;}
double TP06::getO() {return O;}

CellType TP06::getCellType() {return ctype;}