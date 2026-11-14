/*
 * General Code Structure (GCS) for obtaining steady-state single cell initial values
 * 
 * Under Intellectual Property Protection.
 * 
 * 
 * Author		: Shugang Zhang
 * Modified by	: Jorik Stoop, Barton Wei
 * Last update	: Nov-2026
 */


#include "TP06.h"
#include "config.h"

using namespace std;


void run_sim(TP06* cell, string fname, double dt,
			 double stopTime, double stimStrength,
			 double stimDuration, double stimStart)
{
	FILE *initfile = fopen(fname.c_str(),"w+");
	cell->setDt(dt);
	
	// --------start simulation--------
	double time = 0;
	int step = 0;
	for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
	{
		// 1. progress stats
		if(step%4000000 == 0) // 300000 * dt ms = 4.5s    0.005ms*4000000 = 20s
			cout << "Generating initial values, current progress = " << 100.0*time/stopTime << "\%." << endl;

		// Apply stimulus according to user configuration
		if(time - floor(time/BCL)*BCL >= stimStart && 
		   time - floor(time/BCL)*BCL < stimStart + stimDuration)	
		{
			cell->setIstim(stimStrength);
		}
		else
		{
			cell->setIstim(0.0);
		}

		// update all states and currents
		cell->update();

		// write file for each 1ms 
		if(time + dt > stopTime) // 50*dt = 1ms once
		{
			cell->outputAllStates(initfile);
			cout <<"Initial values of all states have been generated!" << endl;
		}
	}
	fclose(initfile);
}

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// Modify following user-dependent parameters.
	double dt = 0.02; //ms
	double stopTime = 200*BCL; //ms
	double stimStrength = -53; // pA/pF
	double stimDuration = 1.0;	// ms
	double stimStart = 0.0; // ms  // indicates the time point of beginning stimulus in a cycle
	// ----------- end user configuration ---------

	// create file names
	const string tag = result_tag();
	string fname_endo;
	string fname_mcell;
	string fname_epi;
	fname_endo = "init_vals/InitialValues_ENDO_" + tag + ".dat";
	fname_mcell = "init_vals/InitialValues_MCELL_" + tag + ".dat";
	fname_epi = "init_vals/InitialValues_EPI_" + tag + ".dat";

	typedef TP06 CellType;
	typedef TP06* CellPointer;
	CellPointer endo = new CellType(ENDO);
	CellPointer mcell = new CellType(MCELL);
	CellPointer epi = new CellType(EPI);

	// run initialization simulation for each cell type
	run_sim(endo, fname_endo, dt, stopTime, stimStrength, stimDuration, stimStart);
	run_sim(mcell, fname_mcell, dt, stopTime, stimStrength, stimDuration, stimStart);
	run_sim(epi, fname_epi, dt, stopTime, stimStrength, stimDuration, stimStart);

	return 0;
}




