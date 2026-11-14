/*
 * General Code Structure (GCS) for one-dimensional simulation
 * 
 * Under Intellectual Property Protection.
 * 
 * 
 * Author		: Shugang Zhang
 * Modified by	: Jorik Stoop, Barton Wei
 * Last update	: Nov-2026
 */


#include "TP06.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "omp.h"
#include "config.h"

using namespace std;

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// 1D parameters for all cell types
	double dx = 0.15; // mm
	double dt = 0.02; // ms
	int sanCellNum = 0;
	int atrialCellNum = 0;
	int epiCellNum = 0;
	int mCellNum = 0;
	int endoCellNum = 0;
	
	// Define number of cells in each layer
	endoCellNum = 74;
	mCellNum = 52;
	epiCellNum = 74;

	//Define diffusion coefficient (mm^2/ms)
	double ventCoeff = 0.154; 
	int numS1 = 15; // number of stimulus before last stim
	double stopTime = numS1*BCL + BCL; //ms
	double stimStrength = -52.0; // pA
	double stimDuration = 1.0;	// ms
	double stimStart = 0.0; // ms  // indicates the time point of beginning stimulus in a cycle
	// -------- end user configuration --------


	// --------start simulation--------	
	// CV calculation
	double cvStartTime = 0;
	double cvEndTime = 0;
	int cvStartFlag = 0; // 0 for not start yet
	int cvEndFlag = 0;
	double cv = 0;

	// parallel
	int coreNum = 8;
	omp_set_num_threads(coreNum);

	// strand initilization, diffusion params
	int cellNum = sanCellNum + atrialCellNum + epiCellNum + mCellNum + endoCellNum; // number of cells in OneD strand
	typedef Cell* CellPointer;
	TP06* strand[cellNum];		// note that constructor contains the initializer
	double coeff[cellNum];		// diffusion parameters for each cell
	double dcoeff_dx[cellNum];	// first order derivative for each cell
	double oldV[cellNum];

	// assign coeff according to cell type
	for(int i = 0; i < cellNum; i++)
	{
		// set coeff to 'ventCoeff' or specify decoupling
		if(i == endoCellNum+mCellNum)
			coeff[i] = 0.2*ventCoeff;
		else
			coeff[i] = ventCoeff;
	}

	// Calculate the dcoeff/dx(i.e. dcoeff_dx in the code) in the 1D strand
	for(int i = 0; i < cellNum; i++)
	{
		if (i == 0) 
			dcoeff_dx[i] = (coeff[i+1] - coeff[i])/dx;
		else if (i == cellNum-1) 
			dcoeff_dx[i] = (coeff[i] - coeff[i-1])/dx;
		else
			dcoeff_dx[i] = (coeff[i+1] - coeff[i-1])/(2.0*dx);
	}


	const string tag = result_tag();
	string fname_endo;
	string fname_mcell;
	string fname_epi;
	fname_endo = "init_vals/InitialValues_ENDO_" + tag + ".dat";
	fname_mcell = "init_vals/InitialValues_MCELL_" + tag + ".dat";
	fname_epi = "init_vals/InitialValues_EPI_" + tag + ".dat";


	#pragma omp parallel for schedule(static)
	for (int i = 0; i < cellNum; i++)
	{
		FILE *initfile_endo;
		FILE *initfile_mcell;
		FILE *initfile_epi;
		// read in stable initial values
		// if the initfile is not available, run Initialization.cc first
		if(i >= 0 && i < endoCellNum)
		{
			strand[i] = new TP06(ENDO);
			initfile_endo = fopen(fname_endo.c_str(), "r");
			strand[i]->readinAllStates(initfile_endo);
			fclose(initfile_endo);
		}
		else if (i < endoCellNum + mCellNum)
		{
			strand[i] = new TP06(MCELL);
			initfile_mcell = fopen(fname_mcell.c_str(), "r");
			strand[i]->readinAllStates(initfile_mcell);
			fclose(initfile_mcell);
		}
		else
		{
			strand[i] = new TP06(EPI);
			// read in stable initial values
			// if the initfile is not available, run the initialization.cc first
			initfile_epi = fopen(fname_epi.c_str(), "r");
			strand[i]->readinAllStates(initfile_epi);
			fclose(initfile_epi);
		}	

		strand[i]->setDt(dt);
	}


	string fname_data = "output/OneDResults_" + tag + ".dat";
	string fname_sodium = "output/OneDsodium_" + tag + ".dat";
	FILE *datafile = fopen(fname_data.c_str(),"w+");
	FILE *sodiumfile = fopen(fname_sodium.c_str(),"w+");


	double time = 0;
	int step = 0;
	for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
	{
		// 1. progress stats
		if(step%25000 == 0) // 25000 * dt ms = 0.125s 
			cout << "Progress = " << 100.0*time/stopTime << "\%." << endl;

		for(int i = 0; i < cellNum; i++)
		{
			oldV[i] = strand[i]->getV();
		}

		
		// 2. cell loop
		// #pragma omp parallel for schedule(static)
		for(int i = 0; i < cellNum; i++)
		{
			// ---stimulation applying or not, based on the time and cell location---
			// initialize all stim to 0.0
			strand[i]->setIstim(0.0);
			// Apply stimulus according to user configuration (about time)
			if(time - floor(time/BCL)*BCL >= stimStart && 
		   	   time - floor(time/BCL)*BCL < stimStart + stimDuration)
			{// first 5 cells get S1 stimulation
		    	if(i < 3 && i >= 0)
				{// cells get stimulation in certain duration
					// as strand[i] is a pointer, A->B should be used instead of A.B
					// ref. http://blog.csdn.net/qq457163027/article/details/54237782
					// cout << "here" << endl;
					strand[i]->setIstim(stimStrength);
				}
			}

			// ---------calculate diffusion, i.e. dVgap---------
		
			double dVgap_dt = 0;
			double first_order;
			double second_order;

			// Step 1: calculate first and second order of membrane potential
			if(i == 0) 
			{
				// use strand[0] instead of "strand[-1]"
				first_order = (oldV[i+1] - oldV[i])/(1.0*dx);
				second_order = (oldV[i+1] + oldV[i] - 2.0*oldV[i])/(dx*dx);
			}
			else if(i > 0 && i < cellNum - 1) 
			{
				// normal case
				first_order = (oldV[i+1] - oldV[i-1])/(2.0*dx);
				second_order = (oldV[i+1] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
			}
			else if(i == cellNum - 1)
			{
				// use strand[cellNum-1] instead of "strand[cellNum]" as the latter is out of index
				first_order = (oldV[i] - oldV[i-1])/(1.0*dx);
				second_order = (oldV[i] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
			}

			// Step 2: calculate dVgap according to equations
			dVgap_dt = dcoeff_dx[i]*first_order + coeff[i]*second_order;

			//if(step%1000 == 0) cout<<strand[i]->getVgap()<<endl;
			strand[i]->setDVgap_dt(dVgap_dt);
				
			// update
			strand[i]->update();
		}// end cell loop

		// 3. output file. Unfortunately, this part cannot run in parallel
		// if( floor(time/BCL) >= numS1 - 20) // output final cycle only
		if(step%10 == 0) // 50*dt = 1 ms once
		{
			for(int j = 0; j < cellNum; j++)
			{
				// write time before each time step
				if(j == 0){
					fprintf(datafile,"%4.10f\t", time);
					fprintf(sodiumfile,"%4.10f\t", time);
				}
				// write volt for each cell
				fprintf(datafile,"%4.10f\t", strand[j]->getV()); // unit: mV
				// write INa for each cell
				fprintf(sodiumfile,"%4.10f\t", strand[j]->getINa()); // unit: mV
				// write enter after each time step
				if(j == cellNum - 1){
					fprintf(datafile,"\n");
					fprintf(sodiumfile,"\n");
				}
			}

		}// end Membrane Potential recording

		// 4. CV stats. Note that this can be only applied to non-selfpacing cells!!!
		// only record the last S1 data
		if (floor(time/BCL) == numS1)
		{ 
			// record propagation start and end time for CV calculation
			if(strand[10]->getV() >= -30 && cvStartFlag == 0)
			{
				cvStartTime = time;
				cout << "start = " << cvStartTime << endl;
				cvStartFlag = 1;
			}
			if(strand[90]->getV() >= -30 && cvEndFlag == 0)
			{
				cvEndTime = time;
				cout << "end = " << cvEndTime << endl;
				cvEndFlag = 1;
				cv = (dx * 80) / (cvEndTime - cvStartTime);
				cout << "duration = " << cvEndTime - cvStartTime << endl;
			}
		} // end CV calculation

	}// end of timeloop
	fclose(datafile);
	fclose(sodiumfile);

	if(cvStartFlag == 1 && cvEndFlag == 1)
		cout << "CV = " << cv << " m/s." << endl;
	else
		cout << "Conduction failure!" << endl;
	printf("All done.\n");

	return 0;
}