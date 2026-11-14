/*
 * General Code Structure (GCS) for single cell simulation 
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

void run_singlecell_sim(TP06* cell, string fname_data, string fname_apd, 
						double dt, double numS1,
			 			double stopTime, double stimStrength,
			 			double stimDuration, double stimStart)
{
	// statistics for single cell
	double apd20;
	double apd25;
	double apd50;
	double apd75;
	double apd80;
	double apd90;
	double dvdt_max;
	double rest_potential;
	double amplitude;
	double overshoot;

	FILE *datafile = fopen(fname_data.c_str(), "w+");
	FILE *apdfile = fopen(fname_apd.c_str(), "w+");

	// --------start simulation--------
	cell->setDt(dt);

	double time = 0;
	int step = 0;
	for(time = 0.0, step = 0; time < stopTime; time += dt, step++)
	{
		// 1. Apply stimulus according to user configuration
		if(time - floor(time/BCL)*BCL >= stimStart && 
		   time - floor(time/BCL)*BCL < stimStart + stimDuration)	
		{
			cell->setIstim(stimStrength);
		}
		else
		{
			cell->setIstim(0.0);
		}

		// 2. update all states and currents
		cell->update();
		//cout << "Finished updating states" << endl;


		// 3. output file, write file each 1ms 
		if(step%(int(0.02/dt)) == 0 && floor(time/BCL) >= numS1 - 2)
		{
			//fprintf(datafile,"%4.10f\t", fmod(time,BCL)); // if only focus on the last cycle
			fprintf(datafile,"%4.10f\t", time - (numS1 - 2)*BCL ); // if focus on the last two

			// for the affected current
			fprintf(datafile,"%4.10f\t", cell->getV()); // unit: mV
			fprintf(datafile,"%4.10f\t", cell->getIKr()); // unit: 
			fprintf(datafile,"%4.10f\t", cell->getICaL()); // unit: 
			//fprintf(datafile,"%4.10f\t", cell->getIK1()); // unit:

			//Markov states
			fprintf(datafile,"%4.10f\t", cell->getC1());
			fprintf(datafile,"%4.10f\t", cell->getC2());
			fprintf(datafile,"%4.10f\t", cell->getC3());
			fprintf(datafile,"%4.10f\t", cell->getI());
			fprintf(datafile,"%4.10f\t", cell->getO());
			fprintf(datafile,"\n");
		}

		// 4. Statistics of single cell, including APD20, APD50, APD90, dvdt_max, resting potential, overshoot, amplitude
		double oldV, oldDvdt, repo20, repo25, repo50, repo75, repo80, repo90;
		bool peakfound;
		// time-stimStart is used here.
		if (floor((time-stimStart)/BCL) == numS1 - 2 && floor((time+dt-stimStart)/BCL) == numS1 - 1)
		{
			rest_potential = cell->getV();
			oldDvdt = cell->getAbsDvdt();// record this to help detect the maximum dv/dt (UpstrokeVelocity_max)
			oldV = cell->getV();// record this to help detect whether the peak potential (overshoot)
			peakfound = false; // record whether peak potential (overshoot) found or not

			// 20190420 test
			dvdt_max = -1;
		}
		if (floor((time-stimStart)/BCL) == numS1 - 1)
		{
			if(cell->getV() > oldV)
			{
				// 20190420 test
				oldV = cell->getV();

				if( cell->getAbsDvdt() > dvdt_max )
				{
					dvdt_max = cell->getAbsDvdt();
					// cout << cell->getINa() << " " << cell->getIstim() << endl; 
				}
			}	
			else if(cell->getV() <= oldV && !peakfound) // peak not found yet
			{
				peakfound = true;
				overshoot = oldV;
				amplitude = overshoot - rest_potential; // should always be a positive value
				repo20 = overshoot - 0.20*amplitude;
				repo25 = overshoot - 0.25*amplitude;
				repo50 = overshoot - 0.50*amplitude;
				repo75 = overshoot - 0.75*amplitude;
				repo80 = overshoot - 0.80*amplitude;
				repo90 = overshoot - 0.90*amplitude;
				oldV = cell->getV();

			}
			else if(cell->getV() <= oldV && peakfound) // peak already found
			{
				if(oldV >= repo20 && cell->getV() <= repo20)
					apd20 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo25 && cell->getV() <= repo25)
					apd25 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo50 && cell->getV() <= repo50)
					apd50 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo75 && cell->getV() <= repo75)
					apd75 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo80 && cell->getV() <= repo80)
					apd80 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo90 && cell->getV() <= repo90)
					apd90 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				oldV = cell->getV();
			}
		}
			
		// 5. output apd statistics to file
		if (time + dt >= stopTime)
		{	// APD20, APD50, APD90, ,  , , 
			fprintf(apdfile,"Resting membrane potential = %.5f mV\n",rest_potential); // unit: mV
			fprintf(apdfile,"Maximum dV/dt (UpstrokeVelocity_max) = %.5f mV/ms\n",dvdt_max); // unit: mV/ms
			fprintf(apdfile,"Overshoot = %.5f mV\n",overshoot); // unit: mV	
			fprintf(apdfile,"Amplitude = %.5f mV\n",amplitude); // unit: mV
			fprintf(apdfile,"APD20 = %.5f ms\n",apd20); // unit: ms					
			fprintf(apdfile,"APD25 = %.5f ms\n",apd25); // unit: ms	
			fprintf(apdfile,"APD50 = %.5f ms\n",apd50); // unit: ms	
			fprintf(apdfile,"APD75 = %.5f ms\n",apd75); // unit: ms		
			fprintf(apdfile,"APD80 = %.5f ms\n",apd80); // unit: ms		
			fprintf(apdfile,"APD90 = %.5f ms\n",apd90); // unit: ms
		}

		// 6. output apd statistics to screen
		if (time + dt >= stopTime)
		{	// APD20, APD50, APD90, ,  , , 
			printf("Resting membrane potential = %.5f mV\n",rest_potential); // unit: mV
			printf("Maximum dV/dt (UpstrokeVelocity_max) = %.5f mV/ms\n",dvdt_max); // unit: mV/ms
			printf("Overshoot = %.5f mV\n",overshoot); // unit: mV	
			printf("Amplitude = %.5f mV\n",amplitude); // unit: mV
			printf("APD20 = %.5f ms\n",apd20); // unit: ms	
			printf("AP20 = %.5f ms\n",repo20); // unit: ms	
			printf("APD25 = %.5f ms\n",apd25); // unit: ms	
			printf("AP25 = %.5f ms\n",repo25); // unit: ms		
			printf("APD50 = %.5f ms\n",apd50); // unit: ms	
			printf("APD75 = %.5f ms\n",apd75); // unit: ms		
			printf("APD80 = %.5f ms\n",apd80); // unit: ms	
			printf("AP90 = %.5f ms\n",repo90); // unit: ms	
			printf("APD90 = %.5f ms\n",apd90); // unit: ms
		}
	}
	fclose(datafile);
	fclose(apdfile);
}

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// Modify following user-dependent parameters.
	double numS1 = 100; //number of stimuli
	double dt = 0.02; //ms
	double stopTime = numS1*BCL; //ms
	double stimStrength = -53;	//pA/pF
	double stimDuration = 1;	// ms
	double stimStart = 50.0;   // indicates the time point of beginning stimulus in a cycle
	// ----------------

	typedef TP06 CellType;
	typedef TP06* CellPointer;

	CellPointer endo = new CellType(ENDO);
	CellPointer mcell = new CellType(MCELL);
	CellPointer epi = new CellType(EPI);

	// read in stable initial values
	// if the initfile is not available, run the initialization.cc first
	const string tag = result_tag();
	string fname_endo;
	string fname_mcell;
	string fname_epi;
	fname_endo = "init_vals/InitialValues_ENDO_" + tag + ".dat";
	fname_mcell = "init_vals/InitialValues_MCELL_" + tag + ".dat";
	fname_epi = "init_vals/InitialValues_EPI_" + tag + ".dat";

	FILE *initfile_endo = fopen(fname_endo.c_str(), "r");
	FILE *initfile_mcell = fopen(fname_mcell.c_str(), "r");
	FILE *initfile_epi = fopen(fname_epi.c_str(), "r");

	endo->readinAllStates(initfile_endo);
	mcell->readinAllStates(initfile_mcell);
	epi->readinAllStates(initfile_epi);

	fclose(initfile_endo);
	fclose(initfile_mcell);
	fclose(initfile_endo);

	string fname_data_endo, fname_data_mcell, fname_data_epi;
	string fname_apd_endo, fname_apd_mcell, fname_apd_epi;
	// output file names
	fname_data_endo = "output/SingleCellResults_ENDO_" + tag + ".dat";
	fname_data_mcell = "output/SingleCellResults_MCELL_" + tag + ".dat";
	fname_data_epi = "output/SingleCellResults_EPI_" + tag + ".dat";

	fname_apd_endo = "output/SingleCellStats_ENDO_" + tag + ".dat";
	fname_apd_mcell = "output/SingleCellStats_MCELL_" + tag + ".dat";
	fname_apd_epi = "output/SingleCellStats_EPI_" + tag + ".dat";

	// run single cell simulation(s)
	// endo
	run_singlecell_sim(endo, fname_data_endo, fname_apd_endo, dt, numS1,
					   stopTime, stimStrength, stimDuration, stimStart);
	// mcell
	run_singlecell_sim(mcell, fname_data_mcell, fname_apd_mcell, dt, numS1,
					   stopTime, stimStrength, stimDuration, stimStart);
	// epi
	run_singlecell_sim(epi, fname_data_epi, fname_apd_epi, dt, numS1,
					   stopTime, stimStrength, stimDuration, stimStart);

	return 0;
}



