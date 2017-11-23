#ifndef TPmcSimH
#define TPmcSimH

#include "TSim.h"
#include "TSimDatabase.h"
#include <vector>
//---------------------------------------------------------------------------

class TPmcSim:public TSim{
	public:
		my_string calibrationFile;
		TSimDatabase** mySimData;
		TOutputVector* myOutputFiles;
		long _Idum;
		int numInterations;
		int particleSize;
		int lastParticleSize;
		int repeatsPerParameterVector;
		int nextDb, currentDb, oldDb;
		double pmc_range_proportion;
		double tolerance;
		int numRetained;
		double** weights;
		double* paramMin;
		double* paramMax;
		SymmetricMatrix Sigma, Sigma_inv;
		Matrix A;
		int numCaliSims;

		TPmcSim(TParameters* gotParameters, my_string gotexedir, ofstream* gotLogFile);
		TPmcSim();
		~TPmcSim(){};

		virtual int runSimulations();
		bool calculcateH();
		double calculcatePi(TPriorVector Prior);
		double calculatePriorProbability(TPrior Prior);
        // double calculateDistance(TDataVector* observed, TDataVector* simulated);
		double chiSquare(double obs, double exp);
		void performCalibratingSimulations();
		void initialCalibration();
		void resetCounters();
		void performSimulation(int s=-1);
		void performSimulationSeveralRepeats(int numRepeatsPerParameterVector, int s=-1);
		void checkAcceptance();
		void performBurnin();
};
#endif
