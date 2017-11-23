#ifndef TMcmcSimH
#define TMcmcSimH

#include "TSim.h"
#include "TSimDatabase.h"
#include <vector>
//---------------------------------------------------------------------------

class TMcmcSim:public TSim{
	public:
		int nbSims;
		my_string calibrationFile;
		TSimDatabase* myCaliData;
		TOutputVector* myOutputFiles;
		long _Idum;
		double threshold;
		double distance, oldDistance;
		int nbSimulationsPerformed;
		int nbSimulationsAccepted;
		int nbSimulationsWhereHAboveOne;
		int nbSimulationsLargerThreshold;
		int burninLength;
		int burninTrials;
		int repeatsPerParameterVector;
		int burninAttemps;
		float thresholdProportion;
		float mcmc_range_proportion;
		double distanceDensity, oldDistanceDensity;
		bool includePriorDens, includeDistDens, includeDistRatio;
		bool stopIfBurninFailed;

		ofstream hFile;

		TMcmcSim(TParameters* gotParameters, my_string gotexedir, ofstream* gotLogFile);
		TMcmcSim();
		~TMcmcSim(){};

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
