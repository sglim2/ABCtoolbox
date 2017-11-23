//---------------------------------------------------------------------------

#ifndef TCaliDatabaseH
#define TCaliDatabaseH
#include "TDataVector.h"
#include "TInputFile.h"
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
//---------------------------------------------------------------------------

class TSimDatabase{
	public:
		TDataVector* myDataObject;
		TInputFileVector* inputFiles;
		ofstream* logFile;
		int num_sims;
		int num_sims_below_threshold;
		int num_sims_in_DB;

		double threshold;
		double* data_mean;
		double* data_var;
		double** data;
		double* distances;
		double* distances_normalized;
		double distance_mean, distance_std;
		double** priors;
		double** standardizedParameters;
		double* paramMin;
		double* paramMax;
		double** distances_ordered;
		int* simsBelowThreshold;
		int smallestDistCaliSim;
		bool thresholdSet, matrixWithSimsBelowThresholdFilled;
		bool distancesCalculated;

		TSimDatabase (int gotNum_cali_sims, TDataVector* gotDataPointer, TInputFileVector* gotinputFiles, ofstream* gotLogFile);
		TSimDatabase(my_string fileName, int gotNum_cali_sims,  TDataVector* gotDataPointer, TInputFileVector* gotinputFiles, ofstream* gotLogFile);
		void empty(int newsize);
		void addSimToDB(TDataVector* gotDataPointer, TInputFileVector* gotinputFiles, double distance=-1);
		void calculateMeanVariances();
		void calculateMinMaxofParameters(double* & min, double* & max);
		void standardizeRetainedParameters(double* min, double* max);
		void calculateWeightedSigmaOfRetainedParameters(SymmetricMatrix & Sigma, double* weights);
		double getPriorDensityOneRetainedSim(int sim);
		void calculateDistances();
		void calculateDistances(TLinearComb* linearComb);
		void fillNormalizedDistances();
		double getThreshold(double thresholdProportion);
		void setThreshold(int numToRetain);
		void writeToFile(my_string fileName);
		void writeDistanceFile(my_string fileName);
		double getDistanceDensity(double dist);
		void setPriorStartingConditionsFromBestSimulation();
		void setPriorStartingConditionsAtRandom();
		void setMcmcRanges(float mcmc_range_proportion);
		void setMcmcRangesOnePrior(int prior, float mcmc_range_proportion);
		void fillArrayOfSimsBelowThreshold();
};
#endif
