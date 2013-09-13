//==============================================================================
// Name        : RBC_CPP.cpp
// Description : Basic RBC model with partial depreciation and endogenous labor
// Date        : July 21, 2013
// Modified    : Sept 12, 2013
//==============================================================================

#include <iostream>
#include <math.h>      // power
#include <cmath>       // abs
#include <ctime>       // time

using namespace std;

int main() {

  clock_t begin = clock();

  ///////////////////////////////////////////////////////////////////////////////////////////
  // 1. Calibration
  ///////////////////////////////////////////////////////////////////////////////////////////

  const double aalpha = 0.33333333333;     // Elasticity of output w.r.t. capital
  const double bbeta  = 0.95;              // Discount factor;
  const double ddelta = 0.09;			   // Rate of Depreciation

  // Productivity values

  double vProductivity[5] ={0.9792, 0.9896, 1.0000, 1.0106, 1.0212};

  // Transition matrix
  double mTransition[5][5] = {
			{0.9727, 0.0273, 0.0000, 0.0000, 0.0000},
			{0.0041, 0.9806, 0.0153, 0.0000, 0.0000},
			{0.0000, 0.0082, 0.9837, 0.0082, 0.0000},
			{0.0000, 0.0000, 0.0153, 0.9806, 0.0041},
			{0.0000, 0.0000, 0.0000, 0.0273, 0.9727}
			};

  ///////////////////////////////////////////////////////////////////////////////////////////
  // 2a. Steady State
  ///////////////////////////////////////////////////////////////////////////////////////////
  double laborSteadyState = 0.33333333333;     // Elasticity of output w.r.t. capital
  double capitalSteadyState = laborSteadyState*pow(((1/bbeta)+ddelta-1)/aalpha,1/(aalpha-1));
  double outputSteadyState  = pow(capitalSteadyState,aalpha)*pow(laborSteadyState,1-aalpha);
  double consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState;

  ///////////////////////////////////////////////////////////////////////////////////////////
  // 2a. Calibrate Disutility of Labor
  ///////////////////////////////////////////////////////////////////////////////////////////
  double ppsi = (1/consumptionSteadyState)*3*(1-aalpha)*pow((((1/bbeta)+ ddelta - 1)/aalpha),(aalpha/(aalpha-1))); //disutility of labor - added 9/12/13

  cout <<"Calibration";
  cout <<"\n";
  cout <<"Psi = " <<ppsi<<"\n";
  cout <<"Steady State Values"<<"\n";
  cout <<"Output = "<<outputSteadyState<<", Capital = "<<capitalSteadyState<<", Consumption = "<<consumptionSteadyState<<"\n";
  cout <<" ";

  // We generate the grid of capital
  int nCapital, nCapitalNextPeriod, gridCapitalNextPeriod, nProductivity, nProductivityNextPeriod;
  const int nGridCapital = 1191, nGridProductivity = 5;
  double vGridCapital[nGridCapital];

  // Now generate the grid for labor
  int nLabor;
  const int nGridLabor = 334;
  double vGridLabor[nGridLabor];

  for (nCapital = 0; nCapital < nGridCapital; ++nCapital){
    vGridCapital[nCapital] = 0.5*capitalSteadyState+0.001*nCapital; //TODO Change Grid to 0.00001
  }

  for (nLabor = 0; nLabor < nGridLabor; ++nLabor){
    vGridLabor[nLabor] = 0.5*laborSteadyState+0.001*nLabor; //TODO Change Grid to 0.00001
  }


  // 3. Required matrices and vectors

  double mLabor[nGridCapital][nGridProductivity];
  double mOutput[nGridCapital][nGridProductivity][nGridLabor];
  double mValueFunction[nGridCapital][nGridProductivity];
  double mValueFunctionNew[nGridCapital][nGridProductivity];
  double mPolicyFunction[nGridCapital][nGridProductivity];
  double expectedValueFunction[nGridCapital][nGridProductivity];

  // 4. We pre-build output for each point in the grid

  for (nProductivity = 0; nProductivity<nGridProductivity; ++nProductivity){
    for (nCapital = 0; nCapital < nGridCapital; ++nCapital){
    	for (nLabor = 0; nLabor < nGridLabor; ++nLabor) {
    		mOutput[nCapital][nProductivity][nLabor] = vProductivity[nProductivity]*pow(vGridCapital[nCapital],aalpha)*pow(vGridLabor[nLabor],1-aalpha);
    	}
    }
  }

  // 5. Main iteration

  double maxDifference = 10.0, diff, diffHighSoFar;
  double tolerance = 0.0000001;
  double valueHighSoFar, valueProvisional, consumption, capitalChoice, laborChoice;

  int iteration = 0;

  while (maxDifference>tolerance){

	//loop to fill in expected value function
    for (nProductivity = 0;nProductivity<nGridProductivity;++nProductivity){
      for (nCapital = 0;nCapital<nGridCapital;++nCapital){
	expectedValueFunction[nCapital][nProductivity] = 0.0;
		for (nProductivityNextPeriod = 0;nProductivityNextPeriod<nGridProductivity;++nProductivityNextPeriod){
			expectedValueFunction[nCapital][nProductivity] += mTransition[nProductivity][nProductivityNextPeriod]*mValueFunction[nCapital][nProductivityNextPeriod];
		}
      }
    }

    //value function iteration loop
    for (nProductivity = 0;nProductivity<nGridProductivity;++nProductivity){

      // We start from previous choice (monotonicity of policy function)
      gridCapitalNextPeriod = 0;

      for (nCapital = 0;nCapital<nGridCapital;++nCapital){

	valueHighSoFar = -1000.0;
	capitalChoice  = vGridCapital[0];

	for (nCapitalNextPeriod = gridCapitalNextPeriod;nCapitalNextPeriod<nGridCapital;++nCapitalNextPeriod){

		double valueHighSoFar2 = -1000.0;
		int gridLabor = 0;

				for (nLabor = 0;nLabor<nGridLabor;++nLabor){

					consumption = mOutput[nCapital][nProductivity][nLabor]-vGridCapital[nCapitalNextPeriod]+(1-ddelta)*vGridCapital[nCapital];
					valueProvisional = (1-bbeta)*log(consumption)-(ppsi/2)*pow(vGridLabor[nLabor],2)+bbeta*expectedValueFunction[nCapitalNextPeriod][nProductivity];


					if (valueProvisional>valueHighSoFar2){
					    valueHighSoFar2 = valueProvisional;
					    laborChoice = vGridLabor[nLabor];
					    gridLabor = nLabor;
						}
					 else{
					      break; //We break when we have achieved the max
					 }

				}

	  if (valueProvisional>valueHighSoFar){
	    valueHighSoFar = valueProvisional;
	    capitalChoice = vGridCapital[nCapitalNextPeriod];
	    gridCapitalNextPeriod = nCapitalNextPeriod;
	  }
	  else{
	    break; // We break when we have achieved the max
	  }

      }
		  mValueFunctionNew[nCapital][nProductivity] = valueHighSoFar;
		  mPolicyFunction[nCapital][nProductivity] = capitalChoice;
		  mLabor[nCapital][nProductivity] = laborChoice;

      }
    }

    diffHighSoFar = -100000.0;
    for (nProductivity = 0;nProductivity<nGridProductivity;++nProductivity){
      for (nCapital = 0;nCapital<nGridCapital;++nCapital){
	diff = std::abs(mValueFunction[nCapital][nProductivity]-mValueFunctionNew[nCapital][nProductivity]);
	if (diff>diffHighSoFar){
	  diffHighSoFar = diff;
	}
	mValueFunction[nCapital][nProductivity] = mValueFunctionNew [nCapital][nProductivity];
      }
    }
    maxDifference = diffHighSoFar;

    iteration = iteration+1;
    if (iteration % 10 == 0 || iteration ==1){
      cout <<"Iteration = "<<iteration<<", Sup Diff = "<<maxDifference<<"\n";
    }
  }

  cout <<"Iteration = "<<iteration<<", Sup Diff = "<<maxDifference<<"\n";
  cout <<" \n";
  cout <<"My check = "<< mPolicyFunction[0][0]<<"\n";
  cout <<" \n";


  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout <<"Elapsed time is "<<elapsed_secs<<" seconds.";
  cout <<" \n";  

  return 0;

}
