%% Basic RBC model with partial depreciation and endogenous labor

% Original by Jesus Fernandez-Villaverde at Haverford, July 31, 2013
% Modifiied by Kory Kantenga at University of Pennsylvania, Sept 12, 2013

%% 0. Housekeeping

clear all
close all
clc

tic

%%  1. Calibration

aalpha = 1/3;     % Elasticity of output w.r.t. capital
bbeta  = 0.95;    % Discount factor
ddelta = 0.09;    % Rate of depreciation - added 9/10/13

% Productivity values
vProductivity = [0.9792; 0.9896; 1.0000; 1.0106; 1.0212]';

% Transition matrix
mTransition   = [0.9727, 0.0273, 0.0000, 0.0000, 0.0000;
                 0.0041, 0.9806, 0.0153, 0.0000, 0.0000;
                 0.0000, 0.0082, 0.9837, 0.0082, 0.0000;
                 0.0000, 0.0000, 0.0153, 0.9806, 0.0041;
                 0.0000, 0.0000, 0.0000, 0.0273, 0.9727];

%% 2a. Steady State

laborSteadyState = 1/3; %added 9/10/13
capitalSteadyState = laborSteadyState*((((1/bbeta)+ddelta-1)/aalpha)^(1/(aalpha-1))); %modified 9/10/13
outputSteadyState = (capitalSteadyState^aalpha)*(laborSteadyState^(1-aalpha)); %modified 9/10/13
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState; %modified 9/10/13

%% 2b. Calibrate Disutility of Labor
ppsi = (1/consumptionSteadyState)*3*(1-aalpha)*((((1/bbeta)+ ddelta - 1)/aalpha)^(aalpha/(aalpha-1))); %disutility of labor - added 9/12/13
fprintf('Calibration\n');
fprintf(' Psi = %2.6f\n', ppsi); 
fprintf('Steady State Values\n');
fprintf(' Output = %2.6f, Labor = %2.6f\n',outputSteadyState, laborSteadyState); 
fprintf(' Capital = %2.6f, Consumption = %2.6f\n',  capitalSteadyState, consumptionSteadyState); 
fprintf('\n');

% We generate the grid of capital
vGridCapital = 0.5*capitalSteadyState:0.001:1.5*capitalSteadyState; %TODO: change step to 0.00001 - modified 9/10/13
vGridLabor = 0.5*laborSteadyState:0.001:1.5*laborSteadyState; %TODO: change step to 0.00001 - modified 9/12/13

nGridLabor = length(vGridLabor); %added 9/12/13
nGridCapital = length(vGridCapital);
nGridProductivity = length(vProductivity);

%% 3. Required matrices and vectors

mOutput           = zeros(nGridCapital,nGridProductivity,nGridLabor); %modified 9/12/13
mLabor            = zeros(nGridCapital,nGridProductivity);
mValueFunction    = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

%% 4. We pre-build output for each point in the grid

for j = 1:nGridLabor
    
    mOutput(:,:,j) = (vGridCapital'.^aalpha)*vProductivity*(vGridLabor(j)^(1-aalpha));
    
end

%% 5. Main iteration

maxDifference = 10.0;
tolerance = 0.0000001;
iteration = 0;

while (maxDifference>tolerance)  
    
    expectedValueFunction = mValueFunction*mTransition';
    
    for nProductivity = 1:nGridProductivity
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nGridCapital
                        
            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                
                
                valueHighSoFar2 = -1000.0;
                
                for nLabor = 1:nGridLabor
                
                consumption = mOutput(nCapital,nProductivity,nLabor)+(1-ddelta)*vGridCapital(nCapital)-vGridCapital(nCapitalNextPeriod);
                valueProvisional = (1-bbeta)*(log(consumption)-ppsi*(vGridLabor(nLabor)^2)/2)+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);
                
                   if (valueProvisional>valueHighSoFar2)
                    valueHighSoFar2 = valueProvisional;
                    laborChoice = vGridLabor(nLabor);
                    gridLabor = nLabor;
                   else
                        break; %We break when we have achieved the max
                    end
                end
                
                
              
               
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                   break; % We break when we have achieved the max
                end    
            end
            
            mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar;
            mPolicyFunction(nCapital,nProductivity) = capitalChoice;
            mLabor(nCapital,nProductivity) = laborChoice; %Labor function choosing next period capital optimally
            
        end %end for capital grid loop
        
    end %end of productivity grid loop
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end

fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
fprintf('\n')

fprintf(' My Check = %2.8f\n', mPolicyFunction(1,1)); 
fprintf('\n')

toc

%% 6. Plotting results

figure(1)

subplot(3,1,1)
plot(vGridCapital,mValueFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Value Function')

subplot(3,1,2)
plot(vGridCapital,mPolicyFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Policy Function')

subplot(3,1,3)
plot(vGridCapital,mLabor)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Labor Function choosing Future Capital Optimally')

