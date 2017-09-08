function [objj,m] = newWay(lambda,Param,datamoments,uniformDrawsLambda,...
    lambdaArray,capRArray,mArray,omega,parLambdaProcess,discountFactor,weightingMat)

capT = Param.T;

numberDrawsLambda = 300;
healthProcessDimensionsRandomness = 2;
% uniformDrawsLambda = rand(numberDrawsLambda,Param.B,capT,healthProcessDimensionsRandomness);

% initialize
vAllPeriods = NaN(Param.A,Param.B,capT+1); %period T+1 will contain terminal value
vChoiceSpecificAllPeriods = NaN(Param.A,Param.B,Param.C,capT);
mOptimalAllPeriods = NaN(Param.A,Param.B,capT);

% health process does not depend on previous period
% can be generalized
lambdaPreviousPeriod = [];

% terminal value
vAllPeriods(:,:,capT+1) = 0;

% backward recursion
for t=capT:-1:1
    % vNext keeps the future value function (the solution to the finite
    % horizon problem) for each value of m possible. This will come up
    % later during interpolation (remember to include).
    vNext = repmat(vAllPeriods(:,:,t+1),1,1,Param.C);
    
    [vAllPeriods(:,:,t),vChoiceSpecificAllPeriods(:,:,:,t),mOptimalAllPeriods(:,:,t)] ...
        = valueFunction(lambdaArray,capRArray,...
        mArray,omega,discountFactor,uniformDrawsLambda,parLambdaProcess,vNext,...
        lambdaPreviousPeriod,t);
    
end
% The grid below turns out to be too sparse. Need to make it finer.

% gridDifference = [repmat(5/2,76,1);repmat(10/2,62,1);repmat(5000/2,6,1)];
% lambdaGridLower = lambda - gridDifference;
% lambdaGridUpper = lambda + gridDifference;

gridd = [0:0.1:1000,1005:5:50000]';
gridDifference = [repmat(0.05,10001,1);repmat(2.5,9800,1)];
upperGrid = gridd + gridDifference;
lowerGrid = gridd - gridDifference;

healthHist = logncdf(upperGrid,parLambdaProcess(2),parLambdaProcess(3))...
    - logncdf(lowerGrid,parLambdaProcess(2),parLambdaProcess(3));

capR = squeeze(capRArray(1,:,1));
[newHealthGrid,newCapRGrid] = ndgrid(gridd,capR);
newLambdaArray = squeeze(lambdaArray(:,:,1));
newCapRArray = squeeze(capRArray(:,:,1));
 
newMOptimalAllPeriods = NaN(size(healthHist,1),size(capRArray,2),Param.T);
for t = 1:Param.T
    newMOptimalAllPeriods(:,:,t) = interpn(newLambdaArray,newCapRArray,...
        mOptimalAllPeriods(:,:,t),newHealthGrid,newCapRGrid);
end

moments = healthMoments(newMOptimalAllPeriods,Param,healthHist);
m1 = datamoments.spendingDiff50 - moments.meanSpending50;
m2 = datamoments.spendingDiff150 - moments.meanSpending150;
m3 = datamoments.spendingDiffBelow150 - moments.meanSpendingBelow150;
m4 = datamoments.meanAbove - moments.meanSpendingAbove;
% m5 = m1*25 + m2*100 + m3*112.5;
% m4 = datamoments.varSpending - moments.varSpending;

m = [m1;m2;m3;m4];

objj = m'*weightingMat*m;





