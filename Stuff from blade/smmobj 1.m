function gmm = smmobj(dataMoments,Param,lambdaArray,capRArray,mArray,omega,uniformDrawsLambda,...
    parLambdaProcess,discountFactor)

% This function calculates the objective function for the SMM estimator. It
% proceeds in the following steps:
% For each value of the parameter,
%   1. Solve the model
%   2. Get optimal consumption choices for simulated health shocks
%   3. Match the moment conditions.
% 
% Inputs: 
% mChoice         :        NxT           This is the observed dataset
% Param           :        Structure     This contains parameters
% lambdaArray     :        107x39x208    Grid to solve model with
% capRArray       :        same          Grid to solve model with
% mArray          :        same          Grid to solve model with
% parLambdaProcess:        3x1           Contains parameters of health
%                                        shocks
% delta           :        scalar        Discount Rate
% 
% Output:
% gmm             :        scalar        Value of the objective function

%% solve model

capT = Param.T;

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
        = valueFunction2(lambdaArray,capRArray,...
        mArray,omega,discountFactor,uniformDrawsLambda,parLambdaProcess,vNext,...
        lambdaPreviousPeriod,t);
    
end
%% Generates optimal consumption paths given the parameters
lambdaArray1 = squeeze(lambdaArray(:,:,1));
capRArray1 = squeeze(capRArray(:,:,1));

[mSimulate,~] = dgp2(lambdaArray1,capRArray1,mOptimalAllPeriods,...
    Param);

simulatedMoments = momentcond(mSimulate,Param);
m1 = dataMoments.spendingDiff50 - simulatedMoments.spendingDiff50;
m2 = dataMoments.spendingDiff150 - simulatedMoments.spendingDiff150;
% m3 = dataMoments.spendingDiff250 - simulatedMoments.spendingDiff250;
m4 = dataMoments.meanBelow - simulatedMoments.meanBelow;
m3 = dataMoments.averageHealthShock - simulatedMoments.averageHealthShock;
m6 = dataMoments.stdHealthShock - simulatedMoments.stdHealthShock;
% m4 = dataMoments.shareOfnoHealthShock - simulatedMoments.shareOfnoHealthShock;
m = [m1;m2;m3;m4;m6];

gmm = m'*m;