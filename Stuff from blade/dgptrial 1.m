function [mChoice,realHealthShocks] = dgptrial(Param,omega,parLambdaProcess,discountFactor) 

% number of time periods within one year
capT = 12;

% deductible
capD = 375;



%% define grids

% health care need
lambda = [(0:5:375)'];

% remaining deductible
capR = (0:5:375)';

% health care expenditures
m = [(0:0.01:375)'];

% dimensions
Param.A = size(lambda,1);
Param.B = size(capR,1);
Param.C = size(m,1);

%   we start with the following variables
%     lambda     Ax1    a grid for lambda
%     capR       Bx1    a grid for the remaining deductible
%     m          Cx1    a grid for the chosen health care expenditures
%     omega      1x1    parameter omega

% now we construct an aligned grid of dimension AxBxC
% All this does is repeat the grid. For example, using the following
% would generate the same for lambdaArray: lambdaArray =
% repmat(lambda,1,B,C); This is done so that the value function can be
% computed for the triple lambda, R, m.
[lambdaArray,capRArray,mArray] = ndgrid(lambda,capR,m);

% Why are we evaluating at every possible combination of capRArray? Once we
% evaluate the function at mArray, wouldn't we get a capRArray by the
% function specified? Maybe some grid points are excluded later in the code
% so that the conditions hold.
%% produce random draws

% draws for health process
healthProcessDimensionsRandomness = 2;
uniformDrawsLambda = rand(numberDrawsLambda,Param.B,capT,healthProcessDimensionsRandomness);
% Why are we drawing lambdas for each value of the remaining deductible?
% Are we allowing for a general model so that we can incorporate stuff in
% later easier?

%% solve model
tic
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
toc
Param.N = 8000;
Param.T = 12;
uniformDrawsLambda = rand(Param.T,Param.N,Param.T,2);
realHealthShocks = simulateLambda2([],uniformDrawsLambda,parLambdaProcess,Param.T)';


mChoice = NaN(Param.N,Param.T);
nextR = [repmat(375,Param.N,1)];
for i = 1:Param.T
    if i == 1
        ind = realHealthShocks(:,i) < 375;
        [realHealthShocksArray,nextRArray] = ndgrid(realHealthShocks(ind(:,i),i),nextR);
        mChoice(ind(:,i),i) = diag(interpn(lambdaArray(:,:,1),capRArray(:,:,1),mOptimalAllPeriods(:,:,i),...
        realHealthShocksArray,nextRArray));
        mChoice(~ind(:,i),i) = realHealthShocks(~ind(:,i),i) + omega;
        nextR = max(nextR - mChoice(:,i),0);

    else
        ind = cumsum(mChoice(:,1:i-1),2) < 375;
        [realHealthShocksArray,nextRArray] = ndgrid(realHealthShocks(ind(:,i-1),i),nextR);
    mChoice(ind(:,i-1),i) = diag(interpn(lambdaArray(:,:,1),capRArray(:,:,1),mOptimalAllPeriods(:,:,i),...
        realHealthShocksArray,nextRArray));
    mChoice(~ind(:,i-1),i) = realHealthShocks(~ind(:,i-1),i) + omega;
    nextR = max(nextR - mChoice(:,i),0);
    end
    
end
