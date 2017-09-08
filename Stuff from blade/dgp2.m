function [mChoice,realHealthShocks] = dgp2(lambdaArray,capRArray,...
    mOptimalAllPeriods,Param)
%   This function generates the data for random draws of health shocks from
%   the function simulateLambda. The primary data generated is the choice
%   of consumption in each period, given the health shock faced in the
%   period. 
%   
%   Input:
%   lambdaArray       :    AxB        Grid for possible values of lambda
%   capRArray         :    AxB        Grid for possible values of remaining
%                                       deductible
%   mOptimalAllPeriods:    AxBxT      Optimal Choice given lambda and capR
%   Param             :    Structure  Contains all parameters required for
%                                       running this code.
%     
%   Output:
%   mChoice           :    NxT        Gives the choice of consumers from
%                                       the finite horizon dynamic problem.
%   realHealthShocks  :    NxT        This is a matrix of health shocks
%     
%   FUNCTIONS CALLED:
%     
%   simulateLambda.m 
%     

%% 

% First, I need to draw the health shocks for all the periods, for all
% individuals. 

uniformDrawsLambda = rand(Param.T,Param.N,Param.T,Param.healthProcessDimensionsRandomness);
realHealthShocks = simulateLambda2([],uniformDrawsLambda,Param.parLambdaProcess,Param.t);


%% 

% Now, I need to calculate the optimal m given the first period health
% shock and the remaining deductible. In the first period, the remaining
% deductible is always 375 (except for those under 18).

mChoice = NaN(Param.N,Param.T);
nextR = [repmat(375,Param.N,1)];

%% 


for i = 1:Param.T
    [realHealthShocksArray,nextRArray] = ndgrid(realHealthShocks(i,:),nextR);
    mChoice(:,i) = diag(interpn(lambdaArray,capRArray,mOptimalAllPeriods(:,:,i),...
        realHealthShocksArray,nextRArray));
    nextR = max(nextR - mChoice(:,i),0);
end







    