function [objj,m] = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
    capRArray,mArray,omega,parLambdaProcess,discountFactor,weightingMat)

lookingAtThis = [20,2.75,1,0.99,
    omega,parLambdaProcess(2),parLambdaProcess(3),discountFactor]

% Solving the model
capT = Param.T;
Param.N = size(mChoice,1);
vAllPeriods = NaN(Param.A,Param.B,capT+1); %period T+1 will contain terminal value
vChoiceSpecificAllPeriods = NaN(Param.A,Param.B,Param.C,capT);
mOptimalAllPeriods = NaN(Param.A,Param.B,capT);

lambdaPreviousPeriod = [];

vAllPeriods(:,:,capT+1) = 0;

for t=capT:-1:1
    vNext = repmat(vAllPeriods(:,:,t+1),1,1,Param.C);
    
    [vAllPeriods(:,:,t),vChoiceSpecificAllPeriods(:,:,:,t),mOptimalAllPeriods(:,:,t)] ...
        = valueFunction2(lambdaArray,capRArray,...
        mArray,omega,discountFactor,uniformDrawsLambda,parLambdaProcess,vNext,...
        lambdaPreviousPeriod,t);
end

% for t = 2:Param.T
%     holding1 = mOptimalAllPeriods(:,:,t-1);
%     holding2 = mOptimalAllPeriods(:,:,t);
%     index = holding2 > holding1;
%     holding1(index) = holding2(index);
%     mOptimalAllPeriods(:,:,t-1) = holding1;
% end
% 
% for j = 1:size(capRArray,2)-1
%     holding2 = mOptimalAllPeriods(:,j,:);
%     holding1 = mOptimalAllPeriods(:,j+1,:);
%     index = holding2 < holding1;
%     holding1(index) = holding2(index);
%     mOptimalAllPeriods(:,j+1,:) = holding1;
% end
% 
% 
% for t = 1:Param.T-1
%     holding1 = mOptimalAllPeriods(:,:,t);
%     holding2 = mOptimalAllPeriods(:,:,Param.T);
%     index = holding2 > holding1;
%     holding1(index) = holding2(index);
%     mOptimalAllPeriods(:,:,t) = holding1;
% end


% Initialize the health grid we want to approximate our health shock
% distribution with.

healthGrid = [0:0.1:1000,1005:5:50000]';
gridDifference = [repmat(0.05,10001,1);repmat(2.5,9800,1)];
upperGrid = healthGrid + gridDifference;
lowerGrid = healthGrid - gridDifference;

% The discretized version of the health shock distribution is computed
% below.

healthHist = logncdf(upperGrid,parLambdaProcess(2),parLambdaProcess(3))...
    - logncdf(lowerGrid,parLambdaProcess(2),parLambdaProcess(3));

% Initialize arrays to be used for interpolation

newCapRArray = squeeze(capRArray(:,:,1));
newLambdaArray = squeeze(lambdaArray(:,:,1));
capR = squeeze(capRArray(1,:,1));
[newHealthGrid,newCapRGrid] = ndgrid(healthGrid,capR);

for t = 2:Param.T
    
    % Initialize the required vectors: Remaining deductible of individuals
    % and the amount they spend at time t.
    
    remainingDeductible50 = otherInfo.actualRemainingDeductible50(:,t-1);
    remainingDeductible50(any(isnan(remainingDeductible50),2)) = [];
    remainingDeductible150 = otherInfo.actualRemainingDeductible150(:,t-1);
    remainingDeductible150(any(isnan(remainingDeductible150),2)) = [];
    remainingDeductibleBelow150 = otherInfo.actualRemainingDeductibleBelow150(:,t-1);
    remainingDeductibleBelow150(any(isnan(remainingDeductibleBelow150),2)) = [];
    
    aboveDeductibleGroup = otherInfo.aboveDeductibleGroup(:,t-1);
    aboveDeductibleGroup(any(isnan(aboveDeductibleGroup),2)) = [];
    deductibleGroup50 = otherInfo.group50(:,t-1);
    deductibleGroup50(any(isnan(deductibleGroup50),2)) = [];
    deductibleGroup150 = otherInfo.group150(:,t-1);
    deductibleGroup150(any(isnan(deductibleGroup150),2)) = [];
    deductibleGroupBelow150 = otherInfo.groupBelow150(:,t-1);
    deductibleGroupBelow150(any(isnan(deductibleGroupBelow150),2)) = [];
    
    % Stacking these variables like this will make life easier later.
    stackedGroup = [deductibleGroup50;deductibleGroup150;deductibleGroupBelow150];
    
    trueRD = [remainingDeductible50;remainingDeductible150;remainingDeductibleBelow150];
    
    % Interpolate over both health grid and remaining deductible. 
    newMOptimalAllPeriods = interpn(newLambdaArray,newCapRArray,...
        mOptimalAllPeriods(:,:,t),newHealthGrid,newCapRGrid);
    % The first period requires special attention. For understanding the
    % code, this can be skipped for now. 
    
    if t == 2
        firstPeriodMOptimal = interp1(newLambdaArray(:,1),mOptimalAllPeriods(:,end-1,t-1)...
            ,newHealthGrid(:,1));
        firstPeriodPredictionError = mChoice(:,t-1) - sum(healthHist.*firstPeriodMOptimal);
        firstPeriodSecondMomentPredictionError = mChoice(:,t-1).^2 - sum(healthHist.*firstPeriodMOptimal.^2);
    end
    
    % Initialize the prediction error used for the second moment: E[y^2] -
    % E_{\lambda}[y_hat]^2 given R_it. 
    
    usedForSecondMoment = sum(newMOptimalAllPeriods.^2.*repmat(healthHist,...
        1,size(newMOptimalAllPeriods,2)));
    forSecond50 = deductibleGroup50.^2 - mean(usedForSecondMoment(1:11));
    forSecond150 = deductibleGroup150.^2 - mean(usedForSecondMoment(11:31));
    forSecondBelow150 = deductibleGroupBelow150.^2 - mean(usedForSecondMoment(31:end));
    forSecondAbove = aboveDeductibleGroup.^2 - usedForSecondMoment(1);
    SecondError(:,t-1) = [forSecondAbove;forSecond50;forSecond150;forSecondBelow150];
    
    % This interpolation is done so that mean spending is evaluated for
    % each individual at her actual remaining deductible. The order of how
    % the vector trueRD is made makes life simpler. 
    
    meanNewMOptimalAllPeriods = sum(newMOptimalAllPeriods...
        .*repmat(healthHist,1,size(newMOptimalAllPeriods,2)));
    actualDeductibleMOptimal = interp1(capR,meanNewMOptimalAllPeriods,trueRD);
    
    aboveDeductiblePredictionError = aboveDeductibleGroup - meanNewMOptimalAllPeriods(1);
    
    % Creating the matrix which contains N(T-1) prediction errors, one for each
    % person and time period.
    
    predictionError(:,t-1) = [aboveDeductiblePredictionError;stackedGroup - actualDeductibleMOptimal];
    
    % Using indicator functions to generate the dummy variables used as
    % covariates. 
    
    indicator1 = trueRD > 0 & trueRD < 50;
    indicator1 = [zeros(size(aboveDeductibleGroup,1),1);indicator1];
    indicator2 = trueRD >50 & trueRD <150;
    indicator2 = [zeros(size(aboveDeductibleGroup,1),1);indicator2];
    indicator3 = trueRD >150;
    indicator3 = [zeros(size(aboveDeductibleGroup,1),1);indicator3];
    indicator4 = [ones(size(aboveDeductibleGroup,1),1);zeros(Param.N - size...
        (aboveDeductibleGroup,1),1)];
    
    % The following lines generate the large covariate matrix Z. 
    
    columnAllocator = 4*(t-2) + 1;
    columnAllocatorEnd = 4*(t-1);
    rowAllocator = Param.N*(t-2) + 1;
    rowAllocatorEnd = Param.N*(t-1);
    Z(rowAllocator:rowAllocatorEnd,columnAllocator:columnAllocatorEnd) = [indicator1,indicator2,indicator3,indicator4];
    
end

% We need to add in the covariates for the first period. 

Z = [zeros(Param.N,44);Z];
firstPeriodInstrument = [ones(Param.N,1);zeros(Param.N*(Param.T-1),1)];

% Creating the large NT X 1 vector of prediction errors.
PredictedError = reshape(predictionError,Param.N*(Param.T-1),1);
PredictedError = [firstPeriodPredictionError;PredictedError];

% Creating the large NT X 1 vector of second moment prediction errors.
SecondPredictedError = reshape(SecondError,Param.N*(Param.T-1),1);
SecondPredictedError = [firstPeriodSecondMomentPredictionError;SecondPredictedError];

% Getting all the covariates together. 
Z1 = [firstPeriodInstrument,Z];
Z2 = [firstPeriodInstrument,Z];
Z2 = Z2(:,1);
% Matrix multiplication should yield our desired result. 
firstMomentConditions = (Z1'*PredictedError)/(Param.N);
secondMomentConditions = (Z2'*SecondPredictedError)/(Param.N);

m = [firstMomentConditions;secondMomentConditions];



objj = m'*weightingMat*m;

% if nargout >2
%     otherInfo.identifier = [(1:8000)',otherInfo.identifier];
%     otherInfo.identifier = reshape(otherInfo.identifier,Param.T*Param.N,1);
%     for i = 1:Param.N
%         j = otherInfo.identifier == i;
%         Z22(1:Param.T,:,i) = Z1(j,:);
%         Z11(1:Param.T,1,i) = Z2(j);
%         u22(1:Param.T,1,i) = PredictedError(j);
%         u32(1:Param.T,1,i) = SecondPredictedError(j);
%         mom(:,:,i) = [Z22(:,:,i)'*u22(:,1,i);Z11(:,:,i)'*u32(:,1,i)];
%         m32(:,:,i) = mom(:,:,i)*mom(:,:,i)';
%     end
% end
    
    
    
    
    
    
    
