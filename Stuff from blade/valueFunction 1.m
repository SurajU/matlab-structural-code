function [v,vChoiceSpecific,mOptimal] = valueFunction(lambdaArray,capRArray,...
    mArray,omega,discountFactor,uniformDrawsLambda,parLambdaProcess,vNext,...
    lambdaPreviousPeriod,t)
%Value function for period t
%   Detailed explanation goes here

[A,B,C] = size(lambdaArray);
numberDrawsLambda = size(uniformDrawsLambda,1);

% draw values of lambda to numerically integrate
% for now assume that these draws do not depend on m in the current period
drawsLambdaNextTemp = simulateLambda(lambdaPreviousPeriod,uniformDrawsLambda,...
    parLambdaProcess,t);
drawsLambdaNextArray = repmat(drawsLambdaNextTemp,1,1,C);

% for each draw of lambda and each value of the remaining deductible and
% health care expenditures, calculate choice specific value

uChoiceSpecific = flowUtility(lambdaArray,capRArray,mArray,omega);
capRNextArrayTemp = max(capRArray - mArray,0);

if numberDrawsLambda<=A
    mNextArray = mArray(1:numberDrawsLambda,:,:);
    capRNextArray = capRNextArrayTemp(1:numberDrawsLambda,:,:);
else
    mNextArray = repmat(mArray(1,:,:),numberDrawsLambda,1,1);
    capRNextArray = repmat(capRNextArrayTemp(1,:,:),numberDrawsLambda,1,1);    
end

continuationValue = interpn(lambdaArray,capRArray,mArray,vNext,...
    drawsLambdaNextArray,capRNextArray,mNextArray,'linear');
vChoiceSpecific = bsxfun(@plus,uChoiceSpecific,discountFactor*mean(continuationValue,1));

% determine the lowest of those values that maximizes the value function
[v,indexMaximizingM] = max(vChoiceSpecific,[],3);

m = squeeze(mArray(1,1,:));
mOptimal = m(indexMaximizingM);

end

