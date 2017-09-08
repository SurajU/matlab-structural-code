function [objj,m] = newWay2(lambda,Param,datamoments,uniformDrawsLambda,...
    lambdaArray,capRArray,mArray,omega,parLambdaProcess,discountFactor,weightingMat,otherInfo,y)

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
        = valueFunction2(lambdaArray,capRArray,...
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


indicatorQ1 = cumsum(healthHist)<=0.25;
indicatorMedian = cumsum(healthHist)<=0.5;
indicatorQ3 = cumsum(healthHist)<=0.75;

% aboveCapR = [0;(50:5:380)']';
% aboveCapR = repmat(aboveCapR,size(gridd,1),1);
% aboveHealthGrid = repmat(gridd,1,size(aboveCapR,2));

newCapRArray = squeeze(capRArray(:,:,1));
newLambdaArray = squeeze(lambdaArray(:,:,1));
capR = squeeze(capRArray(1,:,1));
[newHealthGrid,newCapRGrid] = ndgrid(gridd,capR);

spendingPredictionError50 = NaN(Param.N,Param.T-1);
spendingPredictionErrorAbove = NaN(Param.N,Param.T-1);
spendingPredictionError150 = NaN(Param.N,Param.T-1);
spendingPredictionErrorBelow150 = NaN(Param.N,Param.T-1);
moments.meanSpending501 = NaN(Param.T-1,Param.N);
 
for t = 2:Param.T
    
    s = otherInfo.actualRemainingDeductible50(:,t-1);
    s(any(isnan(s),2)) = [];
    g = otherInfo.actualRemainingDeductible150(:,t-1);
    g(any(isnan(g),2)) = [];
    f = otherInfo.actualRemainingDeductibleBelow150(:,t-1);
    f(any(isnan(f),2)) = [];
    z = otherInfo.aboveDeductibleGroup(:,t-1);
    z(any(isnan(z),2)) = [];
    v = otherInfo.group50(:,t-1);
    v(any(isnan(v),2)) = [];
    x = otherInfo.group150(:,t-1);
    x(any(isnan(x),2)) = [];
    j = otherInfo.groupBelow150(:,t-1);
    j(any(isnan(j),2)) = [];
%     
% %     a = otherInfo.actualRemainingDeductible150(:,t-1);
% %     a(any(isnan(a),2)) = [];
%     
% %     c = [s;a];
%     
%     c = s;
%     [newHealthGrid1,newCapRGrid1] = ndgrid(gridd,c);
%     newMOptimalAllPeriods = interpn(newLambdaArray,newCapRArray,...
%         mOptimalAllPeriods(:,:,t),newHealthGrid1,newCapRGrid1);
% %     moments.meanSpending50(t-1) = mean(sum(newMOptimalAllPeriods(:,1:size(s,1)).*repmat(healthHist,1,size(s,1))));
%     moments.meanSpending501(t-1,1:size(s,1)) = sum(newMOptimalAllPeriods.*repmat(healthHist,1,size(s,1)));
% %     moments.meanSpending51(t-1) = mean(sum(newMOptimalAllPeriods.*repmat(healthHist,1,size(s,1))));
%     moments.varSpending501(t-1) = mean((sum(newMOptimalAllPeriods.^2.*repmat(healthHist,1,size(s,1)))) - moments.meanSpending501(t-1)^2);
% %     spendingPredictionError50(1:size(v,1),t-1) = v - moments.meanSpending50(t-1,1:size(s,1))';
% %     interactedSpendingPredictionError50(t-1) = sum(spendingPredictionError50(1:size(s,1),t-1).*s)/Param.N;
% %     moments.meanSpending150(t-1) = mean(sum(newMOptimalAllPeriods(:,size(s,1)+1:end).*repmat(healthHist,1,size(a,1))));
%     newMOptimalAllPeriods = interpn(newLambdaArray,newCapRArray,...
%         mOptimalAllPeriods(:,:,t),aboveHealthGrid,aboveCapR);
%     moments.meanSpendingAbove1(t-1) = sum(healthHist.*newMOptimalAllPeriods(:,1));
%     moments.varSpendingAbove1(t-1) = sum(healthHist.*(newMOptimalAllPeriods(:,1).^2)) - moments.meanSpendingAbove1(t-1)^2;
% %     spendingPredictionErrorAbove(1:size(z,1),t-1) = z - moments.meanSpendingAbove(t-1);
%     
%     moments.meanSpendingBelow1501(t-1) = mean(sum(newMOptimalAllPeriods(:,22:end).*repmat(healthHist,1,size(newMOptimalAllPeriods(:,22:end),2))));
%     moments.varSpendingBelow1501(t-1) = mean(sum(newMOptimalAllPeriods(:,22:end).^2.*repmat(healthHist,1,size(newMOptimalAllPeriods(:,22:end),2)))) - moments.meanSpendingBelow1501(t-1)^2;
% %     spendingPredictionErrorBelow150(1:size(j,1),t-1) = j - moments.meanSpendingBelow150(t-1);
% %     monthlySPEB150(t-1) = nansum(spendingPredictionErrorBelow150(:,t-1))/Param.N;
% %     interactedSpendingPredictionErrorBelow150(t-1) = sum(spendingPredictionErrorBelow150(1:size(f,1),t-1).*f)/Param.N;
%     
%     moments.meanSpending1501(t-1) = mean(sum(newMOptimalAllPeriods(:,2:22).*repmat(healthHist,1,size(newMOptimalAllPeriods(:,2:22),2))));
%     moments.varSpending1501(t-1) = mean(sum(newMOptimalAllPeriods(:,2:22).^2.*repmat(healthHist,1,size(newMOptimalAllPeriods(:,2:22),2)))) - moments.meanSpending1501(t-1)^2;
% %     spendingPredictionError150(1:size(x,1),t-1) = x - moments.meanSpending150(t-1);
% %     interactedSpendingPredictionError150(t-1) = sum(spendingPredictionError150(1:size(g,1),t-1).*g)/Param.N;


    newMOptimalAllPeriods = interpn(newLambdaArray,newCapRArray,...
         mOptimalAllPeriods(:,:,t),newHealthGrid,newCapRGrid);
     
     firstQuartileSpending = newMOptimalAllPeriods(indicatorQ1,:);
     firstQuartileSpending = firstQuartileSpending(end,:);
     moments.firstQuartileSpending(t-1) = (otherInfo.number50(t-1)'.*mean(firstQuartileSpending(1:11)) + otherInfo.number150(t-1)'.*mean(firstQuartileSpending(11:31)) + otherInfo.numberBelow150(t-1)'.*mean(firstQuartileSpending(31:end)) + otherInfo.numberAbove(t-1).*firstQuartileSpending(1))/Param.N;
     
     medianSpending = newMOptimalAllPeriods(indicatorMedian,:);
     medianSpending = medianSpending(end,:);
     moments.medianSpending(t-1) = (otherInfo.number50(t-1)'.*mean(medianSpending(1:11)) + otherInfo.number150(t-1)'.*mean(medianSpending(11:31)) + otherInfo.numberBelow150(t-1)'.*mean(medianSpending(31:end)) + otherInfo.numberAbove(t-1).*medianSpending(1))/Param.N;

     thirdQuartileSpending = newMOptimalAllPeriods(indicatorQ3,:);
     thirdQuartileSpending = thirdQuartileSpending(end,:);
     moments.thirdQuartileSpending(t-1) = (otherInfo.number50(t-1)'.*mean(thirdQuartileSpending(1:11)) + otherInfo.number150(t-1)'.*mean(thirdQuartileSpending(11:31)) + otherInfo.numberBelow150(t-1)'.*mean(thirdQuartileSpending(31:end)) + otherInfo.numberAbove(t-1).*thirdQuartileSpending(1))/Param.N;

     
     if t == 2
         moments.firstSpending = sum(newMOptimalAllPeriods(:,end).*healthHist);
         moments.varSpendingFirst = sum(newMOptimalAllPeriods(:,end).^2.*healthHist) - moments.firstSpending^2;
         moments.firstQuartileFirst = firstQuartileSpending(end);
     end
     
     % First Quartile for each remaining deductible group
     
     
     mean50 = sum(newMOptimalAllPeriods(:,1:11).*repmat(healthHist,1,size(newMOptimalAllPeriods(:,1:11),2)));
     moments.meanSpending50(t-1) = mean(mean50);
     moments.varSpending50(t-1) = mean((sum(newMOptimalAllPeriods(:,1:11).^2.*repmat(healthHist,1,size(newMOptimalAllPeriods(:,1:11),2)))) - moments.meanSpending50(t-1)^2);
     spendingPredictionError50(1:size(v,1),t-1) = v - moments.meanSpending50(t-1);
     interactedSpendingPredictionError50(t-1) = sum(spendingPredictionError50(1:size(s,1),t-1).*s)/Param.N;
     
     
     mean150 = sum(newMOptimalAllPeriods(:,11:31).*repmat(healthHist,1,size(newMOptimalAllPeriods(:,11:31),2)));
     moments.meanSpending150(t-1) = mean(mean150);
     moments.varSpending150(t-1) = mean((sum(newMOptimalAllPeriods(:,11:31).^2.*repmat(healthHist,1,size(newMOptimalAllPeriods(:,11:31),2)))) - moments.meanSpending150(t-1)^2);
     spendingPredictionError150(1:size(x,1),t-1) = x - moments.meanSpending150(t-1);
     interactedSpendingPredictionError150(t-1) = sum(spendingPredictionError150(1:size(g,1),t-1).*g)/Param.N;

     meanBelow150 = sum(newMOptimalAllPeriods(:,31:end).*repmat(healthHist,1,size(newMOptimalAllPeriods(:,31:end),2)));
     moments.meanSpendingBelow150(t-1) = mean(meanBelow150);
     moments.varSpendingBelow150(t-1) = mean((sum(newMOptimalAllPeriods(:,31:end).^2.*repmat(healthHist,1,size(newMOptimalAllPeriods(:,31:end),2)))) - moments.meanSpendingBelow150(t-1)^2);
     spendingPredictionErrorBelow150(1:size(j,1),t-1) = j - moments.meanSpendingBelow150(t-1);
     monthlySPEB150(t-1) = nansum(spendingPredictionErrorBelow150(:,t-1))/Param.N;
     interactedSpendingPredictionErrorBelow150(t-1) = sum(spendingPredictionErrorBelow150(1:size(f,1),t-1).*f)/Param.N;
     
     meanAbove = sum(newMOptimalAllPeriods(:,1).*healthHist);
     moments.meanSpendingAbove = mean(meanAbove);
     moments.varSpendingAbove(t-1) = (sum(newMOptimalAllPeriods(:,1).^2.*healthHist)) - moments.meanSpendingAbove^2;
     spendingPredictionErrorAbove(1:size(z,1),t-1) = z - moments.meanSpendingAbove;
    
end

moments.meanSpending = (otherInfo.number50'.*moments.meanSpending50 + otherInfo.number150'.*moments.meanSpending150 + otherInfo.numberBelow150'.*moments.meanSpendingBelow150 + otherInfo.numberAbove'.*moments.meanSpendingAbove)/Param.N;
moments.varSpending = (otherInfo.number50'.*moments.varSpending50 + otherInfo.number150'.*moments.varSpending150 + otherInfo.numberBelow150'.*moments.varSpendingBelow150 + otherInfo.numberAbove'.*moments.varSpendingAbove)/Param.N;
% m1 = datamoments.spendingDiff50 - moments.meanSpending50';
% m2 = datamoments.spendingDiff150 - moments.meanSpending150';
% m3 = datamoments.spendingDiffBelow150 - moments.meanSpendingBelow150';
% m4 = datamoments.meanAbove - moments.meanSpendingAbove';
m7 = datamoments.varSpending(2:end)' - moments.varSpending'; 
% m21 = [spendingPredictionError50;spendingPredictionErrorAbove;spendingPredictionErrorBelow150;spendingPredictionError150];
% m31 = nansum(m21,1);
% m31 = sum(m31);
m8 = datamoments.firstSpending - moments.firstSpending;
m9 = datamoments.varSpending(1) - moments.varSpendingFirst;
m30 = datamoments.meanSpending(2:end)' - moments.meanSpending';
m41 = datamoments.medianSpending(2:end)' - moments.medianSpending';
m51 = datamoments.thirdQuartileSpending(2:end)' - moments.thirdQuartileSpending';
% m6 = otherInfo.number50.*m1 + otherInfo.number150.*m2 + otherInfo.numberBelow150.*m3 + otherInfo.numberAbove.*m4;
% m5 = otherInfo.number50.*m1*25 + otherInfo.number150.*m2*100 + otherInfo.numberBelow150.*m3*112.5;
% m4 = datamoments.varSpending - moments.varSpending;
% m5 = m1*25 + m2*100 + m3*112.5;
% m6 = m1 + m2 + m3 + m4;


m = [m8;m9;m7;m30;m41;m51];
% if y == 1 
%     
%     m = [m1;m2;m3;m4;m5;m7;m8;m9];
%     weightingMat = eye(68);
% elseif y == 2
%     m = [m7;m1;m2;m3];
%     weightingMat = eye(44);
% elseif y == 3
%     m = [m1(5:end);m2(5:end);m3(5:end);m4(5:end);m7];
%     weightingMat = eye(39);
% elseif y == 4
%     m11 = mean(interactedSpendingPredictionError150 + interactedSpendingPredictionError50 + interactedSpendingPredictionErrorBelow150);
%     m21 = [spendingPredictionError50;spendingPredictionErrorAbove;spendingPredictionErrorBelow150;spendingPredictionError150];
%     m31 = nansum(m21,1)/Param.N;
%     m31 = sum(m31)/Param.T;
%     m = [m31;m7;m1;m2;m3];
%     weightingMat = eye(45);
% elseif y == 5
%     m11 = mean(interactedSpendingPredictionError150 + interactedSpendingPredictionError50 + interactedSpendingPredictionErrorBelow150);
%     m21 = [spendingPredictionError50;spendingPredictionErrorAbove;spendingPredictionErrorBelow150;spendingPredictionError150];
%     m31 = nansum(m21,1)/Param.N;
%     m31 = sum(m31)/Param.T;
%     m = [m11;m31;m7;m1;m2;m3];
%     weightingMat = eye(46);
% end
%     
% 
% m = [m1;m2'];

objj = m'*weightingMat*m;





