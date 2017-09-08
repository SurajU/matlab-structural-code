function [objj,m] = newWay2(mChoice,lambda,Param,datamoments,uniformDrawsLambda,...
    lambdaArray,capRArray,mArray,omega,parLambdaProcess,discountFactor,weightingMat,otherInfo,y)

capT = Param.T;

yesDeductible = mChoice;
cumYesDeductible = cumsum(yesDeductible,2);

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

% Creating the mid-point of the histogram for the health process.

gridd = [0:0.1:1000,1005:5:50000]';
gridDifference = [repmat(0.05,10001,1);repmat(2.5,9800,1)];
upperGrid = gridd + gridDifference;
lowerGrid = gridd - gridDifference;

healthHist = logncdf(upperGrid,parLambdaProcess(2),parLambdaProcess(3))...
    - logncdf(lowerGrid,parLambdaProcess(2),parLambdaProcess(3));

% This creates the indicators used later to get quantile moments. Has to be
% generalized if the health process is not iid.

indicatorQ1 = cumsum(healthHist)<=0.25;
indicatorMedian = cumsum(healthHist)<=0.5;
indicatorQ3 = cumsum(healthHist)<=0.75;

% Initialize the arrays to be used in interpolation. 

newCapRArray = squeeze(capRArray(:,:,1));
newLambdaArray = squeeze(lambdaArray(:,:,1));
capR = squeeze(capRArray(1,:,1));
[newHealthGrid,newCapRGrid] = ndgrid(gridd,capR);

% Initialize some of the vectors. Most of them are not included here but
% will be later, in order to speed up the code. 

spendingPredictionError50 = NaN(Param.N,Param.T-1);
spendingPredictionErrorAbove = NaN(Param.N,Param.T-1);
spendingPredictionError150 = NaN(Param.N,Param.T-1);
spendingPredictionErrorBelow150 = NaN(Param.N,Param.T-1);
moments.meanSpending501 = NaN(Param.T-1,Param.N);
 
for t = 2:Param.T
    
    % The following just stores the relevant variables into a single
    % alphabet. This will be used primarily for interpolating at true
    % remaining deductibile values and to get interacted predicted errors.
    
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
    
    % bs is the vector of true remaining deductible for each time period.
    % This vector changes every time period since there are different
    % people in each spending region at different time periods.
    
    bs = [s;g;f]';
%{    
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
%}

    % The following lines interpolate optimal consumption for the fine grid
    % of health shocks we have used. 
    
    newMOptimalAllPeriods = interpn(newLambdaArray,newCapRArray,...
         mOptimalAllPeriods(:,:,t),newHealthGrid,newCapRGrid);
     
     % When we're in the first time period, everyone has the same remaining
     % deductible, 375. Thus, this has to be treated differently than the
     % others. The following lines of code just compute the model
     % predictions for mean spending and variance spending for the first
     % period.
     
     if t == 2
         moments.firstSpending = sum(newMOptimalAllPeriods(:,end).*healthHist);
         moments.varSpendingFirst = sum(newMOptimalAllPeriods(:,end)...
             .^2.*healthHist) - moments.firstSpending^2;
     end
     
     % In the following lines, we interpolate mean optimal consumption on
     % the true remaining deductible amounts individual have, which is
     % obtained from the data.
     
     meanNewMOptimalAllPeriods = sum(newMOptimalAllPeriods...
         .*repmat(healthHist,1,size(newMOptimalAllPeriods,2)));
     actualDeductibleMOptimal = interp1(capR,meanNewMOptimalAllPeriods,bs);
     
     % The mean moment conditions are computed immediately from the results
     % of interpolation done in the previous lines. Basically, the mean
     % spending of each individual, in a spending region, is computed and 
     % the mean is taken over all the individual mean spending of that 
     % spending region.
     % The vectors s,g,f represents the individuals in 3 different spending
     % regions respectively.
     
     moments.meanSpending50(t-1) = mean(actualDeductibleMOptimal(1:size(s,1)));
     moments.meanSpending150(t-1) = mean...
         (actualDeductibleMOptimal(size(s,1)+1:size(s,1)+size(g,1)));
     moments.meanSpendingBelow150(t-1) = mean...
         (actualDeductibleMOptimal(size(s,1)+size(g,1)+1:size(f,1)+size(s,1)+size(g,1)));
     
     % There is a slight difference in computing the variance. Here,
     % instead of using the true remaining deductible, I just took the mean
     % over all relevant remaining deductible amounts and computed 
     % conditional variances first. The formulation for this variance is
     % just the variance of a discrete r.v.
     
     moments.varSpending50(t-1) = mean((sum(newMOptimalAllPeriods(:,1:11).^2.*...
         repmat(healthHist,1,size(newMOptimalAllPeriods(:,1:11),2))))...
         - moments.meanSpending50(t-1)^2);
     
     % The following is the prediction error of the people in the relevant
     % spending region. The interacted spending prediction error is just
     % the prediction error interacted with remaining deductible. 
     % v is the vector of true period t spending, for individuals with 50
     % euros left in the deductible. 
     
     spendingPredictionError50(1:size(v,1),t-1) = v - moments.meanSpending50(t-1);
     interactedSpendingPredictionError50(t-1) = sum(spendingPredictionError50(1:size(s,1),t-1).*s)/Param.N;
     
     % Same as above, just for the spending region of 150 euros of hitting
     % the deductible. 
     
     moments.varSpending150(t-1) = mean((sum(newMOptimalAllPeriods(:,11:31).^2.*repmat(healthHist,1,size(newMOptimalAllPeriods(:,11:31),2)))) - moments.meanSpending150(t-1)^2);
     spendingPredictionError150(1:size(x,1),t-1) = x - moments.meanSpending150(t-1);
     interactedSpendingPredictionError150(t-1) = sum(spendingPredictionError150(1:size(g,1),t-1).*g)/Param.N;
     
     % Same as above, just below 150 of hitting the deductible. 
     
     moments.varSpendingBelow150(t-1) = mean((sum(newMOptimalAllPeriods(:,31:end).^2.*repmat(healthHist,1,size(newMOptimalAllPeriods(:,31:end),2)))) - moments.meanSpendingBelow150(t-1)^2);
     spendingPredictionErrorBelow150(1:size(j,1),t-1) = j - moments.meanSpendingBelow150(t-1);
     interactedSpendingPredictionErrorBelow150(t-1) = sum(spendingPredictionErrorBelow150(1:size(f,1),t-1).*f)/Param.N;
     
     % Since we know the remaining deductible for those above, the
     % computations become much easier. The following lines obtain the
     % relevant moment conditions for those above the deductible.
     
     meanAbove = sum(newMOptimalAllPeriods(:,1).*healthHist);
     moments.meanSpendingAbove = mean(meanAbove);
     moments.varSpendingAbove(t-1) = (sum(newMOptimalAllPeriods(:,1).^2.*healthHist)) - moments.meanSpendingAbove^2;
     spendingPredictionErrorAbove(1:size(z,1),t-1) = z - moments.meanSpendingAbove;
     
     % Now to write quantile moment conditions.First, to define the groups
     indicator50 = cumYesDeductible(:,t-1) >= 325 & cumYesDeductible(:,t-1) < 375;
     indicatorAbove = cumYesDeductible(:,t-1) >= 375;
     indicator150 = cumYesDeductible(:,t-1) >= 225 & cumYesDeductible(:,t-1) < 325;
     indicatorBelow150 = cumYesDeductible(:,t-1) < 225;
     
     % This picks out the relevant spending of those with the following
     % remaining deductible amounts.
     
     those50 = cumYesDeductible(indicator50,t)-cumYesDeductible(indicator50,t-1);
     thoseAbove = cumYesDeductible(indicatorAbove,t) - cumYesDeductible(indicatorAbove,t-1);
     those150 = cumYesDeductible(indicator150,t) - cumYesDeductible(indicator150,t-1);
     thoseBelow150 = cumYesDeductible(indicatorBelow150,t) - cumYesDeductible(indicatorBelow150,t-1);
     
     % The first following line just picks out the model prediction for the
     % 1st quartile of the health shock. It must be mentioned that
     % monotonicity in the mapping of health shocks to consumption is
     % assumed. 
     % Recall that indicatorQ1 just picks the exact cell which is the first
     % quartile health shock specified in healthHist.
     
     firstQuartileSpending = newMOptimalAllPeriods(indicatorQ1,:);
     if  size(firstQuartileSpending,1) == 0
         firstQuartileSpending = newMOptimalAllPeriods(1,:);
     end
     firstQuartileSpending = firstQuartileSpending(end,:);
     
     % The following lines create the moment conditions. First, we obtain
     % the model predictions for first quartile spending by just taking the
     % mean over all relevant remaining deductible amounts for the optimal
     % consumption, which was mapped onto the first quartile of the health 
     % shock process. Then, the spending of individuals in the relevant
     % region is compared to it. We would require that 1/4th of the
     % individuals in the relevant spending region to have spending less
     % than or equal to the model predicted first quartile.
     
     moments.firstQuartileSpending50(t-1) = mean(firstQuartileSpending(1:11));
     indi = those50<= moments.firstQuartileSpending50(t-1);
     moments.first50(t-1) = mean(indi);
     
     moments.firstQuartileSpending150(t-1) =  mean(firstQuartileSpending(11:31));
     indi = those150<= moments.firstQuartileSpending150(t-1);
     moments.first150(t-1) = mean(indi);
     
     moments.firstQuartileSpendingBelow150(t-1) = mean(firstQuartileSpending(31:end));
     indi = thoseBelow150 <= moments.firstQuartileSpendingBelow150(t-1);
     moments.firstBelow150(t-1) = mean(indi);
     
     moments.firstQuartileSpendingAbove(t-1) = firstQuartileSpending(1);
     indi = thoseAbove<= moments.firstQuartileSpendingAbove(t-1);
     moments.firstAbove(t-1) = mean(indi);
     
     % Does the same as above, but just for the median.
     
     medianSpending = newMOptimalAllPeriods(indicatorMedian,:);
     if size(medianSpending,1) == 0
         medianSpending = newMOptimalAllPeriods(1,:);
     end
     medianSpending = medianSpending(end,:);
     
     moments.medianSpending50(t-1) = mean(medianSpending(1:11));
     indi = those50 <= moments.medianSpending50(t-1);
     moments.med50(t-1) = mean(indi);
     
     moments.medianSpending150(t-1) = mean(medianSpending(11:31));
     indi = those150 <= moments.medianSpending150(t-1);
     moments.med150(t-1) = mean(indi);
     
     moments.medianSpendingBelow150(t-1) = mean(medianSpending(31:end));
     indi = thoseBelow150 <= moments.medianSpendingBelow150(t-1);
     moments.medBelow150(t-1) = mean(indi);
     
     moments.medianSpendingAbove(t-1) = medianSpending(1);
     indi = thoseAbove <= moments.medianSpendingAbove(t-1);
     moments.medAbove(t-1) = mean(indi);
     
     % Does the same as above, just for the third quartile. 
     
     thirdQuartileSpending = newMOptimalAllPeriods(indicatorQ3,:);
     if size(thirdQuartileSpending,1) == 0
         thirdQuartileSpending = newMOptimalAllPeriods(1,:);
     end
     thirdQuartileSpending = thirdQuartileSpending(end,:);
     
     moments.thirdQuartileSpending50(t-1) = mean(thirdQuartileSpending(1:11));
     indi = those50 <= moments.thirdQuartileSpending50(t-1);
     moments.third50(t-1) = mean(indi);
     
     moments.thirdQuartileSpending150(t-1) = mean(thirdQuartileSpending(11:31));
     indi = those150 <= moments.thirdQuartileSpending150(t-1);
     moments.third150(t-1) = mean(indi);
     
     moments.thirdQuartileSpendingBelow150(t-1) = mean(thirdQuartileSpending(31:end));
     indi = thoseBelow150 <= moments.thirdQuartileSpendingBelow150(t-1);
     moments.thirdBelow150(t-1) = mean(indi);
     
     moments.thirdQuartileSpendingAbove(t-1) = thirdQuartileSpending(1);
     indi = thoseAbove <= moments.thirdQuartileSpendingAbove(t-1);
     moments.thirdAbove(t-1) = mean(indi);
     
     % When in the first period, again, everyone has the same remaining
     % deductible, so quartile moment conditions become easy to compute.
     
     if t == 2 
         firstQ1 = mChoice(:,t-1) <= firstQuartileSpending(end);
         moments.FQ1 = mean(firstQ1);
         
         firstmed = mChoice(:,t-1) <= medianSpending(end);
         moments.FM = mean(firstmed);
         
         firstQ3 = mChoice(:,t-1) <= thirdQuartileSpending(end);
         moments.FQ3 = mean(firstQ3);
     end
         
     
end

Param.g = Param.N*Param.T;
moments.varSpending = (otherInfo.number50'.*moments.varSpending50 + otherInfo.number150'.*moments.varSpending150 + otherInfo.numberBelow150'.*moments.varSpendingBelow150 + otherInfo.numberAbove'.*moments.varSpendingAbove)/Param.N;
secondTermForVariance = (otherInfo.number50'.*moments.meanSpending50.^2.*(Param.N-otherInfo.number50') + otherInfo.number150'.*moments.meanSpending150.^2.*(Param.N-otherInfo.number150') + otherInfo.numberBelow150'.*moments.meanSpendingBelow150.^2.*(Param.N-otherInfo.numberBelow150') + otherInfo.numberAbove'.*moments.meanSpendingAbove.^2.*(Param.N-otherInfo.numberAbove'))/Param.N^2;
thirdTermForVariance = (otherInfo.number50'.*moments.meanSpending50.*otherInfo.number150'.*moments.meanSpending150 + otherInfo.number50'.*moments.meanSpending50.*otherInfo.numberBelow150'.*moments.meanSpendingBelow150 + otherInfo.number50'.*moments.meanSpending50.*otherInfo.numberAbove'.*moments.meanSpendingAbove + otherInfo.number150'.*moments.meanSpending150.*otherInfo.numberBelow150'.*moments.meanSpendingBelow150 + otherInfo.number150'.*moments.meanSpending150.*otherInfo.numberAbove'.*moments.meanSpendingAbove + otherInfo.numberBelow150'.*moments.meanSpendingBelow150.*otherInfo.numberAbove'.*moments.meanSpendingAbove)/Param.N^2;
moments.varSpending = moments.varSpending + secondTermForVariance - 2*thirdTermForVariance;

% The following just lists the moment conditions and chooses them. 

m1 = datamoments.spendingDiff50 - moments.meanSpending50';
m2 = datamoments.spendingDiff150 - moments.meanSpending150';
m3 = datamoments.spendingDiffBelow150 - moments.meanSpendingBelow150';
m4 = datamoments.meanAbove - moments.meanSpendingAbove';
m7 = datamoments.varSpending(2:end)' - moments.varSpending'; 
m21 = [spendingPredictionError50;spendingPredictionErrorAbove;spendingPredictionErrorBelow150;spendingPredictionError150];
m31 = nansum(m21,1);
m31 = sum(m31);
m8 = datamoments.firstSpending - moments.firstSpending;
m9 = datamoments.varSpending(1) - moments.varSpendingFirst;
m11 = sum(interactedSpendingPredictionError150 + interactedSpendingPredictionError50 + interactedSpendingPredictionErrorBelow150);
v1 = datamoments.varAbove' - moments.varSpendingAbove';
v2 = datamoments.var50' - moments.varSpending50';
v3 = datamoments.var150' - moments.varSpending150';
v4 = datamoments.varBelow150' - moments.varSpendingBelow150';
% m6 = otherInfo.number50.*m1 + otherInfo.number150.*m2 + otherInfo.numberBelow150.*m3 + otherInfo.numberAbove.*m4;
% m5 = otherInfo.number50.*m1*25 + otherInfo.number150.*m2*100 + otherInfo.numberBelow150.*m3*112.5;
% m4 = datamoments.varSpending - moments.varSpending;
m5 = otherInfo.number50.*m1*25 + otherInfo.number150.*m2*100 + otherInfo.numberBelow150.*m3*112.5;
m6 = m1 + m2 + m3 + m4;
% if y == 1 
%     
%     m = [otherInfo.number50.*m1/Param.N;otherInfo.number150.*m2/Param.N;otherInfo.numberBelow150.*m3/Param.N;otherInfo.numberAbove.*m4/Param.N;m5;m7;m8;m9];
%     weightingMat = eye(68);
% elseif y == 2
%     m = [m7;m1;m2;m3];
%     weightingMat = eye(44);
% elseif y == 3
%     m = [m1(5:end);m2(5:end);m3(5:end);m4(5:end);m7];
%     weightingMat = eye(39);
% elseif y == 4
%     m = [m31;m7;m1;m2;m3];
%     weightingMat = eye(45);
% elseif y == 5
%      m = [m11;m31;m7;otherInfo.number50.*m1/Param.g;otherInfo.number150.*m2/Param.g;otherInfo.numberBelow150.*m3/Param.g;otherInfo.numberAbove.*m4/Param.g];
%     weightingMat = eye(57);

% Quantile moments
% First Quartile
q0 = moments.FQ1 - 0.25;
q1 = moments.first50 - 0.25;
q2 = moments.first150 - 0.25;
q3 = moments.firstBelow150 - 0.25;
q4 = moments.firstAbove - 0.25;

% Second Quartile
q01 = moments.FM - 0.5;
q5 = moments.med50 - 0.5;
q6 = moments.med150 - 0.5;
q7 = moments.medBelow150 - 0.5;
q8 = moments.medAbove - 0.5;

% Third Quartile
q02 = moments.FQ3 - 0.75;
q9 = moments.third50 - 0.75;
q10 = moments.third150 - 0.75;
q11 = moments.thirdBelow150 - 0.75;
q12 = moments.thirdAbove - 0.75;

if y == 1
    m = [m11;m31;m8];
    weightingMat = eye(3);
elseif y == 2
    m = [m11;m31;m8;m1;m2;m3;m4];
    weightingMat = eye(47);
elseif y == 3
    m = [m11;otherInfo.number50.*m1/Param.g;otherInfo.number150.*m2/Param.g;otherInfo.numberBelow150.*m3/Param.g;otherInfo.numberAbove.*m4/Param.g;m8;m7;m9];
    weightingMat = eye(58);
elseif y == 4
    m = [m1;m2;m3;m4;m8];
    weightingMat = eye(45);
elseif y == 5
    m = [otherInfo.number50.*m1/Param.g;otherInfo.number150.*m2/Param.g;otherInfo.numberBelow150.*m3/Param.g;otherInfo.numberAbove.*m4/Param.g;m8];
    weightingMat = eye(45);
elseif y == 6
    m = [m11;otherInfo.number50.*m1/Param.g;otherInfo.number150.*m2/Param.g;otherInfo.numberBelow150.*m3/Param.g;otherInfo.numberAbove.*m4/Param.g;m8];
    weightingMat = eye(46);
elseif y == 7
    m = [m11;m1;m2;m3;m4;m8];
    weightingMat = eye(46);
elseif y == 8
    m = [q0;q01;q02;q1';q2';q3';q4';q5';q6';q7';q8';q9';q10';q11';q12';m9;m1;m2;m3;m4;v1;v2;v3;v4;m8];
    weightingMat = eye(225);
elseif y == 9
    m = [m8;m1;m2;m3;m4;v1;v2;v3;v4;m9];
    weightingMat = eye(90);
elseif y == 10
    m = [m8;m1;m2;m3;m4;m7;m9];
    weightingMat = eye(57);
elseif y == 11
    m = [m11;m8;m1;m2;m3;m4;m7;m9];
    weightingMat = eye(58);
elseif y == 12
    m = [m11;otherInfo.number50.*m1/Param.g;otherInfo.number150.*m2/Param.g;otherInfo.numberBelow150.*m3/Param.g;otherInfo.numberAbove.*m4/Param.g;m8];
    weightingMat = eye(4);
end
    

% 
% m = [m1;m2'];

objj = m'*weightingMat*m;





