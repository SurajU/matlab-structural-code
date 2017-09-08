function [moment,otherInfo] = momentcond(mChoice,Param)

%% This section generates descriptive graphs for the generate dataset

% First, we construct the groups: those exposed to the deductible and those
% that aren't. This part of the code is suited for MC and not for an actual
% dataset. 
k = 1:8000;
yesDeductible = mChoice;
Param.N = size(mChoice,1);
% Now, use the B-G quantile matching strategy. 
% Generate cumulative spending for noDeductible group so that ranking based
% on that will be easier.

cumYesDeductible = cumsum(yesDeductible,2);
sortedYesDeductible = sort(cumYesDeductible);

aboveDeductibleIncrementalSpending = NaN(Param.T-1,1);
belowDeductibleIncrementalSpending = NaN(Param.T-1,1);
moment.meanAbove = NaN(Param.T-1,1);
moment.meanBelow = NaN(Param.T-1,1);
moment.spendingDiff50 = NaN(Param.T-1,1);
moment.meanSpendingDiff50 = NaN(Param.T-1,1);
moment.spendingDiff150 = NaN(Param.T-1,1);
moment.spendingDiffBelow = NaN(Param.T-1,1);
moment.spendingDiffBelow150 = NaN(Param.T-1,1);
otherInfo.actualRemainingDeductible50 = NaN(Param.N,Param.T-1);
otherInfo.actualRemainingDeductible150 = NaN(Param.N,Param.T-1);
otherInfo.actualRemainingDeductibleBelow150 = NaN(Param.N,Param.T-1);
otherInfo.number50 = NaN(Param.T-1,1);
otherInfo.number150 = NaN(Param.T-1,1);
otherInfo.numberBelow150 = NaN(Param.T-1,1);
otherInfo.aboveDeductibleGroup = NaN(Param.N,Param.T-1);
otherInfo.group50 = NaN(Param.N,Param.T-1);
otherInfo.group150 = NaN(Param.N,Param.T-1);
otherInfo.groupBelow150 = NaN(Param.N,Param.T-1);

for j = 1:Param.T-1
    indicator = cumYesDeductible(:,j) >= 375;
    otherInfo.aboveDeductible(j) = sum(indicator);
    
    % Second group within 50 euros of hitting the deductible.
    indicator2 = cumYesDeductible(:,j) >= 325 & cumYesDeductible(:,j) < 375;
    otherInfo.number50(j) = sum(indicator2);
    
    % The following line gives me the remaining deductible for the group of
    % individuals 50 euros of hitting the deductible at time j+1.
    otherInfo.actualRemainingDeductible50(1:size(cumYesDeductible(indicator2,j),1),j) = max(375 - cumYesDeductible(indicator2,j),0);
    
    % Third group [150,50] of hitting the deductible
    indicator3 = cumYesDeductible(:,j) >= 225 & cumYesDeductible(:,j) < 325;
    otherInfo.number150(j) = sum(indicator3);
    otherInfo.actualRemainingDeductible150(1:size(cumYesDeductible(indicator3,j),1),j) = max(375 - cumYesDeductible(indicator3,j),0);
    
    % Fourth group below 150 euros of hitting the deductible
    indicator4 = cumYesDeductible(:,j) < 225;
    otherInfo.numberBelow150(j) = sum(indicator4);
    otherInfo.actualRemainingDeductibleBelow150(1:size(cumYesDeductible(indicator4,j),1),j) = max(375 - cumYesDeductible(indicator4,j),0);
    otherInfo.identifier(:,j) = [k(indicator)';k(indicator2)';k(indicator3)';k(indicator4)'];
    
    % For the simple above-below deductible mean spending.
    aboveDeductibleTreatmentGroup = cumYesDeductible(indicator,:);
    belowDeductibleTreatmentGroup = cumYesDeductible(~indicator,:);
    otherInfo.aboveDeductibleGroup(1:size(aboveDeductibleTreatmentGroup,1),j) = aboveDeductibleTreatmentGroup(:,j+1) - aboveDeductibleTreatmentGroup(:,j);
    aboveDeductibleIncrementalSpending(j) = nanmean(aboveDeductibleTreatmentGroup(:,j+1) - aboveDeductibleTreatmentGroup(:,j));
    belowDeductibleIncrementalSpending(j) = nanmean(belowDeductibleTreatmentGroup(:,j+1) - belowDeductibleTreatmentGroup(:,j));
    moment.meanAbove(j) = nanmean(aboveDeductibleTreatmentGroup(:,j+1) - aboveDeductibleTreatmentGroup(:,j));
    moment.meanBelow(j) = nanmean(belowDeductibleTreatmentGroup(:,j+1) - belowDeductibleTreatmentGroup(:,j));
    moment.varAbove(j) = nanvar(aboveDeductibleTreatmentGroup(:,j+1) - aboveDeductibleTreatmentGroup(:,j));
    moment.varBelow(j) = nanvar(belowDeductibleTreatmentGroup(:,j+1) - belowDeductibleTreatmentGroup(:,j));
    
    % For spending means for those in different groups below the 
    % deductible.
    treatmentGroup50 = cumYesDeductible(indicator2,:);
    otherInfo.group50(1:size(treatmentGroup50,1),j) = treatmentGroup50(:,j+1) - treatmentGroup50(:,j);
    moment.spendingDiff50(j) = nanmean(treatmentGroup50(:,j+1)-treatmentGroup50(:,j));
    moment.var50(j) = nanvar(treatmentGroup50(:,j+1)-treatmentGroup50(:,j));
    
    % For spending of those with [150,50].
    
    treatmentGroup150 = cumYesDeductible(indicator3,:);
    otherInfo.group150(1:size(treatmentGroup150,1),j) = treatmentGroup150(:,j+1) - treatmentGroup150(:,j);
    moment.spendingDiff150(j) = nanmean(treatmentGroup150(:,j+1)-treatmentGroup150(:,j));
    moment.var150(j) = nanvar(treatmentGroup150(:,j+1)-treatmentGroup150(:,j));
    
    % For those below 150.
    treatmentGroupBelow150 = cumYesDeductible(indicator4, :);
    otherInfo.groupBelow150(1:size(treatmentGroupBelow150,1),j) = treatmentGroupBelow150(:,j+1) - treatmentGroupBelow150(:,j);
    moment.spendingDiffBelow150(j) = nanmean(treatmentGroupBelow150(:,j+1)-treatmentGroupBelow150(:,j));
    moment.varBelow150(j) = nanmean(treatmentGroupBelow150(:,j+1)-treatmentGroupBelow150(:,j));

    
end

% Some of the moments are used from below to get mean spending for first
% period. 
moment.meanSpending = mean(yesDeductible);
moment.firstSpending = mean(yesDeductible(:,1));
moment.firstQuartileSpending = prctile(yesDeductible,25);
moment.medianSpending = median(yesDeductible);
moment.thirdQuartileSpending = prctile(yesDeductible,75);
moment.meanBelow = mean(belowDeductibleIncrementalSpending);
indicatorHealth = cumYesDeductible(:,11)<375 & (cumYesDeductible(:,12)-cumYesDeductible(:,11)==0);
indicatorTotal = cumYesDeductible(:,11)<375;
moment.shareOfnoHealthShock = sum(indicatorHealth)/sum(indicatorTotal);
positiveHealthShock = cumYesDeductible;
positiveHealthShock(indicatorHealth,:) = [];
moment.averageHealthShock = nanmean(positiveHealthShock(:,12)-positiveHealthShock(:,11));
moment.stdHealthShock = std(positiveHealthShock(:,12)-positiveHealthShock(:,11));
otherInfo.numberAbove = size(mChoice,1) - (otherInfo.number50 + otherInfo.number150 + otherInfo.numberBelow150);
moment.varSpending = var(mChoice);

% Using the full vector to help identify discount factor
% Do not use the mean as moment conditions since health shocks have high
% variance == not good information. Use other quantiles instead.
% meanAbove = mean(aboveDeductibleIncrementalSpending);