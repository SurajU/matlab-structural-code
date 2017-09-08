yesDeductible = mChoice(1:Param.N/2,:);
noDeductible = mChoice(Param.N/2+1:end,:);

cumNoDeductible = cumsum(noDeductible,2);
sortedNoDeductible = sort(cumNoDeductible);
cumYesDeductible = cumsum(yesDeductible,2);

indicator = cumYesDeductible(:,1) > 375;
    newIndicator = repmat(indicator, 1, Param.T);
    
    aboveDeductibleTreatGroup = newIndicator.*cumYesDeductible;
    belowDeductibleTreatGroup = (1 - newIndicator).*cumYesDeductible;
    
    belowDeductibleTreatGroupSpending = sum(belowDeductibleTreatGroup(:,2)-belowDeductibleTreatGroup(:,1))/sum(1-indicator);
    aboveDeductibleTreatGroupSpending = sum(aboveDeductibleTreatGroup(:,2)-aboveDeductibleTreatGroup(:,1))/sum(indicator);
    
    aboveDeductibleControlGroup = sortedNoDeductible(end-sum(indicator)+1:end,1);
    
    [~,~,indNoDeductible] = intersect(aboveDeductibleControlGroup,cumNoDeductible(:,1),'rows');
    
    aboveDeductibleControlGroup = cumNoDeductible(indNoDeductible,:);
    belowDeductibleControlGroup = cumNoDeductible;
    belowDeductibleControlGroup(indNoDeductible,:) = [];
    
    aboveDeductibleIncrementalSpending(1) = aboveDeductibleTreatGroupSpending - mean(aboveDeductibleControlGroup(:,2)-aboveDeductibleControlGroup(:,1));
    belowDeductibleIncrementalSpending(1) = belowDeductibleTreatGroupSpending - mean(belowDeductibleControlGroup(:,2)-belowDeductibleControlGroup(:,1));
