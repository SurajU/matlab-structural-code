function moments = healthMoments(mOptimalAllPeriods,Param,healthHist)
% Here, I attempt to compute the histograms for a given lognormal health
% process. First, I get the corresponding grid and then compute the
% histogram based on those grids. 


% creating the first spending moments, given how far the consumers are from
% hitting the deductible. 

% first off, we create the mean spending moment. This should not be used in
% the final version since health processes are extremely noisy. 

% mean spending for those 50 euros from hitting the deductible.
moments.medianSpending50 = NaN(Param.T-1,1);
moments.medianSpending150 = NaN(Param.T-1,1);
moments.varSpending = NaN(Param.T-1,1);

for t = 2:Param.T
    % this is the mean.
    % first off, choose the relevant optimal spending values. In our case,
    % we want those within 50 euros of exceeding the deductible. This
    % corresponds to columns 1:11 in mOptimalAllPeriods
    % this is the median for 50 euros near the deductible.
    indicator = cumsum(healthHist)<=0.5;
    % need to add a condition where there's nothing here. 
    relevantGroupMedian = mOptimalAllPeriods(indicator,1:11,t);
    if size(relevantGroupMedian,1)<1
        relevantGroupMedian = mOptimalAllPeriods(1,1:11,t);
    end
    moments.medianSpending50(t-1) = mean(relevantGroupMedian(end,:));
    
    % median for 150 euros below the deductible. WE ARE ASSUMING NO
    % DIFFERENCE IN HEALTH SHOCK DISTRIBUTION FOR PEOPLE AT DIFFERENT
    % SEGMENTS OF THE DEDUCTIBLE!!
    relevantGroupMedian150 = mOptimalAllPeriods(indicator,11:31,t);
    if size(relevantGroupMedian150,1)<1
        relevantGroupMedian150 = mOptimalAllPeriods(1,11:31,t);
    end
    moments.medianSpending150(t-1) = mean(relevantGroupMedian150(end,:));

%     moments.meanSpendingAll(t-1) = sum(mOptimalAllPeriods(:,1,t).*healthHist);
% the formula for variance is wrong!
    moments.varSpending(t-1) = sum(mean(mOptimalAllPeriods(:,:,t),2).^2.*healthHist)-(sum(mean(mOptimalAllPeriods(:,:,t),2).*healthHist)^2);
end

% Remember, might need to generalize these moments and include more time
% variation in mean and variance of spending for more complicated health
% processes. 
% the moment condition for meanSpendingAll seems to be pretty lousy. Need
% to modify this. 
moments.meanSpending50 = squeeze(mean((sum(mOptimalAllPeriods(:,1:11,:).*repmat(healthHist,1,size(mOptimalAllPeriods(:,1:11,:),2),Param.T))),2));
moments.meanSpending50 = moments.meanSpending50(2:end);

moments.meanSpending150 = squeeze(mean((sum(mOptimalAllPeriods(:,11:31,:).*repmat(healthHist,1,size(mOptimalAllPeriods(:,11:31,:),2),Param.T))),2));
moments.meanSpending150 = moments.meanSpending150(2:end);

moments.meanSpendingAbove = squeeze(sum(mOptimalAllPeriods(:,1,:).*repmat(healthHist,1,1,Param.T)));
moments.meanSpendingAbove = moments.meanSpendingAbove(2:end);

moments.meanSpendingBelow150 = squeeze(mean(sum(mOptimalAllPeriods(:,31:end,:).*repmat(healthHist,1,size(mOptimalAllPeriods(:,31:end,:),2),Param.T)),2));
moments.meanSpendingBelow150 = moments.meanSpendingBelow150(2:end);

moments.varSpending = mean(moments.varSpending);

% check the above moments with the corresponding "actual" values. 




    