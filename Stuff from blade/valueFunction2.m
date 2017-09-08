function [v,vChoiceSpecific,mOptimal] = valueFunction2(lambdaArray,capRArray,...
    mArray,omega,discountFactor,uniformDrawsLambda,parLambdaProcess,vNext,...
    lambdaPreviousPeriod,t)
%Value function for period t
%   Detailed explanation goes here

[A,B,C] = size(lambdaArray);

% draw values of lambda to numerically integrate
% for now assume that these draws do not depend on m in the current period

% drawsLambdaNextArray = repmat(drawsLambdaNextTemp,1,1,C);

% for each draw of lambda and each value of the remaining deductible and
% health care expenditures, calculate choice specific value

uChoiceSpecific = flowUtility(lambdaArray,capRArray,mArray,omega);
capRNextArray = max(capRArray - mArray,0);


% gridDifference = [repmat(2.5,size(0:5:375,2),1);repmat(5,size(390:10:1000,2),1);repmat(2500,size(5000:5000:25000,2),1);249987500];
% % gridDifference = repmat(2.5,76,1);
% upperGrid = lambdaArray(:,1,1) + gridDifference;
% lowerGrid = lambdaArray(:,1,1) - gridDifference;

healthGrid = (0:1:390)';
gridDifference = 0.5;
upperGrid = healthGrid + gridDifference;
lowerGrid = healthGrid - gridDifference;
healthGridArray = repmat(healthGrid,1,size(capRArray,2),size(mArray,3));
mNextArray = repmat(mArray(1,:,:),length(healthGrid),1,1);
capRNextArray = repmat(capRNextArray(1,:,:),length(healthGrid),1,1);
continuationValue = interpn(lambdaArray,capRArray,mArray,vNext,...
    healthGridArray,capRNextArray,mNextArray,'linear');
% 
healthHist = logncdf(upperGrid,parLambdaProcess(2),parLambdaProcess(3))...
    - logncdf(lowerGrid,parLambdaProcess(2),parLambdaProcess(3));
healthHist = [healthHist(1:end-1);1-sum(healthHist(1:end-1))];
% mNextArray = repmat(mArray(1,:,:),length(newHealthGrid),1,1);
% capRNextArray = repmat(capRNextArrayTemp(1,:,:),length(newHealthGrid),1,1);
% newHealthGridArray = repmat(newHealthGrid,1,size(capRArray,2),size(mArray,3));
% if numberDrawsLambda<=A
%     mNextArray = mArray(1:numberDrawsLambda,:,:);
%     capRNextArray = capRNextArrayTemp(1:numberDrawsLambda,:,:);
% else
%     mNextArray = repmat(mArray(1,:,:),numberDrawsLambda,1,1);
%     capRNextArray = repmat(capRNextArrayTemp(1,:,:),numberDrawsLambda,1,1);    
% end

% continuationValue = interpn(lambdaArray,capRArray,mArray,vNext,...
%     lambdaArray,capRNextArray,mArray,'linear');


% continuationValue = interpn(lambdaArray,capRArray,mArray,vNext,...
%    drawsLambdaNextArray,capRNextArray,mNextArray,'linear');
% vChoiceSpecific = bsxfun(@plus,uChoiceSpecific,discountFactor*mean(continuationValue,1));
vChoiceSpecific = bsxfun(@plus,uChoiceSpecific,discountFactor*sum(repmat(healthHist,1,size(capRArray,2),size(mArray,3)).*continuationValue));

% determine the lowest of those values that maximizes the value function

j=isnan(vChoiceSpecific);
x = sum(sum(sum(j)));

[v,indexMaximizingM] = max(vChoiceSpecific,[],3);

m = squeeze(mArray(1,1,:));
mOptimal = m(indexMaximizingM);
allocator = size(capRArray,2);
% 
% 
% stepSize = 0.01;
% dist = size(-5:stepSize:5,2);
% newVec = 0:dist-1;
% newM = repmat(max(mOptimal(1:allocator,1:end-1)-5,0),1,1,dist);
% toBeKron = ones(allocator,allocator);
% j = kron(newVec,toBeKron);
% jR = reshape(j,allocator,allocator,dist);
% mB = newM + stepSize*jR;
% lambdaNew = repmat(lambdaArray(1:allocator,1:end-1,1),1,1,dist);
% capRNew = repmat(capRArray(1:allocator,1:end-1,1),1,1,dist);
% 
% newVChoice = interpn(lambdaArray,capRArray,mArray,vChoiceSpecific,...
%     lambdaNew,capRNew,mB,'spline');
% 
% [v(1:allocator,1:allocator),ind] = max(newVChoice,[],3);
% indr = repmat((1:allocator)',1,allocator);
% indc = repmat(1:allocator,allocator,1);
% linInd = sub2ind(size(mB),indr,indc,ind);
% mOptimal(1:allocator,1:allocator) = mB(linInd);

% tic
% if x == 0
% lambdaB = squeeze(lambdaArray(1:75,1,1));
% capRB = squeeze(capRArray(1,2:75,1));
% mB = [0:0.25:450];
% [lambdaBArray,capRBArray,mBArray] = ndgrid(lambdaB,capRB,mB);
% newVChoice = interpn(lambdaArray,capRArray,mArray,vChoiceSpecific,...
%     lambdaBArray,capRBArray,mBArray,'spline');
% 
% [vBetter,index] = max(newVChoice,[],3);
% newMOptimal = mB(index);
% mOptimal(1:75,2:75) = newMOptimal;
% v(1:75,2:75) = vBetter;
% end
% toc


exit = NaN(allocator,allocator);
for i = 1:allocator-1
    for j = 1:allocator
        [~,ind] = max(vChoiceSpecific(i,j,:),[],3);
        F = griddedInterpolant(m,squeeze(vChoiceSpecific(i,j,:)),'spline');
        [mOptimal(i,j),~,exit(i,j)] = arrayfun(@(xx)fminbnd(@(m)-F(m),max(xx-5,0),xx+5),m(ind));
        v(i,j) = F(mOptimal(i,j));
    end
end





for j = 1:length(squeeze(capRArray(1,:,1)))
    mOptimal(j:end,j) = squeeze(lambdaArray(j:end,1,1)) + omega;
    v(j:end,j) = -(j-1)*5 + (t<=12).*discountFactor^abs(t-12)*(omega/2) + (t<=11).*discountFactor^abs(t-11)*(omega/2)+...
        (t<=10).*discountFactor^abs(t-10)*(omega/2) + (t<=9).*discountFactor^abs(t-9)*(omega/2)+...
        (t<=8).*discountFactor^abs(t-8)*(omega/2) + (t<=7).*discountFactor^abs(t-7)*(omega/2)+...
        (t<=6).*discountFactor^abs(t-6)*(omega/2) + (t<=5).*discountFactor^abs(t-5)*(omega/2)+...
        (t<=4).*discountFactor^abs(t-4)*(omega/2) + (t<=3).*discountFactor^abs(t-3)*(omega/2)+...
        (t<=2).*discountFactor^abs(t-2)*(omega/2) + (t<=1).*discountFactor^abs(t-1)*(omega/2);
    
end

mOptimal = abs(mOptimal);
end

