% structural model, main file


% %% clean up and set seed, set number random draws
clear;

numberDrawsLambda = 300;
%% define parameters, some of them will be estimated later
rng('default');
rng(1232);
% number of time periods within one year
capT = 12;

% deductible
capD = 375;

% moral hazard parameter
omega = 20;

% discount factor
discountFactor = 0.99;

% health process
% first element: parameter for probability that it is zero, modeled
%                as 1/(1+exp(-parLambdaProcess(1))
% second and third element: mean and std of underlying normal
%                           distribution if positive

parLambdaProcess = [0.1;2.75;1];


%% define grids for data generation

% health care need
lambda = [(0:5:375)';(390:10:1000)';(5000:5000:25000)';500000000];

% remaining deductible
capR = (0:5:380)';

% health care expenditures
m = [(0:5:380)';(390:10:1000)';(5000:5000:25000)';50000;1000000000];

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

% Why are we drawing lambdas for each value of the remaining deductible?
% Are we allowing for a general model so that we can incorporate stuff in
% later easier?
healthProcessDimensionsRandomness = 2;
uniformDrawsLambda = rand(numberDrawsLambda,Param.B,capT,healthProcessDimensionsRandomness);
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

% Just to make solutions consistent with theory
% for t = 2:capT
%     holding1 = mOptimalAllPeriods(:,:,t-1);
%     holding2 = mOptimalAllPeriods(:,:,t);
%     index = holding2 > holding1;
%     holding1(index) = holding2(index);
%     mOptimalAllPeriods(:,:,t-1) = holding1;
% end
% % 
% for j = 1:size(capRArray,2)-1
%     holding2 = mOptimalAllPeriods(:,j,:);
%     holding1 = mOptimalAllPeriods(:,j+1,:);
%     index = holding2 < holding1;
%     holding1(index) = holding2(index);
%     mOptimalAllPeriods(:,j+1,:) = holding1;
% end
% 
% for t = 1:capT-1
%     holding1 = mOptimalAllPeriods(:,:,t);
%     holding2 = mOptimalAllPeriods(:,:,capT);
%     index = holding2 > holding1;
%     holding1(index) = holding2(index);
%     mOptimalAllPeriods(:,:,t) = holding1;
% end

    Param.N = 8000;
    Param.T = capT;
    Param.parLambdaProcess = parLambdaProcess;
    Param.t = 1;
    Param.numberDrawsLambda = 10;
    Param.healthProcessDimensionsRandomness = 2;
    
%% This section generates a Monte Carlo dataset
parfor i = 1:20
    lambdaArray1 = squeeze(lambdaArray(:,:,1));
    capRArray1 = squeeze(capRArray(:,:,1));
    
    [mChoice,realHealthShocks] = dgp2(lambdaArray1,capRArray1,mOptimalAllPeriods,...
        Param);

    [dataMoments,otherInfo] = momentcond(mChoice,Param);
        
    % tic
    % gentrial = smmobj(dataMoments,Param,lambdaArray,capRArray,mArray,omega,uniformDrawsLambda,...
    %     parLambdaProcess,discountFactor);
    % toc
    
    %     startvalues = [15,1.5,0.65,0.7*discountFactor];
    %     lb = [-Inf,-Inf,0,0];
    %     ub = [Inf,Inf,Inf,1];
    %     objFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
    %         capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],estimates(4));
    %     options = optimset('Display','iter','TolFun',1e-10,'TolX',1e-6,...
    %         'Algorithm','interior-point','AlwaysHonorConstraints','bounds');
    %
    %     estimatett(:,i) = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
    
    
    % discountFactorGrid2 = [0.1:0.01:1];
    % %     parfor k = 1:length(discountFactorGrid2)
    % %         evaluation(k) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
    % %             capRArray,mArray,estimatted(1,7,3),[0.1;estimatted(2,7,3);estimatted(3,7,3)],discountFactorGrid2(k));
    % %     end
    %
    % parfor k = 1:length(discountFactorGrid2)
    %     evaluation3(k) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
    %         capRArray,mArray,omega,parLambdaProcess,discountFactorGrid2(k));
    % end
    %
    % plot(discountFactorGrid2,evaluation3)
    %
    %
%     weightingMat = eye(90);
weightingMat = eye(46);
    tic
    eval = properMoments(mChoice,Param,otherInfo,lambdaArray,...
            capRArray,mArray,omega,parLambdaProcess,0.5,weightingMat);
    toc
        %     startvalues = [15,1.5,0.65];
    %     lb = [-Inf,-Inf,0];
    %     ub = [Inf,Inf,Inf];
    %     options = optimset('Display','iter','TolFun',1e-10,'TolX',1e-6,...
    %         'Algorithm','interior-point','AlwaysHonorConstraints','bounds');
    %     discountFactorGrid2 = [0.1:0.05:0.95];
    %     parfor j = 1:length(discountFactorGrid2)
    %         objFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
    %             capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],discountFactorGrid2(j));
    %         [fminconEstimatesTrial(:,j),objFunValTrial(:,j)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
    %     end
    %
    weightingMat = eye(46);
    objectiveFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
        capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],estimates(4),weightingMat);
    startvalues = [15,1.5,0.65,0.7];
    lb = [-Inf,-Inf,0,0];
    ub = [Inf,Inf,Inf,1];
    options = optimset('Display','iter','UseParallel',true);
    tic
    estt(:,i) = fmincon(objectiveFun,startvalues,[],[],[],[],lb,ub,[],options);
    toc
end
%     estt2 = knitromatlab(objectiveFun,startvalues,[],[],[],[],lb,ub,[],[],options);
%     toc
    
%     [~,m] = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%              capRArray,mArray,estt(1),[0.1,estt(2),estt(3)],estt(4),weightingMat);
% %     j = inv(m'*m);
%     g = j(1,1);
%     weightingMat = [j,zeros(length(j),1);zeros(1,length(j)),g];
    
%     weightinMat = inv(m*m');     
%     startvalues = estt;
%     objectiveFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%         capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],estimates(4),weightingMat);
%     twoStepEst2 = fmincon(objectiveFun,startvalues,[],[],[],[],lb,ub,[],options);
%     
% tic
%     func = @(parameters)momentConditions(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%         capRArray,mArray,parameters(1),[0.1;parameters(2);parameters(3)],parameters(4),weightingMat);
%     jac = jacobianest(func,estt);
% toc
% %     %     
%     m4 = -inv(jac'*jac)*jac';
%     varcovMatrix1 = inv(jac'*weightingMat*jac)*(jac'*weightingMat*m1*m1'*weightingMat*jac)*inv(jac'*weightingMat*jac);
%     varcov = m4*(m1*m1')*m4';
%     stdErrors = diag(varcovMatrix);
%     optVarcov = inv(jac'*weightingMat*jac)'
%     optStd = diag(optVarcov);
%     
%% Bootstrapping standard errors
parfor i = 1:20
    j = datasample(mChoice,Param.N,1);
    [dataMoments,otherInfo] = momentcond(j,Param);
    
    weightingMat = eye(46);
    objectiveFun = @(estimates)properMoments(j,Param,otherInfo,lambdaArray,...
        capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],estimates(4),weightingMat);
    startvalues = [15,1.5,0.65,0.7];
    lb = [-Inf,-Inf,0,0];
    ub = [Inf,Inf,Inf,1];
    options = optimset('Display','iter','UseParallel',true);
    tic
    bootEst(:,i) = fmincon(objectiveFun,startvalues,[],[],[],[],lb,ub,[],options);
    toc
end
% %
%     end
%
%     estimatted(:,:,i) = fminconEstimatesTrial;
%     objFunVal(:,:,i) = objFunValTrial;
%%

% omegaGrid = [6:1:25]';
% startvalues = [1.5,0.65,0.6];
% options = optimset('Display','iter','TolFun',1e-10,'TolX',1e-6,...
%     'Algorithm','interior-point','AlwaysHonorConstraints','bounds');
% lb = [0,0,0];
% ub = [Inf,Inf,1];
%
% parfor i = 1:length(omegaGrid)
%     objFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omegaGrid(i),[0.1;estimates(1);estimates(2)],estimates(3));
%     [fminconEstimatesProper1(:,i),objFunValProper1(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
% end
%
%
discountFactorGrid = [0.1:0.05:1];
startvalues = [15,1.5,0.65];
lb = [-Inf,-Inf,0];
ub = [Inf,Inf,Inf];

parfor i = 1:length(discountFactorGrid)
    objFun = @(estimates)properMoments(mChoice,Param,otherInfo,lambdaArray,...
    capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],discountFactorGrid(i),weightingMat);
    [fminconEstimatesProper3(:,i),objFunValProper3(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
end

% discountFactorGrid2 = [0.1:0.01:1];
    parfor i = 1:length(discountFactorGrid)
        evall2(i) = properMoments(mChoice,Param,otherInfo,lambdaArray,...
            capRArray,mArray,fminconEstimatesProper3(1,10),[0.1;fminconEstimatesProper3(2,10);fminconEstimatesProper3(3,10)],discountFactorGrid(i),weightingMat);
    end

parfor i = 1:length(discountFactorGrid)
        evall3(i) = properMoments(mChoice,Param,otherInfo,lambdaArray,...
            capRArray,mArray,fminconEstimatesProper3(1,5),[0.1;fminconEstimatesProper3(2,5);fminconEstimatesProper3(3,5)],discountFactorGrid(i),weightingMat);
    end    
lambdaArray1 = squeeze(lambdaArray(:,:,1));
capRArray1 = squeeze(capRArray(:,:,1));

for i = 1:5    
    [mChoice((i-1)*Param.N + 1:i*Param.N,:),realHealthShocks] = dgp2(lambdaArray1,capRArray1,mOptimalAllPeriods,...
        Param);
end
[dataMoments,otherInfo] = momentcond(mChoice,Param);


weightingMat = eye(90);
objectiveFun = @(estimates)properMoments(mChoice,Param,otherInfo,lambdaArray,...
    capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],estimates(4),weightingMat);
startvalues = [15,1.5,0.65,0.7];
lb = [-Inf,-Inf,0,0];
ub = [Inf,Inf,Inf,1];
options = optimset('Display','iter','UseParallel',true);
estt2 = fmincon(objectiveFun,startvalues,[],[],[],[],lb,ub,[],options);
% weightingMat = eye(90);
% DFGrid = [0.1:0.05:1];
% parfor i = 1:length(DFGrid)
%     evaluation(i) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omega,parLambdaProcess,DFGrid(i),weightingMat);
% end
% end
%
% parfor i = 1:length(omegaGrid)
%     evaluation1(i) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omegaGrid(i),[0.1;fminconEstimatesProper1(1,20);fminconEstimatesProper1(2,20)],fminconEstimatesProper1(3,20));
% end
%
% discountFactorGrid = [0.1:0.01:0.98];
%
% parfor i = 1:length(discountFactorGrid)
%     evaluation2(i) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,fminconEstimatesProper3(1,9),[0.1;fminconEstimatesProper3(2,9);fminconEstimatesProper3(3,9)],discountFactorGrid(i));
% end
%
% parfor i = 1:length(discountFactorGrid)
%     evaluation3(i) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,fminconEstimatesProper3(1,10),[0.1;fminconEstimatesProper3(2,10);fminconEstimatesProper3(3,10)],discountFactorGrid(i));
% end
%
% evaluate1 = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,20,parLambdaProcess,0.4);
% %% Profiling Using Proper Moments.
%
% omegaGrid = [6:1:25]';
% startvalues = [1.5,0.65,0.6];
% options = optimset('Display','iter','TolFun',1e-10,'TolX',1e-6,...
%     'Algorithm','interior-point','AlwaysHonorConstraints','bounds');
% lb = [0,0,0];
% ub = [Inf,Inf,1];
% 
% parfor i = 1:length(omegaGrid)
%     objFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omegaGrid(i),[0.1;estimates(1);estimates(2)],estimates(3));
%     [fminconEstimatesProper1(:,i),objFunValProper1(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
% end
% 
% [~,index] = min(objFunValProper);
% omegaGrid = [5:0.25:30];
% parfor i = 1:length(omegaGrid)
%     valueFunc(i) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omegaGrid(i),[0.1;fminconEstimatesProper1(1,14);fminconEstimatesProper1(2,14)],fminconEstimatesProper1(3,14));
% end
% 
% plot(omegaGrid,objFunValProper)
% hold on
% plot(omegaGrid,valueFunc,'--')
% tic
% evaluate1 = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omegaGrid(93),[0.1,fminconEstimatesProper1(1,14),fminconEstimatesProper1(2,14)],fminconEstimatesProper1(3,14));
% toc
% evaluate2 = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omegaGrid(94),[0.1,fminconEstimatesProper1(1,14),fminconEstimatesProper1(2,14)],fminconEstimatesProper1(3,14));
% 
% parfor i = 1:length(omegaGrid)
%     valueFunc1(i) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omegaGrid(i),[0.1;fminconEstimatesProper(1,81);fminconEstimatesProper(2,81)],fminconEstimatesProper(3,81));
% end
% 
% discountFactorGrid = [0.1:0.05:1];
% startvalues = [15,1.5,0.65];
% lb = [-Inf,-Inf,0];
% ub = [Inf,Inf,Inf];
% 
% parfor i = 1:length(discountFactorGrid)
%     objFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],discountFactorGrid(i));
%     [fminconEstimatesProper3(:,i),objFunValProper3(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
% end
% 
% meanHealthShockGrid = [0.5:0.1:4];
% startvalues = [15,0.65,0.65];
% lb = [-Inf,0,0];
% ub = [Inf,Inf,1];
% 
% parfor i = 1:length(meanHealthShockGrid)
%     objFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,estimates(1),[0.1;meanHealthShockGrid(i);estimates(2)],estimates(3));
%     [fminconEstimatesProper(:,i),objFunValProper1(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
% end
% 
% varianceGrid = [0.5:0.1:3];
% startvalues = [15,1.5,0.65];
% lb = [-Inf,-Inf,0];
% ub = [Inf,Inf,1];
% 
% parfor i = 1:length(varianceGrid)
%     objFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,estimates(1),[0.1;estimates(2);varianceGrid(i)],estimates(3));
%     [fminconEstimatesProper(:,i),objFunValProper2(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
% end
% 
% % %% This section profiles the objective function.
% % 
% % % y and weightingMat are two inputs specific to the function newWay.
% % 
% % omegaGrid = [5:0.25:30]';
% % weightingMat = [];
% % y = 12;
% % 
% % parfor k = 1:length(omegaGrid)
% %     objFunVal5(k) = newWay2(mChoice,lambda,Param,dataMoments,uniformDrawsLambda,...
% %     lambdaArray,capRArray,mArray,omegaGrid(k),parLambdaProcess,discountFactor,weightingMat,otherInfo,y);
% % end
% % 
% % parfor k = 1:length(omegaGrid)
% %     objFunVal5000(k) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
% %     capRArray,mArray,omegaGrid(k),parLambdaProcess,discountFactor);
% % end
% % 
% % y = 3; 
% % parfor k = 1:length(omegaGrid)
% %     objFunVal15(k) = newWay2(mChoice,lambda,Param,dataMoments,uniformDrawsLambda,...
% %     lambdaArray,capRArray,mArray,omegaGrid(k),parLambdaProcess,discountFactor,weightingMat,otherInfo,y);
% % end
% % 
% % y =8;
% % parfor k = 1:length(omegaGrid)
% %     objFunVal51(k) = newWay2(mChoice,lambda,Param,dataMoments,uniformDrawsLambda,...
% %     lambdaArray,capRArray,mArray,omegaGrid(k),parLambdaProcess,discountFactor,weightingMat,otherInfo,y);
% % end
% % %% Profiling for Omega
% % 
% % omegaGrid = [5:1:30];
% % startvalues = [1.5,0.65,0.65];
% % y = 11;
% % options = optimset('TolFun',1e-10,'TolX',1e-6,...
% %     'Algorithm','interior-point','AlwaysHonorConstraints','bounds','UseParallel',false);
% % lb = [-Inf,0,0];
% % ub = [Inf,Inf,1];
% % 
% % parfor i = 1:length(omegaGrid)
% %     objFun = @(estimates)newWay2(mChoice,lambda,Param,dataMoments,uniformDrawsLambda,...
% %             lambdaArray,capRArray,mArray,omegaGrid(i),[0.1;estimates(1);estimates(2)],estimates(3),weightingMat,otherInfo,y);
% %     [fminconEstimates110(:,i),objFunVal110(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
% % end
% 
% %% I will try out a sequential optimizing procedure here.
% 
% % The plan is to use something similar to DellaVigna, List and Malmendier
% % (2012). Run fmincon many times, choosing a starting value that is drawn
% % from a uniform distribution across the permitted parameter space. 
% options = optimset('TolFun',1e-10,'TolX',1e-6,...
%     'Algorithm','interior-point','AlwaysHonorConstraints','bounds');
% lb = [-Inf,-Inf,0];
% ub = [Inf,Inf,1];
% startvalues = [15,1.5,0.6];
% 
% 
% objFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],estimates(4));
% parfor j = 1:24
%     startvalues = [5+rand,1.5+rand,0.65+rand,...
%         rand];
%     [seqFminconEstimates4(:,j),seqObjFunVal4(:,j)] = knitromatlab(objFun,...
%     startvalues,[],[],[],[],lb,ub,[],[],options);
% end 
% 
% omegaGrid = [5:0.25:30]';
% startvalues = [1.5,0.65,0.6];
% options = optimset('TolFun',1e-10,'TolX',1e-6,...
%     'Algorithm','interior-point','AlwaysHonorConstraints','bounds');
% lb = [-Inf,0,0];
% ub = [Inf,Inf,1];
% 
% parfor i = 1:length(omegaGrid)
%     objFun = @(estimates)properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omegaGrid(i),[0.1;estimates(1);estimates(2)],estimates(3));
%     [knitroEstimatesProper(:,i),objFunValKnitro(:,i)] = knitromatlab(objFun,startvalues,[],[],[],[],lb,ub,[],[],options);
% end
% % 
% % %% This one is profiling for other parameters
% % 
% % discountFactorGrid = [0.1:0.05:1];
% % startvalues = [15,0.5,0.65];
% % weightingMat = [];
% % y = 10;
% % options = optimset('TolFun',1e-10,'TolX',1e-6,...
% %     'Algorithm','interior-point','AlwaysHonorConstraints','bounds','UseParallel',false);
% % lb = [-Inf,-Inf,0];
% % ub = [Inf,Inf,Inf];
% %     
% % parfor i = 1:length(discountFactorGrid)
% %     objFun = @(estimates)newWay2(mChoice,lambda,Param,dataMoments,uniformDrawsLambda,...
% %             lambdaArray,capRArray,mArray,estimates(1),[0.1;estimates(2);estimates(3)],discountFactorGrid(i),weightingMat,otherInfo,y);
% %     [fminconEstimates1(:,i),objFunVal1(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
% % end
% % 
% % startvalues = [15,0.5,0.65]; 
% % lb = [-Inf,0,0];
% % ub = [Inf,Inf,1];
% % meanHealthShockGrid = [0.1:0.1:5];
% % 
% % parfor i = 1:length(meanHealthShockGrid)
% %     objFun = @(estimates)newWay2(mChoice,lambda,Param,dataMoments,uniformDrawsLambda,...
% %             lambdaArray,capRArray,mArray,estimates(1),[0.1;meanHealthShockGrid(i);estimates(2)],estimates(3),weightingMat,otherInfo,y);
% %     [fminconEstimates2(:,i),objFunVal2(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
% % end
% % 
% % startvalues = [15,0.5,0.65];
% % lb = [-Inf,-Inf,0];
% % ub = [Inf,Inf,1];
% % varianceHealthShockGrid = [0.1:0.1:4];
% % 
% % parfor i = 1:length(varianceHealthShockGrid)
% %     objFun = @(estimates)newWay2(mChoice,lambda,Param,dataMoments,uniformDrawsLambda,...
% %             lambdaArray,capRArray,mArray,estimates(1),[0.1;estimates(2);varianceHealthShockGrid(i)],estimates(3),weightingMat,otherInfo,y);
% %     [fminconEstimates3(:,i),objFunVal3(:,i)] = fmincon(objFun,startvalues,[],[],[],[],lb,ub,[],options);
% % end
% % 
% 
% %% This section evaluates the objective function for given delta
% 
% omegaGrid = [5:0.5:30];
% meanGrid = [0.5:0.25:3.5];
% varianceGrid = [0.5:0.1:1.5];
% objFunValGivenDelta = NaN(length(omegaGrid),length(meanGrid),length(varianceGrid));
% 
% for k = 1:length(varianceGrid)
%     for j = 1:length(meanGrid)
%         parfor i = 1:length(omegaGrid)
%             objFunValGivenDelta(i,j,k) = properMoments(mChoice,Param,otherInfo,uniformDrawsLambda,lambdaArray,...
%     capRArray,mArray,omegaGrid(i),[0.1;meanGrid(j);varianceGrid(k)],0.95);
%         end
%     end
% end
% 
% [oArray,meanArray,varArray] = ndgrid(omegaGrid,meanGrid,varianceGrid);
  