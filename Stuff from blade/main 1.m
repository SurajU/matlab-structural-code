% structural model, main file


%% clean up and set seed, set number random draws
clear;
rng('default');
rng(1);

% number of random draws to be used for health process
numberDrawsLambda = 300;


%% parameters for figures

set(gca,'fontsize',24); %to be checked

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);



%% define parameters, some of them will be estimated later

% number of time periods within one year
capT = 12;

% deductible
capD = 375;

% moral hazard parameter
omega = 100;

% discount factor
discountFactor = 0.9;

% health process
% first element: parameter for probability that it is zero, modeled
%                as 1/(1+exp(-parLambdaProcess(1))
% second and third element: mean and variance of underlying normal
%                           distribution if positive

parLambdaProcess = [2.5;2;4];


%% define grids

% health care need
lambda = [(0:10:1000)';(5000:5000:25000)';50000];

% remaining deductible
capR = (0:10:400)';

% health care expenditures
m = [(0:5:1000)';(5000:5000:25000)';50000;100000];

% dimensions
A = size(lambda,1);
B = size(capR,1);
C = size(m,1);

%   we start with the following variables
%     lambda     Ax1    a grid for lambda
%     capR       Bx1    a grid for the remaining deductible
%     m          Cx1    a grid for the chosen health care expenditures
%     omega      1x1    parameter omega

% now we construct an aligned grid of dimension AxBxC
[lambdaArray,capRArray,mArray] = ndgrid(lambda,capR,m);


%% produce random draws

% draws for health process
healthProcessDimensionsRandomness = 2;
uniformDrawsLambda = rand(numberDrawsLambda,B,capT,healthProcessDimensionsRandomness);


%% solve model

% initialize
vAllPeriods = NaN(A,B,capT+1); %period T+1 will contain terminal value
vChoiceSpecificAllPeriods = NaN(A,B,C,capT);
mOptimalAllPeriods = NaN(A,B,capT);

% health process does not depend on previous period
% can be generalized
lambdaPreviousPeriod = [];

% terminal value
vAllPeriods(:,:,capT+1) = 0;

% backward recursion
for t=capT:-1:1
    
    vNext = repmat(vAllPeriods(:,:,t+1),1,1,C);
    
    [vAllPeriods(:,:,t),vChoiceSpecificAllPeriods(:,:,:,t),mOptimalAllPeriods(:,:,t)] ...
        = valueFunction(lambdaArray,capRArray,...
        mArray,omega,discountFactor,uniformDrawsLambda,parLambdaProcess,vNext,...
        lambdaPreviousPeriod,t);
    
end


%% present results

%figure plotting m against R and lambda, at different times
for counter=1:4
    whichT = counter*3;
    subplot(2,2,counter)
    X = squeeze(lambdaArray(1:51,:,1));
    Y = squeeze(capRArray(1:51,:,1));
    Z = squeeze(mOptimalAllPeriods(1:51,:,whichT));
    surf(X,Y,Z)
    xlabel('\lambda_t')
    ylabel('R_t')
    zlabel('m_t')
    titleString = ['t=' num2str(whichT)];
    title(titleString)
end
% print('-dpdf','../../output/figureMatlabOptimalMContour.pdf')

%figure plotting choice specific value function for given R and lambda
clf
whichCapR = 41;
Y = [squeeze(mOptimalAllPeriods(1:51,whichCapR,12)+12) ...
        squeeze(mOptimalAllPeriods(1:51,whichCapR,11)+11) ...
        squeeze(mOptimalAllPeriods(1:51,whichCapR,10)+10) ...
        squeeze(mOptimalAllPeriods(1:51,whichCapR,9)+9) ...
        squeeze(mOptimalAllPeriods(1:51,whichCapR,8)+8) ...
        squeeze(mOptimalAllPeriods(1:51,whichCapR,7)+7) ...
        squeeze(mOptimalAllPeriods(1:51,whichCapR,6)+6)];
plot(lambda(1:51),Y);
xlabel('\lambda')
ylabel('m when R=400')
print('-dpdf','../../output/figureMatlabOptimalMAgainstLambda.pdf');

%plot(m(1:201),squeeze(vChoiceSpecificAllPeriods(4,4,1:201,3)))
%plot(lambda(1:101),squeeze(mOptimalAllPeriods(1:101,4,12)))

