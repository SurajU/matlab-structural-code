function drawsLambda = simulateLambda(lambdaPreviousPeriod,uniformDrawsLambda,parLambdaProcess,t)
% Produces draws of lambda in the next period
%
%   INPUTS
%     lambdaPreviousPeriod     Ax1      state variable, for now not used
%     uniformDrawsLambda       DxYxZ    where Y=capT
%                                       and Z=healthProcessDimensionsRandomness
%     parLambdaProcess         Ex1
%     t                        1x1      indicator for time period
%
%   OUTPUTS
%     drawsLambda              AxD  

probNoNeed = 1/(1+exp(-parLambdaProcess(1)));

indicatorNoNeed = uniformDrawsLambda(:,:,t,1) <= probNoNeed;  % 1 because use first dimension of random draws

% latent need is square of a normal distribution
normalDraws = norminv(uniformDrawsLambda(:,:,t,2),parLambdaProcess(2),parLambdaProcess(3));
latentNeed = normalDraws.^2;

% draw of lambda is either zero or latent need
drawsLambda = indicatorNoNeed*0 + (1-indicatorNoNeed).*latentNeed;

end

