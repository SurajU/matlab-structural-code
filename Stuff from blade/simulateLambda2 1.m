function drawsLambda = simulateLambda2(lambdaPreviousPeriod,uniformDrawsLambda,parLambdaProcess,t)
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

% latent need is square of a normal distribution
normalDraws = exp(norminv(uniformDrawsLambda(:,:,t,2),parLambdaProcess(2),parLambdaProcess(3)));


% draw of lambda is either zero or latent need
drawsLambda = normalDraws;

end

