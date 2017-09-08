function u = flowUtility(lambdaArray,capRArray,mArray,omega)
% Computes flow utility
%
%   INPUTS
%     lambdaArray     AxBxC    a grid along the first dimension for lambda
%     capRArray       AxBxC    a grid along the second dimension for the
%                              remaining deductible
%     mArray          AxBxC    a grid along the third dimension for the 
%                              chosen health care expenditures
%     omega      1x1    parameter omega
%
%   OUTPUTS
%     u          AxBxC  an array containing the flow utiltiy for different
%                       combinations of lambda, capR and m

% out of pocket expenditures
outOfPocket = min(mArray,capRArray);

u = (mArray-lambdaArray) - (1/(2*omega)) * (mArray-lambdaArray).^2 - outOfPocket;

end

