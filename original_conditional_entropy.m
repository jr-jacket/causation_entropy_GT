function [ conditionalEntropy ] = original_conditional_entropy(X)
[m,n] = size(X);
h_full =(4/(m+2))^(1/(m+4))*n^(-1/(m+4));  
h_subset =(4/((m-1)+2))^(1/((m-1)+4))*n^(-1/((m-1)+4));     % kernel width for X_subset and X_full

X_subset=X(2:end,:);    %determine data for two joint probabilities needed for relevant conditional probability
X2_full=X;

sumEstimator = 0;
for i = 1:n
    pyz = p_mkde(X_subset(:,i),X_subset,h_subset);
    pxyz = p_mkde(X2_full(:,i),X2_full,h_full);
    sumEstimator = sumEstimator + log(pyz/pxyz);
    
end
conditionalEntropy = sumEstimator/n;

end
