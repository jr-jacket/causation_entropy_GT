%Function takes as inputs the state trajectories (xdata), the corresponding
%time stream (if unnecessary to dynamics can be removed everywhere), and an
%anonymous function of the function space to be considered by the CEM
%(potentialFunctions) and the number of permutations to be used in the
%permutation test (perm). 

%Note that the first line evaluates and creates a matrix of the form
%[f_1(x);...f_n(x)] where each column is the value of the row at a given
%time step. If this matrix has already been computed, the inputs could be
%easily changed to replace potentialFunctions with it and set evaluatedFunctions equal to said input. 
function [rawCEM,adjustedCEM,Tmatrix]=calccem(xdata,time,potentialFunctions,perm)

evaluatedFunctions=potentialFunctions(xdata,time); %Potential function space evaluated at data points

[mfunctions,nfunctions]=size(evaluatedFunctions);   %determine number of potential functions
[mdata,ndata]=size(xdata);  %determine number of states 

dataEarly=evaluatedFunctions(:,1:end-1); %Early time shifted values of potential functions

dataLate=evaluatedFunctions(:,2:end); %Early time shifted values of potential functions

S=cell(mfunctions,1);   %define potential functions cell array

for i=1:length(S)
    for j=1:mfunctions
        if i~=j
            S{i}=[S{i};dataEarly(j,:)];  %make each entry of the cell the set of all functions except for the one to be considered S
        end
    end
end

cem=zeros(mdata,mfunctions);
Tmatrix=zeros(mdata,mfunctions);

for i=1:mdata
    parfor j=1:mfunctions
        cem(i,j)=original_conditional_entropy([dataLate(i,:);S{j}]) - original_conditional_entropy([dataLate(i,:); S{j}; dataEarly(j,:)]);
        [Tmatrix(i,j)] = original_permutation_test_kde(0.99, perm, dataLate(i,:),dataEarly(j,:), S{j});%%%%% X,Z,S
    end
end
rawCEM=cem;

for i=1:mdata
    for j=1:mfunctions
        if (Tmatrix(i,j)>=cem(i,j))|(cem(i,j)<0)
            cem(i,j)=0;                                 %Compare rawCEM values to computed Tmatrix and 0 to see if statistically relvant and larger than 0.
        end
    end
end
adjustedCEM=cem;
return