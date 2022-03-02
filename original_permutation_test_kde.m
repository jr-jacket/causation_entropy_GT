function [threshold] = original_permutation_test_kde(percentage,perm,Xt1,Xt,S)
data_length = length(Xt);
for i = 1:perm
    temp = randperm(data_length);       % randomly permute data entries
    for j = 1:data_length               %assemble necessary vectors from randomly permuted order
        xt1(j) = Xt1(j);
        xt(j) = Xt(temp(j));
    end
    CE1 = original_conditional_entropy([xt1; S]) - original_conditional_entropy([xt1; S; xt]); %compute conditional entropy for random permutation  
    CE(i) = CE1;%save permuted CE value
end
sorted_ce = sort(CE); %Order randomly permuted CE according to magnitude
threshold = sorted_ce(round(percentage*perm));  %Find the corresponding CE value that corresponds to the statistically significant level of randomized CE 
return
