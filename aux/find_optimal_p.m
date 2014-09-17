function [ optimal_p ] = find_optimal_p(output,l)
%find the optimal p in different metrics
%output: an output struct from differnt algorithms
%l: choose lth row
%outimal_p: a struct consist of optimal values for different metrics
sigma_k = [];
froerr = [];
froerr_k = [];
specerr = [];
specerr_k = [];
p_length = size(output,2);
for i=1:p_length
    sigma_k = [sigma_k; output{l,i}.sigma_k];
    froerr = [froerr;mean(output{l,i}.froerr(1,:))];
    froerr_k = [froerr_k;mean(output{l,i}.froerr(2,:))];
    specerr = [specerr; mean(output{l,i}.specerr(1,:))];
    specerr_k = [specerr_k; mean(output{l,i}.specerr(2,:))];
end
optimal_p.sigma_k = max(sigma_k);
optimal_p.froerr = min(froerr);
optimal_p.froerr_k = min(froerr_k);
optimal_p.specerr = min(specerr);
optimal_p.specerr_k = min(specerr_k);

end

