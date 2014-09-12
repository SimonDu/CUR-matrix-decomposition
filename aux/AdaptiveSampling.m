function [idx] = AdaptiveSampling(prob, r)

n = length(prob);

iter = 0;
flag = false;
while true
    iter = iter + 1;
    % Remark: here the sampling strategy is not the same as that in the paper
    %		According to the paper, one should sample one index at a time using the function "mnrnd(1, prob, 1)", and repeat it for r times.
    %		However, the MATLAB function "mnrnd" appears to have bug, so here the implementation is different from the paper.
    selected = binornd(1, min(1, r * prob));
    selected = (selected == 1);
    if sum(selected) >= r
        break;
    end
    if iter > 20
        flag = true;
        break;
    end
end

index = 1:n;
idx = index(~~selected);

if flag == false
    idx = idx(1:r); % more than r columns are selected
end

end