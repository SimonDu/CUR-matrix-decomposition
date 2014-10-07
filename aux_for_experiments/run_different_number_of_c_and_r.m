function savedata = run_different_number_of_c_and_r(in,methods,c_values,r_values)
%
% savedata = generate_dataset(in)
%
% Stores the error and timing information from running chosen Nystrom
% methods on a given dataset with given parameters using a specified range
% of column samples
%
% in is a structure with (at least) the following fields:
% - descr, a string describing this dataset
% - datasetbasename, a base filename for this dataset
% - datasetdir, the directory to store this dataset in
% - methods, a cell array of the Nystrom methods to apply; see below for
%   valid methods
% - A, a matrix
% - k, the target rank of the approximation
% - p, a vector of the oversampling parameter for deterministic and
% subspace sampling method
% - q, the number of times to repeat each random algorithms
% - sigma_k, 1 if we want output contain sigma_k, 0 otherwise
% - froerr, 1 if we want output contain froerr, 0 otherwise
% - froerr_k, 1 if we want output contain froerr_k, 0 otherwise
% - specerr, 1 if we want output contain specerr, 0 otherwise
% - specerr_k, 1 if we want output contain specerr_k, 0 otherwise

%
% Other fields may be required to be present, depending on the Nystrom
% methods specified. Valid methods (case insensitive):
% 'simple' - uniform column sampling without replacement
% 'srft' - SRFT mixture-based
% 'gaussian' - Gaussian mixture-based
% 'levscore' - leverage score based column sampling
% 'froblev' - frobenius sketch approximate leverage score sampling
%  requires additional field vareps, a number between 0 and 1
% 'speclev' - spectral sketch approximate leverage score sampling
%  requires additional fields vareps and chunk, how many iterations to take
%  before reorthogonalization
% 'approxlev' - power method approximate leverage score sampling
%  requires additional fields vareps and chunk
% 'mixedprobs' - sampling based on mixture of uniform and levscore
%  requires addtional field p, a number between 0 and 1 that indicates how
%  much of the mixture is the leverage score probabilities
%
% The return struct savedata contains all the information saved from the
% Nystrom experiments run.
%

in.methods = lower(methods);
wantq = @(methodname) any(strcmp(methodname, methods));
savedata.in = in;

% store the outputs for each methods

if wantq('near_optimal')
    fprintf('...near_optimal\n');
    if(in.froerr)
        savedata.near_optimal.froerr = zeros(1,length(c_values));
    end
    if(in.froerr_k)
        savedata.near_optimal.froerr_k = zeros(1,length(c_values));
    end
    if(in.specerr)
        savedata.near_optimal.specerr = zeros(1,length(c_values));
    end
    if(in.specerr_k)
        savedata.near_optimal.specerr_k = zeros(1,length(c_values));
    end
    if(in.sigma_k)
        savedata.near_optimal.sigma_k = zeros(1,length(c_values));
    end
    for i =1:length(c_values)
        in.c = c_values(i);
        in.r = r_values(i);
        fprintf('...running c=%d, r = %d using subspace near_optimal\n',in.c, in.r);
        output = near_optimal(in);
        if(in.froerr)
            savedata.near_optimal.froerr(i) = mean(output.froerr);
        end
        if(in.froerr_k)
            savedata.near_optimal.froerr_k(i) = mean(output.froerr_k);
        end
        if(in.specerr)
            savedata.near_optimal.specerr(i) = mean(output.specerr);
        end
        if(in.specerr_k)
            savedata.near_optimal.specerr_k(i) = mean(output.specerr_k);
        end
        if(in.sigma_k)
            savedata.near_optimal.sigma_k(i) = mean(output.sigma_k);
        end
    end
end


if wantq('subspace_approxlevscore_gaussian')
    fprintf('...subspace_approxlevscore_gaussian\n');
    if(in.froerr)
        savedata.subspace_approxlevscore_gaussian.froerr = zeros(1,length(c_values));
    end
    if(in.froerr_k)
        savedata.subspace_approxlevscore_gaussian.froerr_k = zeros(1,length(c_values));
    end
    if(in.specerr)
        savedata.subspace_approxlevscore_gaussian.specerr = zeros(1,length(c_values));
    end
    if(in.specerr_k)
        savedata.subspace_approxlevscore_gaussian.specerr_k = zeros(1,length(c_values));
    end
    if(in.sigma_k)
        savedata.subspace_approxlevscore_gaussian.sigma_k = zeros(1,length(c_values));
    end
    for i =1:length(c_values)
        in.c = c_values(i);
        in.r = r_values(i);
        fprintf('...running c=%d, r=%d using approxlevscore_gaussian subspace sampling\n',in.c,in.r);
        output = subspace_approxlevscore_gaussian(in);
        if(in.froerr)
            savedata.subspace_approxlevscore_gaussian.froerr(i) = mean(output.froerr);
        end
        if(in.froerr_k)
            savedata.subspace_approxlevscore_gaussian.froerr_k(i) = mean(output.froerr_k);
        end
        if(in.specerr)
            savedata.subspace_approxlevscore_gaussian.specerr(i) = mean(output.specerr);
        end
        if(in.specerr_k)
            savedata.subspace_approxlevscore_gaussian.specerr_k(i) = mean(output.specerr_k);
        end
        if(in.sigma_k)
            savedata.subspace_approxlevscore_gaussian.sigma_k(i) = mean(output.sigma_k);
        end
    end
end

if wantq('subspace_approxlevscore_powermethod')
    fprintf('...subspace_approxlevscore_powermethod\n');
    if(in.froerr)
        savedata.subspace_approxlevscore_powermethod.froerr = zeros(1,length(c_values));
    end
    if(in.froerr_k)
        savedata.subspace_approxlevscore_powermethod.froerr_k = zeros(1,length(c_values));
    end
    if(in.specerr)
        savedata.subspace_approxlevscore_powermethod.specerr = zeros(1,length(c_values));
    end
    if(in.specerr_k)
        savedata.subspace_approxlevscore_powermethod.specerr_k = zeros(1,length(c_values));
    end
    if(in.sigma_k)
        savedata.subspace_approxlevscore_powermethod.sigma_k = zeros(1,length(c_values));
    end
    for i =1:length(c_values)
        in.c = c_values(i);
        in.r = r_values(i);
        fprintf('...running c=%d. r=%d using approxlevscore_powermethod subspace sampling\n',in.c,in.r);
        output = subspace_approxlevscore_powermethod(in);
        if(in.froerr)
            savedata.subspace_approxlevscore_powermethod.froerr(i) = mean(output.froerr);
        end
        if(in.froerr_k)
            savedata.subspace_approxlevscore_powermethod.froerr_k(i) = mean(output.froerr_k);
        end
        if(in.specerr)
            savedata.subspace_approxlevscore_powermethod.specerr(i) = mean(output.specerr);
        end
        if(in.specerr_k)
            savedata.subspace_approxlevscore_powermethod.specerr_k(i) = mean(output.specerr_k);
        end
        if(in.sigma_k)
            savedata.subspace_approxlevscore_powermethod.sigma_k(i) = mean(output.sigma_k);
        end
    end
end

if wantq('deterministic')
    fprintf('...deterministic\n');
    if(in.froerr)
        savedata.deterministic.froerr = zeros(1,length(c_values));
    end
    if(in.froerr_k)
        savedata.deterministic.froerr_k = zeros(1,length(c_values));
    end
    if(in.specerr)
        savedata.deterministic.specerr = zeros(1,length(c_values));
    end
    if(in.specerr_k)
        savedata.deterministic.specerr_k = zeros(1,length(c_values));
    end
    if(in.sigma_k)
        savedata.deterministic.sigma_k = zeros(1,length(c_values));
    end
    for i =1:length(c_values)
        in.c = c_values(i);
        in.r = r_values(i);
        fprintf('...running c=%d,r=%d using deterministic\n',in.c,in.r);
        output = deterministic(in);
        if(in.froerr)
            savedata.deterministic.froerr(i) = mean(output.froerr);
        end
        if(in.froerr_k)
            savedata.deterministic.froerr_k(i) = mean(output.froerr_k);
        end
        if(in.specerr)
            savedata.deterministic.specerr(i) = mean(output.specerr);
        end
        if(in.specerr_k)
            savedata.deterministic.specerr_k(i) = mean(output.specerr_k);
        end
        if(in.sigma_k)
            savedata.deterministic.sigma_k(i) = mean(output.sigma_k);
        end
    end
end

if wantq('subspace_expected')
    fprintf('...subspace_expected\n');
    if(in.froerr)
        savedata.subspace_expected.froerr = zeros(1,length(c_values));
    end
    if(in.froerr_k)
        savedata.subspace_expected.froerr_k = zeros(1,length(c_values));
    end
    if(in.specerr)
        savedata.subspace_expected.specerr = zeros(1,length(c_values));
    end
    if(in.specerr_k)
        savedata.subspace_expected.specerr_k = zeros(1,length(c_values));
    end
    if(in.sigma_k)
        savedata.subspace_expected.sigma_k = zeros(1,length(c_values));
    end
    for i =1:length(c_values)
        in.c = c_values(i);
        in.r = r_values(i);
        fprintf('...running c=%d, r = %d using subspace subspace_expected\n',in.c, in.r);
        output = subspace_expected(in);
        if(in.froerr)
            savedata.subspace_expected.froerr(i) = mean(output.froerr);
        end
        if(in.froerr_k)
            savedata.subspace_expected.froerr_k(i) = mean(output.froerr_k);
        end
        if(in.specerr)
            savedata.subspace_expected.specerr(i) = mean(output.specerr);
        end
        if(in.specerr_k)
            savedata.subspace_expected.specerr_k(i) = mean(output.specerr_k);
        end
        if(in.sigma_k)
            savedata.subspace_expected.sigma_k(i) = mean(output.sigma_k);
        end
    end
end


if wantq('uniform_sampling')
    fprintf('...uniform_sampling\n');
    if(in.froerr)
        savedata.uniform_sampling.froerr = zeros(1,length(c_values));
    end
    if(in.froerr_k)
        savedata.uniform_sampling.froerr_k = zeros(1,length(c_values));
    end
    if(in.specerr)
        savedata.uniform_sampling.specerr = zeros(1,length(c_values));
    end
    if(in.specerr_k)
        savedata.uniform_sampling.specerr_k = zeros(1,length(c_values));
    end
    if(in.sigma_k)
        savedata.uniform_sampling.sigma_k = zeros(1,length(c_values));
    end
    for i =1:length(c_values)
        in.c = c_values(i);
        in.r = r_values(i);
        fprintf('...running c=%d, r = %d using uniform_sampling\n',in.c, in.r);
        output = uniform_sampling(in);
        if(in.froerr)
            savedata.uniform_sampling.froerr(i) = mean(output.froerr);
        end
        if(in.froerr_k)
            savedata.uniform_sampling.froerr_k(i) = mean(output.froerr_k);
        end
        if(in.specerr)
            savedata.uniform_sampling.specerr(i) = mean(output.specerr);
        end
        if(in.specerr_k)
            savedata.uniform_sampling.specerr_k(i) = mean(output.specerr_k);
        end
        if(in.sigma_k)
            savedata.uniform_sampling.sigma_k(i) = mean(output.sigma_k);
        end
    end
end

if wantq('randomized_unweighted')
    fprintf('...randomized_unweighted\n');
    if(in.froerr)
        savedata.randomized_unweighted.froerr = zeros(1,length(c_values));
    end
    if(in.froerr_k)
        savedata.randomized_unweighted.froerr_k = zeros(1,length(c_values));
    end
    if(in.specerr)
        savedata.randomized_unweighted.specerr = zeros(1,length(c_values));
    end
    if(in.specerr_k)
        savedata.randomized_unweighted.specerr_k = zeros(1,length(c_values));
    end
    if(in.sigma_k)
        savedata.randomized_unweighted.sigma_k = zeros(1,length(c_values));
    end
    for i =1:length(c_values)
        in.c = c_values(i);
        in.r = r_values(i);
        fprintf('...running c=%d, r = %d using subspace randomized_unweighted\n',in.c, in.r);
        output = randomized_unweighted(in);
        if(in.froerr)
            savedata.randomized_unweighted.froerr(i) = mean(output.froerr);
        end
        if(in.froerr_k)
            savedata.randomized_unweighted.froerr_k(i) = mean(output.froerr_k);
        end
        if(in.specerr)
            savedata.randomized_unweighted.specerr(i) = mean(output.specerr);
        end
        if(in.specerr_k)
            savedata.randomized_unweighted.specerr_k(i) = mean(output.specerr_k);
        end
        if(in.sigma_k)
            savedata.randomized_unweighted.sigma_k(i) = mean(output.sigma_k);
        end
    end
end



end
