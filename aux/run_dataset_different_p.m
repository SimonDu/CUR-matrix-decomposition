function savedata = run_dataset_different_p(in,methods,p_values)
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
% - c : number of columns
% - r : number of rows
% - q, the number of times to repeat each random algorithms
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
if wantq('deterministic')
    fprintf('...deterministic\n');
    savedata.deterministic.froerr = zeros(1,length(p_values));
    savedata.deterministic.froerr_k = zeros(1,length(p_values));
%     savedata.deterministic.specerr = zeros(1,length(p_values));
%     savedata.deterministic.specerr_k = zeros(1,length(p_values));
    savedata.deterministic.sigma_k = zeros(1,length(p_values));
    for i =1:length(p_values)
        in.p = p_values(i);
        fprintf('...running p=%d\n',in.p);
        output = deterministic(in);
        savedata.deterministic.froerr(i) = mean(output.froerr);
        savedata.deterministic.froerr_k(i) = mean(output.froerr_k);
%         savedata.deterministic.specerr(i) = mean(output.specerr);
%         savedata.deterministic.specerr_k(i) = mean(output.froerr_k);
        savedata.deterministic.sigma_k(i) = mean(output.sigma_k);
    end
end


if wantq('subspace_expected')
    fprintf('...subspace_expected\n');  
    savedata.subspace_expected.froerr = zeros(1,length(p_values));
    savedata.subspace_expected.froerr_k = zeros(1,length(p_values));
%     savedata.subspace_expected.specerr = zeros(1,length(p_values));
%     savedata.subspace_expected.specerr_k = zeros(1,length(p_values));
    savedata.subspace_expected.sigma_k = zeros(1,length(p_values));
    for i =1:length(p_values)
        in.p = p_values(i);
        fprintf('...running p=%d\n',in.p);
        output = subspace_expected(in);
        savedata.subspace_expected.froerr(i) = mean(output.froerr);
        savedata.subspace_expected.froerr_k(i) = mean(output.froerr_k);
%         savedata.subspace_expected.specerr(i) = mean(output.specerr);
%         savedata.subspace_expected.specerr_k(i) = mean(output.froerr_k);
        savedata.subspace_expected.sigma_k(i) = mean(output.sigma_k);
    end
end



if wantq('subspace_approxlevscore_gaussian')
    fprintf('...subspace_approxlevscore_gaussian\n');  
    savedata.subspace_approxlevscore_gaussian.froerr = zeros(1,length(p_values));
    savedata.subspace_approxlevscore_gaussian.froerr_k = zeros(1,length(p_values));
%     savedata.subspace_approxlevscore_gaussian.specerr = zeros(1,length(p_values));
%     savedata.subspace_approxlevscore_gaussian.specerr_k = zeros(1,length(p_values));
    savedata.subspace_approxlevscore_gaussian.sigma_k = zeros(1,length(p_values));
    for i =1:length(p_values)
        in.p = p_values(i);
        fprintf('...running p=%d\n',in.p);
        output = subspace_approxlevscore_gaussian(in);
        savedata.subspace_approxlevscore_gaussian.froerr(i) = mean(output.froerr);
        savedata.subspace_approxlevscore_gaussian.froerr_k(i) = mean(output.froerr_k);
%         savedata.subspace_approxlevscore_gaussian.specerr(i) = mean(output.specerr);
%         savedata.subspace_approxlevscore_gaussian.specerr_k(i) = mean(output.froerr_k);
        savedata.subspace_approxlevscore_gaussian.sigma_k(i) = mean(output.sigma_k);
    end
end


if wantq('subspace_approxlevscore_powermethod')
    fprintf('...subspace_approxlevscore_powermethod\n');  
    savedata.subspace_approxlevscore_powermethod.froerr = zeros(1,length(p_values));
    savedata.subspace_approxlevscore_powermethod.froerr_k = zeros(1,length(p_values));
%     savedata.subspace_approxlevscore_powermethod.specerr = zeros(1,length(p_values));
%     savedata.subspace_approxlevscore_powermethod.specerr_k = zeros(1,length(p_values));
    savedata.subspace_approxlevscore_powermethod.sigma_k = zeros(1,length(p_values));
    for i =1:length(p_values)
        in.p = p_values(i);
        fprintf('...running p=%d\n',in.p);
        output = subspace_approxlevscore_powermethod(in);
        savedata.subspace_approxlevscore_powermethod.froerr(i) = mean(output.froerr);
        savedata.subspace_approxlevscore_powermethod.froerr_k(i) = mean(output.froerr_k);
%         savedata.subspace_approxlevscore_powermethod.specerr(i) = mean(output.specerr);
%         savedata.subspace_approxlevscore_powermethod.specerr_k(i) = mean(output.froerr_k);
        savedata.subspace_approxlevscore_powermethod.sigma_k(i) = mean(output.sigma_k);
    end
end



% save(datasetfname, 'savedata');

end
