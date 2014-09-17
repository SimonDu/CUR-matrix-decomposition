function savedata = generate_dataset(in)
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
% - number_of_c_and_r: a 2d matrix (contains column, row) pairs
% - q, the number of times to repeat each Nystrom method for each number of
%  column samples
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

in.methods = lower(in.methods);
wantq = @(methodname) any(strcmp(methodname, in.methods));
savedata.in = in;
datasetfname = fullfile(in.datasetdir, in.datasetbasename);
% compute the leverage scores and optimal rank-k approximation errors
savedata.properties = datadescription(in.A);

% store the outputs for each methods

if wantq('uniform_sampling')
    fprintf('...simple\n');
    savedata.uniform_sampling_output = uniform_sampling(in);
end

if wantq('deterministic')
    fprintf('...simple\n');
    savedata.deterministic_output = deterministic(in);
end

if wantq('subspace_expected')
    fprintf('...simple\n');
    savedata.subspace_expected_output = subspace_expected(in);
end



save(datasetfname, 'savedata');

end
