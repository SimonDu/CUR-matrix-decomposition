function savedata = run_dataset(in)
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
% - p, the rank of first two partition matrix
% - number_of_c_and_r, a two row matrix specifying the numbers of column and rows 
%   samples to use   
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
% CUR experiments run.
%

in.methods = lower(in.methods);
wantq = @(methodname) any(strcmp(methodname, in.methods));
savedata.in = in;
datasetfname = fullfile(in.datasetdir, in.datasetbasename);
% compute the leverage scores and optimal rank-k approximation errors
if in.linearkernelflag == 0 % A is a PSD matrix
    tic
       [U, Sigma] = orderedeigs(in.A, in.k+1);
       U1t = U(:, 1:in.k)';
       savedata.levscores = sum(U1t.*U1t);
    in.levscorecomputationtime = toc;
    in.levscoreprobs = savedata.levscores/in.k;
    savedata.topspectrum = diag(Sigma(1:in.k,1:in.k));
    savedata.optspecerr = Sigma(in.k+1,in.k+1);
    savedata.optfroerr = sqrt(norm(in.A, 'fro')^2 - sum(savedata.topspectrum.^2));
    savedata.opttrerr = trace(in.A) - sum(savedata.topspectrum);
else % the actual PSD matrix is AA^T
    tic
       [U, Sigma, ~] = svds(in.A, in.k+1);
       U1t = U(:, 1:in.k)';
       savedata.levscores = sum(U1t.*U1t);
    in.levscoreprobs = savedata.levscores/in.k;
    in.levscorecomputationtime = toc;
    savedata.topspectrum = diag(Sigma(1:in.k, 1:in.k)).^2;
    savedata.optspecerr = Sigma(in.k+1,in.k+1)^2;
    A = in.A*in.A';
    savedata.optfroerr = sqrt(norm(A, 'fro')^2 - sum(savedata.topspectrum.^2));
    savedata.opttrerr = trace(A) - sum(savedata.topspectrum);
end

% if we use the tall thin leverage score algorithm to compute Nystrom
% approximations, we should compare it to using the QR algorithm to compute
% the standard leverage score extension for fairness
if wantq('tallthinlevscoreapprox')
    tic
        [Q,~] = qr(in.A,0);
        in.qrlevscoreprobs = sum(Q.*Q,2)/size(in.A,2);
    in.qrlevscorecomputationtime = toc;
end

% store the errors of the specified Nystrom methods
for lidx = 1:length(in.lvals)
    in.l = in.lvals(lidx);
    
    fprintf('Evaluating errors for %s, l = %d (%d of %d values)\n', in.datasetbasename, in.l, lidx, length(in.lvals));
    
    if wantq('testmethod')
        fprintf('...test method');
        testmethodData(lidx) = test_Nystrom(in);
    end
    
    if wantq('simple')
        fprintf('...simple\n');
        simpleData(lidx) = simple_Nystrom(in);
    end
    
    if wantq('srft')
        fprintf('...srft\n');
        srftData(lidx) = srft_Nystrom(in);
    end
    
    if wantq('gaussian')
        fprintf('...gaussian\n');
        gaussianData(lidx) = gaussian_Nystrom(in);
    end
    
    if wantq('levscore')
        fprintf('...levscore\n');
        levscoreData(lidx) = levscore_Nystrom(in);
    end
    
    if wantq('froblev')
        fprintf('...froblev\n');
        froblevData(lidx) = froblev_Nystrom(in);
    end
    
    if wantq('speclev')
        fprintf('...speclev\n');
        speclevData(lidx) = speclev_Nystrom(in);
    end
    
    if wantq('approxlev')
        fprintf('...approxlev\n');
        approxlevData(lidx) = approxlev_Nystrom(in);
    end
    
    if wantq('mixedprobs')
        fprintf('...mixedprobs\n');
        mixedprobsData(lidx) = mixedprobs_Nystrom(in);
    end

    if wantq('tallthinlevscoreapprox')
        fprintf('...tallthinlevscoreapprox\n');
            tallthinData(lidx) = tallthin_Nystrom(in);
        fprintf('...qrlevscore for comparison to tallthinalg\n');
            qrlevscoreData(lidx) = qrlevscore_Nystrom(in);
    end
end

if wantq('testmethod')
    savedata.testmethodData = testmethodData;
end

if wantq('simple')
    savedata.simpleData = simpleData;
end

if wantq('srft')
    savedata.srftData = srftData;
end

if wantq('gaussian')
    savedata.gaussianData = gaussianData;
end

if wantq('levscore')
    savedata.levscoreData = levscoreData;
end

if wantq('froblev')
    savedata.froblevData = froblevData;
end

if wantq('speclev')
    savedata.speclevData = speclevData;
end

if wantq('approxlev')
    savedata.approxlevData = approxlevData;
end

if wantq('tallthinlevscoreapprox')
    savedata.tallthinData = tallthinData;
    savedata.qrlevscoreData = qrlevscoreData;
end

save(datasetfname, 'savedata');

end
