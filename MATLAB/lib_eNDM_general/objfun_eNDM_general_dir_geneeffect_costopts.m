function [f, newxt, newpath] = ...
    objfun_eNDM_general_dir_geneeffect_costopts( ...
    param, seed_location, pathology, ts, C_, U_, solvetype_, ...
    volcorrect_, costfun_, excltpts_costfun_, ...
    exclseed_costfun_, gene_exp_)
% OBJFUN_ENDM_GENERAL_DIR_GENEEFFECT_COSTOPTS
% Objective function for the directional gene-modulated eNDM.
%
% Parameter order:
%   param(1)                             - Initial seeding rescale factor.
%   param(2)                             - Amplification/clearance rate.
%   param(3)                             - Spreading rate.
%   param(4)                             - Directionality, s.
%   param(5)                             - Synaptic-locus parameter, lambda.
%   param(6:(n_types+5))                 - Source coefficients; not used.
%   param((n_types+6):(2*n_types+5))     - Transmission coefficients; not used.
%   param((2*n_types+6):(3*n_types+5))   - Amplification coefficients; not used.
%
% The current project primarily uses costfun_ = 'log_LinR_sum'. Other
% options are retained for compatibility with the original model code.

LinRcalc = @(x, y) ...
    2 * corr(x, y) * std(x) * std(y) / ...
    (std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

n_types = size(U_, 2);


%% Extract model parameters

x0_ = param(1) * seed_location;
alpha_ = param(2);
beta_ = param(3);
s_ = param(4);
outgoingpara_ = param(5);
a_ = param(6:(n_types + 5));
b_ = param((n_types + 6):(2 * n_types + 5));
p_ = param((2 * n_types + 6):(3 * n_types + 5));


%% Remove excluded time points

ts(excltpts_costfun_) = [];
pathology(:, excltpts_costfun_) = [];


%% Calculate gene-modulated eNDM predictions

y = eNDM_general_dir_geneeffect( ...
    x0_, ts, C_, U_, alpha_, beta_, s_, outgoingpara_, ...
    a_, b_, p_, gene_exp_, solvetype_, volcorrect_);

newxt = y;
newpath = pathology;


%% Optionally exclude the injection region

if logical(exclseed_costfun_)
    seedbin = logical(seed_location);
    y(seedbin, :) = NaN;
    pathology(seedbin, :) = NaN;
    newxt = y;
    newpath = pathology;
end


%% Calculate the selected objective function

if strcmp(costfun_, 'sse_sum')

    f = nansum(nansum((y - pathology).^2));

elseif strcmp(costfun_, 'rval_sum')

    Rvalues = zeros(1, length(ts));
    for jj = 1:length(ts)
        Rvalues(jj) = corr( ...
            y(:, jj), pathology(:, jj), 'rows', 'complete');
    end
    f = length(ts) - sum(Rvalues);

elseif strcmp(costfun_, 'sse_end')

    f = nansum(nansum((y(:, end) - pathology(:, end)).^2));

elseif strcmp(costfun_, 'rval_end')

    Rvalues = zeros(1, length(ts));
    for jj = 1:length(ts)
        Rvalues(jj) = corr( ...
            y(:, jj), pathology(:, jj), 'rows', 'complete');
    end
    f = 1 - Rvalues(end);

elseif strcmp(costfun_, 'LinR')

    Rvalues = zeros(1, length(ts));
    naninds = isnan(prod(pathology, 2));
    newxt = y;
    newxt(naninds, :) = [];
    newpath = pathology;
    newpath(naninds, :) = [];

    for jj = 1:length(ts)
        Rvalues(jj) = LinRcalc(newxt(:, jj), newpath(:, jj));
    end

    f = length(ts) - sum(Rvalues);

elseif strcmp(costfun_, 'LinR_end')

    Rvalues = zeros(1, length(ts));
    naninds = isnan(prod(pathology, 2));
    newxt = y;
    newxt(naninds, :) = [];
    newpath = pathology;
    newpath(naninds, :) = [];

    for jj = 1:length(ts)
        Rvalues(jj) = LinRcalc(newxt(:, jj), newpath(:, jj));
    end

    f = 1 - Rvalues(end);

elseif strcmp(costfun_, 'log_rval_sum')

    Rvalues = zeros(1, length(ts));
    finiteinds = isfinite(sum(log(pathology), 2));
    newxt = y;

    for jj = 1:length(ts)
        Rvalues(jj) = corr( ...
            log(y(finiteinds, jj)), ...
            log(pathology(finiteinds, jj)), ...
            'rows', 'complete');
    end

    f = length(ts) - sum(Rvalues);

elseif strcmp(costfun_, 'New_log_rval_sum')

    finiteinds = isfinite(sum(log(pathology), 2));
    newxt = y;

    Rvalues = corr( ...
        reshape(log(y(finiteinds, :)), [], 1), ...
        reshape(log(pathology(finiteinds, :)), [], 1), ...
        'rows', 'complete');

    f = 1 - Rvalues;

elseif strcmp(costfun_, 'log_rval_end')

    Rvalues = zeros(1, length(ts));
    finiteinds = isfinite(sum(log(pathology), 2));
    newxt = y;

    for jj = 1:length(ts)
        Rvalues(jj) = corr( ...
            log(y(finiteinds, jj)), ...
            log(pathology(finiteinds, jj)), ...
            'rows', 'complete');
    end

    f = 1 - Rvalues(end);

elseif strcmp(costfun_, 'log_LinR_sum')

    Rvalues = zeros(1, length(ts));
    finiteinds = isfinite(sum(log(pathology), 2));
    newxt = y;

    for jj = 1:length(ts)
        Rvalues(jj) = LinRcalc( ...
            log(y(finiteinds, jj)), ...
            log(pathology(finiteinds, jj)));
    end

    f = length(ts) - sum(Rvalues);

else

    error('Unknown cost function: %s', costfun_);

end

end
