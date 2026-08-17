function [f, newxt, newpath] = objfun_eNDM_general_dir_costopts( ...
    param, seed_location, pathology, ts, C_, U_, solvetype_, ...
    volcorrect_, costfun_, excltpts_costfun_, exclseed_costfun_)
% OBJFUN_ENDM_GENERAL_DIR_COSTOPTS
% Objective function for the directional eNDM model.
%
% Parameter order:
%   param(1)                         = Initial seeding rescale factor
%   param(2)                         = Global amplification/clearance rate
%   param(3)                         = Global spreading rate
%   param(4)                         = Directionality
%   param(5:(n_types+4))             = Source-term coefficients
%                                      Not used in the current project
%   param((n_types+5):(2*n_types+4)) = Transmission-modulation coefficients
%                                      Not used in the current project
%   param((2*n_types+5):(3*n_types+4)) = Regional amplification coefficients
%                                        Not used in the current project
%
% Current PFF fitting uses:
%   costfun_ = 'log_LinR_sum'
%
% Other cost-function options are retained for compatibility but are not
% used in the primary PFF analysis.

LinRcalc = @(x, y) ...
    2 * corr(x, y) * std(x) * std(y) / ...
    (std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

n_types = size(U_, 2);


%% Extract model parameters

x0_ = param(1) * seed_location;

alpha_ = param(2);
beta_ = param(3);
s_ = param(4);

a_ = param(5:(n_types + 4));
b_ = param((n_types + 5):(2 * n_types + 4));
p_ = param((2 * n_types + 5):(3 * n_types + 4));


%% Remove excluded time points

% No time points are excluded under the current default PFF settings.
ts(excltpts_costfun_) = [];
pathology(:, excltpts_costfun_) = [];


%% Calculate eNDM predictions

y = eNDM_general_dir( ...
    x0_, ...
    ts, ...
    C_, ...
    U_, ...
    alpha_, ...
    beta_, ...
    s_, ...
    a_, ...
    b_, ...
    p_, ...
    solvetype_, ...
    volcorrect_);


%% Optionally exclude the injection region

% Seed-region exclusion is not used under the current default settings.
if logical(exclseed_costfun_)

    seedbin = logical(seed_location);

    y(seedbin, :) = NaN;
    pathology(seedbin, :) = NaN;

end


%% Calculate the selected objective function

if strcmp(costfun_, 'sse_sum')

    f = nansum(nansum((y - pathology).^2));


elseif strcmp(costfun_, 'rval_sum')

    Rvalues = zeros(1, length(ts));

    for jj = 1:length(ts)

        Rvalues(jj) = corr( ...
            y(:, jj), ...
            pathology(:, jj), ...
            'rows', 'complete');

    end

    f = length(ts) - sum(Rvalues);


elseif strcmp(costfun_, 'sse_end')

    f = nansum(nansum( ...
        (y(:, end) - pathology(:, end)).^2));


elseif strcmp(costfun_, 'rval_end')

    Rvalues = zeros(1, length(ts));

    for jj = 1:length(ts)

        Rvalues(jj) = corr( ...
            y(:, jj), ...
            pathology(:, jj), ...
            'rows', 'complete');

    end

    f = 1 - Rvalues(end);


elseif strcmp(costfun_, 'LinR')

    Rvalues = zeros(1, length(ts));

    % Remove regions with missing pathology at any fitted time point.
    naninds = isnan(prod(pathology, 2));

    newxt = y;
    newxt(naninds, :) = [];

    newpath = pathology;
    newpath(naninds, :) = [];

    for jj = 1:length(ts)

        Rvalues(jj) = LinRcalc( ...
            newxt(:, jj), ...
            newpath(:, jj));

    end

    f = length(ts) - sum(Rvalues);


elseif strcmp(costfun_, 'LinR_end')

    Rvalues = zeros(1, length(ts));

    % Remove regions with missing pathology at any fitted time point.
    naninds = isnan(prod(pathology, 2));

    newxt = y;
    newxt(naninds, :) = [];

    newpath = pathology;
    newpath(naninds, :) = [];

    for jj = 1:length(ts)

        Rvalues(jj) = LinRcalc( ...
            newxt(:, jj), ...
            newpath(:, jj));

    end

    f = 1 - Rvalues(end);


elseif strcmp(costfun_, 'log_rval_sum')

    Rvalues = zeros(1, length(ts));

    % Retain regions with finite log-transformed pathology values at all
    % fitted time points.
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

    % Primary objective function used for the current PFF analysis.
    Rvalues = zeros(1, length(ts));

    % A region is included only when its observed pathology is positive and
    % finite at every fitted time point.
    finiteinds = isfinite(sum(log(pathology), 2));

    newxt = y;

    for jj = 1:length(ts)

        Rvalues(jj) = LinRcalc( ...
            log(y(finiteinds, jj)), ...
            log(pathology(finiteinds, jj)));

    end

    f = length(ts) - sum(Rvalues);

end

end