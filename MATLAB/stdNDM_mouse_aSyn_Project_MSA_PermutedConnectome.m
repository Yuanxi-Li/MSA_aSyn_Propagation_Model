function outputs_all = stdNDM_mouse_aSyn_Project_MSA_PermutedConnectome(varargin)
% STDNDM_MOUSE_ASYN_PROJECT_MSA_PERMUTEDCONNECTOME
% Fit the MSA-alpha-Syn model using element-wise permuted connectomes.
%
% In each permutation, all elements of the original connectome are randomly
% reassigned to new matrix locations. This preserves the complete weight
% distribution and the number of zero elements, while disrupting the
% original connectome organization.
%
% Parameters fitted:
%   param(1) - Initial seeding rescale factor, gamma.
%   param(2) - Global amplification/clearance rate, alpha.
%   param(3) - Global spreading rate, beta.
%   param(4) - Directionality parameter, s.
%
% Parameters retained for compatibility but not used:
%   param(5) - Source coefficient, a; fixed to zero.
%   param(6) - Transmission-modulation coefficient, b; fixed to zero.
%   param(7) - Regional amplification-modulation coefficient, p; fixed to zero.


%% Define defaults

study_ = 'aSyn_Project';
costfun_ = 'log_LinR_sum';
solvetype_ = 'analytic';

% Retained for compatibility but not used in the current project.
volcorrect_ = 0;

exclseed_costfun_ = 0;
excltpts_costfun_ = [];

% Retained for interface and output compatibility. Pathology is not
% transformed by this main function.
normtype_ = 'log';

w_dir_ = 0;

% Directionality starts from the neutral value s = 0.5.
param_init_ = [nan, 0, 1, 0.5];

ub_ = [0.000057, Inf, Inf, 1];

lb_ = zeros(1, 4);
lb_(2) = -Inf;

algo_ = 'sqp';

opttol_ = 1e-12;
fxntol_ = 1e-12;
steptol_ = 1e-12;
maxeval_ = 100000;

% Retained for interface compatibility. The current permutation analysis is
% implemented by the outer 1000-iteration loop.
bootstrapping_ = 0;
resample_rate_ = 0.8;
niters_ = 100;

verbose_ = 0;
fmindisplay_ = 0;

% Retained for interface compatibility but not used.
flowthresh_ = 99.93;


%% Parse optional inputs

ip = inputParser;

validScalar = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validBoolean = @(x) isscalar(x) && (x == 0 || x == 1);
validChar = @(x) ischar(x);
validStudy = @(x) ismember(x, {'aSyn_Project'});
validST = @(x) ismember(x, {'analytic', 'numeric'});
validParam = @(x) length(x) == 4;

addParameter(ip, 'study', study_, validStudy);
addParameter(ip, 'costfun', costfun_, validChar);
addParameter(ip, 'solvetype', solvetype_, validST);
addParameter(ip, 'volcorrect', volcorrect_, validBoolean);
addParameter(ip, 'exclseed_costfun', exclseed_costfun_, validBoolean);
addParameter(ip, 'excltpts_costfun', excltpts_costfun_);
addParameter(ip, 'normtype', normtype_, validChar);
addParameter(ip, 'w_dir', w_dir_, validBoolean);
addParameter(ip, 'param_init', param_init_, validParam);
addParameter(ip, 'ub', ub_, validParam);
addParameter(ip, 'lb', lb_, validParam);
addParameter(ip, 'opttol', opttol_, validScalar);
addParameter(ip, 'steptol', steptol_, validScalar);
addParameter(ip, 'fxntol', fxntol_, validScalar);
addParameter(ip, 'algo', algo_, validChar);
addParameter(ip, 'maxeval', maxeval_, validScalar);
addParameter(ip, 'bootstrapping', bootstrapping_, validBoolean);
addParameter(ip, 'resample_rate', resample_rate_, validScalar);
addParameter(ip, 'niters', niters_, validScalar);
addParameter(ip, 'verbose', verbose_, validBoolean);
addParameter(ip, 'fmindisplay', fmindisplay_, validBoolean);
addParameter(ip, 'flowthresh', flowthresh_, validScalar);

parse(ip, varargin{:});

ipR = ip.Results;

global C


%% Load pathology, seeding, and connectome data

data_dir = fullfile(getenv('HOME'), 'DataInput');

load(fullfile(data_dir, 'Pathology_Data_Input.mat'), ...
    'Pathology_Data', 'Seed_Data', 'tpts');

load(fullfile(data_dir, 'Connectome.mat'), ...
    'Connectome');

C_original = Connectome.raw;


%% Define fixed model inputs

% Only the early 0-3 MPI phase is modeled.
time_stamps = tpts.MSA_Average(1);

pathology = median( ...
    Pathology_Data.MSA_Individual_3_MPI, ...
    2);

seed_location = Seed_Data;


%% Configure fmincon

if logical(ipR.fmindisplay)

    options = optimoptions( ...
        @fmincon, ...
        'Display', 'final-detailed', ...
        'Algorithm', ipR.algo, ...
        'MaxFunctionEvaluations', ipR.maxeval, ...
        'OptimalityTolerance', ipR.opttol, ...
        'FunctionTolerance', ipR.fxntol, ...
        'StepTolerance', ipR.steptol);

else

    options = optimoptions( ...
        @fmincon, ...
        'Algorithm', ipR.algo, ...
        'MaxFunctionEvaluations', ipR.maxeval, ...
        'OptimalityTolerance', ipR.opttol, ...
        'FunctionTolerance', ipR.fxntol, ...
        'StepTolerance', ipR.steptol);

end


%% Perform element-wise connectome permutations

outputs_all = struct;

tic

for PermutationNum = 1:1000

    try

        fprintf( ...
            'Element-wise connectome permutation %d/1000\n', ...
            PermutationNum);


        %% Generate a reproducible permuted connectome

        rng(PermutationNum);
        rng_state = rng;

        random_order = randperm(numel(C_original));

        C = reshape( ...
            C_original(random_order), ...
            size(C_original));

        % Preserve the original maximum-eigenvalue normalization.
        C = C / max(eig(C));


        %% Define regional modulation input

        % Regional feature modulation is not used in the current model.
        U = zeros(size(C, 1), 1);


        %% Initialize model parameters

        if isnan(ipR.param_init(1))

            ipR.param_init(1) = ...
                nansum(pathology(:, 1)) / nnz(seed_location);

        end

        % The permuted connectome has no prior directional organization.
        % Directionality starts from s = 0.5 and remains freely optimized
        % over the complete interval [0,1].
        if ~logical(ipR.w_dir)

            ipR.param_init(4) = 0.5;
            ipR.ub(4) = 1;
            ipR.lb(4) = 0;

        end


        %% Count fitted parameters for information criteria

        mordervec = zeros(1, length(ipR.param_init));

        for m = 1:length(ipR.param_init)

            inclparam = ipR.lb(m) ~= ipR.ub(m);
            mordervec(m) = inclparam;

        end

        morder = 1 + sum(mordervec);


        %% Append fixed unused parameters

        param_init = [ipR.param_init, zeros(1, 3)];
        lb = [ipR.lb, zeros(1, 3)];
        ub = [ipR.ub, zeros(1, 3)];


        %% Define the optimization objective

        objfun_handle = @(param) ...
            objfun_eNDM_general_dir_costopts( ...
            param, ...
            seed_location, ...
            pathology, ...
            time_stamps, ...
            C, ...
            U, ...
            ipR.solvetype, ...
            ipR.volcorrect, ...
            ipR.costfun, ...
            ipR.excltpts_costfun, ...
            ipR.exclseed_costfun);


        %% Fit the model

        if logical(ipR.bootstrapping)

            error([ ...
                'Internal bootstrapping is not used in the current ', ...
                'connectome-permutation analysis.']);

        end

        [param_num, fval_num] = fmincon( ...
            objfun_handle, ...
            param_init, ...
            [], ...
            [], ...
            [], ...
            [], ...
            lb, ...
            ub, ...
            'objfun_eNDM_general_dir_nlcon_MSA', ...
            options);


        %% Generate predictions

        x0_num = seed_location * param_num(1);

        alpha_num = param_num(2);
        beta_num = param_num(3);
        s_num = param_num(4);

        % These parameters remain fixed at zero.
        a_num = param_num(5);
        b_num = param_num(6);
        p_num = param_num(7);

        ynum = eNDM_general_dir( ...
            x0_num, ...
            time_stamps, ...
            C, ...
            U, ...
            alpha_num, ...
            beta_num, ...
            s_num, ...
            a_num, ...
            b_num, ...
            p_num, ...
            ipR.solvetype, ...
            ipR.volcorrect);


        %% Store model inputs and fitted outputs

        outputs.ndm = struct;

        outputs.ndm.Full.data = pathology;
        outputs.ndm.Full.time_stamps = time_stamps;
        outputs.ndm.Full.predicted = ynum;
        outputs.ndm.Full.param_fit = param_num;
        outputs.ndm.Full.fval = fval_num;

        outputs.ndm.Full.init.study = ipR.study;
        outputs.ndm.Full.init.solvetype = ipR.solvetype;
        outputs.ndm.Full.init.volcorrect = ipR.volcorrect;
        outputs.ndm.Full.init.normtype = ipR.normtype;
        outputs.ndm.Full.init.costfun = ipR.costfun;
        outputs.ndm.Full.init.exclseed_costfun = ...
            ipR.exclseed_costfun;
        outputs.ndm.Full.init.excltpts_costfun = ...
            ipR.excltpts_costfun;
        outputs.ndm.Full.init.w_dir = ipR.w_dir;
        outputs.ndm.Full.init.param_init = ipR.param_init;
        outputs.ndm.Full.init.ub = ipR.ub;
        outputs.ndm.Full.init.lb = ipR.lb;
        outputs.ndm.Full.init.bootstrapping = ipR.bootstrapping;
        outputs.ndm.Full.init.resample_rate = ipR.resample_rate;
        outputs.ndm.Full.init.niters = ipR.niters;

        % Save the information required to reconstruct the permutation
        % without storing 1000 complete connectome matrices.
        outputs.ndm.Full.init.permutation_type = ...
            'element-wise';

        outputs.ndm.Full.init.permutation_number = ...
            PermutationNum;

        outputs.ndm.Full.init.random_seed = ...
            PermutationNum;

        outputs.ndm.Full.init.rng_state = ...
            rng_state;

        outputs.ndm.Full.init.connectome_size = ...
            size(C_original);

        outputs.ndm.Full.init.connectome_normalization = ...
            'maximum_eigenvalue';

        outputs.ndm.Full.fmincon.optimality_tolerance = ipR.opttol;
        outputs.ndm.Full.fmincon.function_tolerance = ipR.fxntol;
        outputs.ndm.Full.fmincon.step_tolerance = ipR.steptol;
        outputs.ndm.Full.fmincon.algorithm = ipR.algo;
        outputs.ndm.Full.fmincon.max_evaluations = ipR.maxeval;


        %% Calculate model-performance metrics

        Rvalues = zeros(1, length(time_stamps));
        LogRvalues = zeros(1, length(time_stamps));
        LinRvalues = zeros(1, length(time_stamps));
        LogLinRvalues = zeros(1, length(time_stamps));
        sse_individual = zeros(1, length(time_stamps));

        finiteinds = isfinite(sum(pathology, 2));
        Logfiniteinds = isfinite(sum(log(pathology), 2));

        LinRcalc = @(x, y) ...
            2 * corr(x, y) * std(x) * std(y) / ...
            (std(x)^2 + std(y)^2 + ...
            (mean(x) - mean(y))^2);

        for jj = 1:length(time_stamps)

            Rvalues(jj) = corr( ...
                ynum(:, jj), ...
                pathology(:, jj), ...
                'rows', 'complete');

            LogRvalues(jj) = corr( ...
                log(ynum(Logfiniteinds, jj)), ...
                log(pathology(Logfiniteinds, jj)), ...
                'rows', 'complete');

            LinRvalues(jj) = LinRcalc( ...
                ynum(finiteinds, jj), ...
                pathology(finiteinds, jj));

            LogLinRvalues(jj) = LinRcalc( ...
                log(ynum(Logfiniteinds, jj)), ...
                log(pathology(Logfiniteinds, jj)));

            sse_individual(jj) = ...
                sum( ...
                (ynum(:, jj) - pathology(:, jj)).^2, ...
                'all', ...
                'omitnan') / ...
                length( ...
                find( ...
                ~isnan( ...
                (ynum(:, jj) - pathology(:, jj)).^2) == 1));

        end

        NewRLog = corr( ...
            reshape(log(ynum(Logfiniteinds, :)), [], 1), ...
            reshape(log(pathology(Logfiniteinds, :)), [], 1), ...
            'rows', 'complete');


        %% Store model-performance metrics

        outputs.ndm.Full.results.data_means = ...
            mean(pathology, 1, 'omitnan');

        outputs.ndm.Full.results.Corrs = Rvalues;
        outputs.ndm.Full.results.Corrs_Mean = mean(Rvalues);

        outputs.ndm.Full.results.LogCorrs = LogRvalues;
        outputs.ndm.Full.results.LogCorrs_Mean = ...
            mean(LogRvalues);

        outputs.ndm.Full.results.LinR = LinRvalues;
        outputs.ndm.Full.results.LinR_Mean = mean(LinRvalues);

        outputs.ndm.Full.results.New_R_Log = NewRLog;

        outputs.ndm.Full.results.LogLinR = LogLinRvalues;
        outputs.ndm.Full.results.LogLinR_Mean = ...
            mean(LogLinRvalues);

        outputs.ndm.Full.results.sse_all = ...
            sum( ...
            (ynum - pathology).^2, ...
            'all', ...
            'omitnan') / ...
            length( ...
            find( ...
            ~isnan((ynum - pathology).^2) == 1));

        outputs.ndm.Full.results.sse_individual = ...
            sse_individual;


        %% Calculate regression-based model statistics

        P = reshape(pathology, [], 1);
        Y = reshape(ynum, [], 1);

        numObs1 = length(P(~isnan(P)));

        lm_endm = fitlm(Y, P);
        logL = lm_endm.LogLikelihood;

        outputs.ndm.Full.results.lm_LogL = logL;

        outputs.ndm.Full.results.lm_AIC = ...
            -2 * logL + 2 * morder;

        outputs.ndm.Full.results.lm_BIC = ...
            -2 * logL + log(numObs1) * morder;

        outputs.ndm.Full.results.lm_intercept = ...
            lm_endm.Coefficients.Estimate(1);

        outputs.ndm.Full.results.lm_pval = ...
            lm_endm.Coefficients.pValue(1);

        outputs.ndm.Full.results.lm_Rsquared_ord = ...
            lm_endm.Rsquared.Ordinary;

        outputs.ndm.Full.results.lm_Rsquared_adj = ...
            lm_endm.Rsquared.Adjusted;


        %% Display fitted results

        if logical(ipR.verbose)

            disp('--------------------------------------------------')
            disp(['Permutation = ' num2str(PermutationNum)])
            disp(['Optimal seed rescale value = ' ...
                num2str(param_num(1))])
            disp(['Optimal alpha = ' num2str(param_num(2))])
            disp(['Optimal beta = ' num2str(param_num(3))])
            disp(['Optimal s = ' num2str(param_num(4))])
            disp(['Log CCC = ' ...
                num2str(outputs.ndm.Full.results.LogLinR_Mean)])
            disp(['Cost = ' num2str(fval_num)])
            disp(' ')

        end


        %% Store the current permutation

        fldname = sprintf( ...
            'Permutation_%d', ...
            PermutationNum);

        outputs_all.(fldname) = outputs.ndm;

        disp([ ...
            fldname, ...
            ' Finished']);

    catch ME

        fprintf( ...
            'Permutation %d failed: %s\n', ...
            PermutationNum, ...
            ME.message);

        continue

    end

end

toc

end