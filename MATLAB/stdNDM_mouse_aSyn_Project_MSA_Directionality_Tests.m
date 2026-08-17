function outputs_all = ...
    stdNDM_mouse_aSyn_Project_MSA_Directionality_Tests(varargin)
% STDNDM_MOUSE_ASYN_PROJECT_MSA_DIRECTIONALITY_TESTS
% Evaluate MSA-alpha-Syn model performance across fixed directionality
% values from fully anterograde to fully retrograde spreading.
%
% For each fixed directionality value, the model independently optimizes:
%   param(1) - Initial seeding rescale factor, gamma.
%   param(2) - Global amplification/clearance rate, alpha.
%   param(3) - Global spreading rate, beta.
%
% The directionality parameter is fixed separately at each value:
%   param(4) - Directionality, s.
%
% Parameters retained for compatibility but not used:
%   param(5) - Source coefficient, a; fixed to zero.
%   param(6) - Transmission-modulation coefficient, b; fixed to zero.
%   param(7) - Regional amplification-modulation coefficient, p; fixed to zero.
%
% Under the current convention:
%   s = 0 corresponds to fully anterograde spreading.
%   s = 1 corresponds to fully retrograde spreading.
%
% MSA-alpha-Syn pathology declined between 3 and 6 MPI. Therefore, this
% analysis focuses on the early 0-3 MPI phase, when pathology remained in
% the net accumulation stage.


%% Define defaults

study_ = 'aSyn_Project';
costfun_ = 'log_LinR_sum';
solvetype_ = 'analytic';

% Retained for compatibility but not used in the current project.
volcorrect_ = 0;

exclseed_costfun_ = 0;
excltpts_costfun_ = [];

% Retained for interface and output compatibility. Pathology is passed
% directly to the objective function.
normtype_ = 'log';

% Retained for compatibility. Directionality is explicitly fixed by the
% directionality scan below.
w_dir_ = 0;

param_init_ = [nan, 0, 1, 1];

% The PFF global model yielded an optimal initial seeding scale of
% approximately 5.7e-3. Based on the experimentally defined difference in
% injected material, the MSA-alpha-Syn upper bound is set 100-fold lower.
ub_ = [0.000057, Inf, Inf, 1];

lb_ = zeros(1, 4);
lb_(2) = -Inf;

algo_ = 'sqp';

opttol_ = 1e-12;
fxntol_ = 1e-12;
steptol_ = 1e-12;

maxeval_ = 100000;

% Retained for interface compatibility. Internal bootstrapping is not used
% in the current directionality analysis.
bootstrapping_ = 0;
resample_rate_ = 0.8;
niters_ = 100;

verbose_ = 0;
fmindisplay_ = 0;

% Retained for interface compatibility. Flow analysis is not used.
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


%% Define and normalize the connectome

C = Connectome.raw;

% Preserve the original normalization based on the maximum eigenvalue.
C = C / max(eig(C));


%% Evaluate fixed directionality values

outputs.ndm = struct;

% Direct_1 corresponds to s = 0, Direct_2 to s = 0.01, and
% Direct_101 to s = 1.
directionality = 0:0.01:1;

for directionality_n = 1:length(directionality)

    tic

    if ~logical(ipR.bootstrapping)

        fixed_s = directionality(directionality_n);

        fprintf( ...
            'Fitting MSA-alpha-Syn model at fixed s = %.2f\n', ...
            fixed_s);


        %% Select MSA-alpha-Syn pathology data

        % Only the 3 MPI pathology map is modeled.
        time_stamps = tpts.MSA_Average(1);

        pathology = median( ...
            Pathology_Data.MSA_Individual_3_MPI, ...
            2);

        seed_location = Seed_Data;


        %% Define unused regional modulation input

        % Cell-type and regional feature modulation are not used in the
        % current global model.
        U = zeros(size(C, 1), 1);


        %% Initialize the seeding scale

        if isnan(ipR.param_init(1))

            ipR.param_init(1) = ...
                nansum(pathology(:, 1)) / nnz(seed_location);

        end


        %% Fix the directionality parameter

        % For each model, s is fixed while gamma, alpha, and beta are
        % optimized.
        ipR.param_init(4) = fixed_s;
        ipR.ub(4) = fixed_s;
        ipR.lb(4) = fixed_s;


        %% Count fitted parameters for information criteria

        mordervec = zeros(1, length(ipR.param_init));

        for m = 1:length(ipR.param_init)

            inclparam = ipR.lb(m) ~= ipR.ub(m);
            mordervec(m) = inclparam;

        end

        % Include one additional parameter for the regression intercept.
        morder = 1 + sum(mordervec);


        %% Append fixed dummy parameters

        % a, b, and p are not used in the current project and remain fixed
        % at zero.
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


        %% Fit the model

        % Use the MSA-specific nonlinear constraint at the 3 MPI time point.
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


        %% Extract fitted parameters

        x0_num = seed_location * param_num(1);

        alpha_num = param_num(2);
        beta_num = param_num(3);
        s_num = param_num(4);

        % These parameters are fixed at zero and are not fitted.
        a_num = param_num(5);
        b_num = param_num(6);
        p_num = param_num(7);


        %% Generate model predictions

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

        outputs.ndm.Full.data = pathology;
        outputs.ndm.Full.time_stamps = time_stamps;
        outputs.ndm.Full.predicted = ynum;
        outputs.ndm.Full.param_fit = param_num;
        outputs.ndm.Full.fval = fval_num;

        outputs.ndm.Full.init.C = C;
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

            squared_error = ...
                (ynum(:, jj) - pathology(:, jj)).^2;

            sse_individual(jj) = ...
                sum(squared_error, 'all', 'omitnan') / ...
                length(find(~isnan(squared_error) == 1));

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

        squared_error_all = (ynum - pathology).^2;

        outputs.ndm.Full.results.sse_all = ...
            sum(squared_error_all, 'all', 'omitnan') / ...
            length(find(~isnan(squared_error_all) == 1));

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
            disp(['Fixed directionality s = ' num2str(fixed_s)])
            disp(' ')

            disp(['Optimal seed rescale value = ' ...
                num2str(param_num(1))])

            disp(['Optimal alpha = ' ...
                num2str(param_num(2))])

            disp(['Optimal beta = ' ...
                num2str(param_num(3))])

            disp(['Fixed s = ' ...
                num2str(param_num(4))])

            disp(' ')
            disp('CCC values at each time stamp')
            disp(LinRvalues)
            disp(' ')

            disp(ipR.costfun)
            disp(fval_num)

        end

    else

        error([ ...
            'Internal bootstrapping is not used in the current ', ...
            'directionality analysis. Set bootstrapping to 0.']);

    end


    %% Store the current fixed-directionality model

    Name_Result = append( ...
        'Direct_', ...
        num2str(directionality_n));

    disp([ ...
        Name_Result, ...
        ': fixed s = ', ...
        num2str(fixed_s), ...
        ', CCC = ', ...
        num2str(outputs.ndm.Full.results.LogLinR_Mean)]);

    outputs_all.(Name_Result) = outputs.ndm.Full;

    toc

end

end