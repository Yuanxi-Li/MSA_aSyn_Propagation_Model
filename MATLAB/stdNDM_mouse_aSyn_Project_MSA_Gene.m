function outputs_all = stdNDM_mouse_aSyn_Project_MSA_Gene(varargin)
% STDNDM_MOUSE_ASYN_PROJECT_MSA_GENE
% Fit the synaptic-locus-specific gene-modulated eNDM to MSA-alpha-Syn data.
%
% For each gene, regional expression is mapped to the 410 modeled regions
% and incorporated into the connectome as:
%
%   C_gene = G^lambda * C * G^(1-lambda)
%
% Parameters fitted:
%   param(1) - Initial seeding rescale factor, gamma.
%   param(2) - Global amplification/clearance rate, alpha.
%   param(3) - Global spreading rate, beta.
%   param(4) - Directionality, s.
%   param(5) - Synaptic-locus parameter, lambda.
%
% Parameters retained for compatibility but not used:
%   param(6) - Source coefficient, a; fixed to zero.
%   param(7) - Transmission coefficient, b; fixed to zero.
%   param(8) - Regional amplification coefficient, p; fixed to zero.
%
% MSA-alpha-Syn is modeled at 3 MPI only. Directionality is fixed at s = 0
% in this gene-modulated analysis, while lambda is fitted within [0,1].


%% Define defaults

study_ = 'aSyn_Project';
costfun_ = 'log_LinR_sum';
solvetype_ = 'analytic';

% Retained for compatibility; volume correction is not used.
volcorrect_ = 0;

exclseed_costfun_ = 0;
excltpts_costfun_ = [];
normtype_ = 'log';
w_dir_ = 0;

% Parameter order: gamma, alpha, beta, s, lambda.
param_init_ = [NaN, 0, 1, 1, 0.5];

% Preserve the original MSA-alpha-Syn parameter bounds.
ub_ = [0.000057, Inf, Inf, 1, 1];
lb_ = zeros(1, 5);
lb_(2) = -Inf;

algo_ = 'sqp';

% Preserve the original MSA-alpha-Syn optimization settings.
opttol_ = 1e-6;
fxntol_ = 1e-6;
steptol_ = 1e-12;
maxeval_ = 10000;

% Retained for interface compatibility. Internal bootstrapping is not used.
bootstrapping_ = 0;
resample_rate_ = 0.8;
niters_ = 100;

verbose_ = 0;
fmindisplay_ = 0;

% Retained for interface compatibility; flow analysis is not used.
flowthresh_ = 99.93;


%% Parse optional inputs

ip = inputParser;

validScalar = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validBoolean = @(x) isscalar(x) && (x == 0 || x == 1);
validChar = @(x) ischar(x);
validStudy = @(x) ismember(x, {'aSyn_Project'});
validST = @(x) ismember(x, {'analytic', 'numeric'});
validParam = @(x) length(x) == 5;

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
global Gene_expression_diag


%% Load pathology, gene-expression, seeding, and connectome data

data_dir = fullfile(getenv('HOME'), 'DataInput');

load(fullfile(data_dir, 'Pathology_Data_Input.mat'), ...
    'Pathology_Data', 'Seed_Data', 'tpts', ...
    'regvgene_mean', 'BrainRegionReorderMat', 'gene_names_trans');

load(fullfile(data_dir, 'Connectome.mat'), ...
    'Connectome');


%% Fit the gene-modulated model for each gene

outputs_all = struct;
tic

for gene_num = 1:length(gene_names_trans)

    try

        C = Connectome.raw;


        %% Map gene expression from the original 426-region atlas to 410 model regions

        % Gene-expression data are defined in the original 426-region space
        % (213 regions per hemisphere). The current pathology and connectome model
        % uses 410 regions (205 per hemisphere) after selected regions were merged
        % or excluded.
        %
        % BrainRegionReorderMat defines how the original unilateral 213-region
        % atlas is mapped to the unilateral 205-region model space. Each row
        % corresponds to one model region and may contain one or more original
        % region indices. When multiple original regions map to the same model
        % region, their gene-expression values are averaged.
        %
        % The bilateral mapping is constructed by appending the corresponding
        % right-hemisphere indices, offset by 213.
        Gene_expression_raw = regvgene_mean(:, gene_num);

        % Extend the original hemisphere-specific mapping to both
        % hemispheres while preserving the original region order.
        BrainRegionReorderMat_New = [ ...
            BrainRegionReorderMat; ...
            BrainRegionReorderMat + 213];

        [row_reorder, col_reorder] = ...
            size(BrainRegionReorderMat_New);

        Gene_expression_Mat = ...
            nan(size(BrainRegionReorderMat_New));

        for i = 1:row_reorder
            for j = 1:col_reorder
                if isnan(BrainRegionReorderMat_New(i, j))
                    Gene_expression_Mat(i, j) = NaN;
                else
                    Gene_expression_Mat(i, j) = ...
                        Gene_expression_raw( ...
                        BrainRegionReorderMat_New(i, j));
                end
            end
        end

        Gene_expression = ...
            mean(Gene_expression_Mat, 2, 'omitnan');

        % Replace regions without mapped expression using the whole-map
        % mean, preserving the original imputation procedure.
        Gene_expression(isnan(Gene_expression)) = ...
            mean(Gene_expression, 'omitnan');

        Gene_expression_diag = diag(Gene_expression);


        %% Solve and store results

        outputs.ndm = struct;

        if ~logical(ipR.bootstrapping)

            fprintf( ...
                'Fitting MSA-alpha-Syn gene %d/%d\n', ...
                gene_num, length(gene_names_trans));

            % Only the early 0-3 MPI phase is modeled for MSA-alpha-Syn.
            time_stamps = tpts.MSA_Average(1);


            pathology = median(Pathology_Data.MSA_Individual_3_MPI, 2);

            seed_location = Seed_Data;
            U = zeros(size(C, 1), 1);


            %% Initialize the seeding scale

            if isnan(ipR.param_init(1))
                ipR.param_init(1) = ...
                    nansum(pathology(:, 1)) / nnz(seed_location);
            end


            %% Fix MSA-alpha-Syn directionality at s = 0

            % Preserve the original MSA gene-model setting. Lambda remains
            % freely fitted within its original [0,1] bounds.
            if ~logical(ipR.w_dir)
                ipR.param_init(4) = 0;
                ipR.ub(4) = 0;
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

            % a, b, and p are not used and remain fixed at zero.
            param_init = [ipR.param_init, zeros(1, 3)];
            lb = [ipR.lb, zeros(1, 3)];
            ub = [ipR.ub, zeros(1, 3)];


            %% Define the optimization objective

            objfun_handle = @(param) ...
                objfun_eNDM_general_dir_geneeffect_costopts( ...
                param, seed_location, pathology, time_stamps, C, U, ...
                ipR.solvetype, ipR.volcorrect, ipR.costfun, ...
                ipR.excltpts_costfun, ipR.exclseed_costfun, ...
                Gene_expression_diag);


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

            [param_num, fval_num] = fmincon( ...
                objfun_handle, param_init, [], [], [], [], lb, ub, ...
                'objfun_eNDM_general_dir_nlcon_gene_MSA', options);


            %% Generate predictions using the fitted parameters

            x0_num = seed_location * param_num(1);
            alpha_num = param_num(2);
            beta_num = param_num(3);
            s_num = param_num(4);
            outgoingpara_num = param_num(5);

            % These parameters are fixed at zero.
            a_num = param_num(6);
            b_num = param_num(7);
            p_num = param_num(8);

            ynum = eNDM_general_dir_geneeffect( ...
                x0_num, time_stamps, C, U, alpha_num, beta_num, ...
                s_num, outgoingpara_num, a_num, b_num, p_num, ...
                Gene_expression_diag, ipR.solvetype, ipR.volcorrect);


            %% Store model inputs, gene information, and fitted outputs

            outputs.ndm.Full.data = pathology;
            outputs.ndm.Full.time_stamps = time_stamps;
            outputs.ndm.Full.predicted = ynum;
            outputs.ndm.Full.param_fit = param_num;
            outputs.ndm.Full.fval = fval_num;

            outputs.ndm.Full.init.C = C;
            outputs.ndm.Full.init.study = ipR.study;
            outputs.ndm.Full.init.gene_index = gene_num;
            outputs.ndm.Full.init.gene_name = ...
                gene_names_trans(gene_num);
            outputs.ndm.Full.init.gene_expression = Gene_expression;
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

            outputs.ndm.Full.fmincon.optimality_tolerance = ...
                ipR.opttol;
            outputs.ndm.Full.fmincon.function_tolerance = ...
                ipR.fxntol;
            outputs.ndm.Full.fmincon.step_tolerance = ...
                ipR.steptol;
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
                    ynum(:, jj), pathology(:, jj), ...
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




            %% Store model-performance metrics

            outputs.ndm.Full.results.data_means = ...
                mean(pathology, 1, 'omitnan');

            outputs.ndm.Full.results.Corrs = Rvalues;
            outputs.ndm.Full.results.Corrs_Mean = mean(Rvalues);
            outputs.ndm.Full.results.LogCorrs = LogRvalues;
            outputs.ndm.Full.results.LogCorrs_Mean = mean(LogRvalues);
            outputs.ndm.Full.results.LinR = LinRvalues;
            outputs.ndm.Full.results.LinR_Mean = mean(LinRvalues);
            outputs.ndm.Full.results.LogLinR = LogLinRvalues;
            outputs.ndm.Full.results.LogLinR_Mean = ...
                mean(LogLinRvalues);
            outputs.ndm.Full.results.synaptic_locus_lambda = ...
                outgoingpara_num;

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
                disp(['Gene index = ' num2str(gene_num)])
                disp(['Optimal seed rescale value = ' ...
                    num2str(param_num(1))])
                disp(['Optimal alpha = ' num2str(param_num(2))])
                disp(['Optimal beta = ' num2str(param_num(3))])
                disp(['Directionality s = ' num2str(param_num(4))])
                disp(['Synaptic-locus lambda = ' ...
                    num2str(param_num(5))])
                disp(['Mean log-CCC = ' ...
                    num2str(mean(LogLinRvalues))])
                disp(['Cost = ' num2str(fval_num)])
                disp(' ')
            end

        else

            error([ ...
                'Internal bootstrapping is not used in the current ', ...
                'gene-modulated analysis. Set bootstrapping to 0.']);

        end


        %% Store the current gene result

        fldname = sprintf('Gene_%d', gene_num);
        outputs_all.(fldname) = outputs.ndm;

        disp([fldname ' Finished'])

    catch ME

        fprintf('Gene %d failed: %s\n', gene_num, ME.message);
        continue

    end

end

toc


%% Data normalization

function normdata = normalizer(data, ntype)

    if strcmp(ntype, 'sum')
        normdata = data / nansum(data(:, 1));
    elseif strcmp(ntype, 'mean')
        normdata = data / nanmean(data(:, 1));
    elseif strcmp(ntype, 'norm2')
        normdata = data / norm(data(:, 1), 2);
    elseif strcmp(ntype, 'log')
        normdata = log(data + 1);
    else
        normdata = data;
    end

end

end
