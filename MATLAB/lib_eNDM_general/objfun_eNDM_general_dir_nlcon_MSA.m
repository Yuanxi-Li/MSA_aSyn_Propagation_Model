function [c, ceq] = objfun_eNDM_general_dir_nlcon_MSA(param)
% OBJFUN_ENDM_GENERAL_DIR_NLCON_MSA
% Define nonlinear constraints for fitting the directional eNDM to MSA.
%
% The constraint requires all predicted pathology values at the modeled
% MSA time point to remain no greater than 1.
%
% Parameter order:
%   param(1) - Initial seeding rescale factor
%   param(2) - Global amplification/clearance rate, alpha
%   param(3) - Global spreading rate, beta
%   param(4) - Directionality parameter, s
%   param(5) - Source-term coefficient, a
%              Not used in the current project; fixed to zero
%   param(6) - Transmission-modulation coefficient, b
%              Not used in the current project; fixed to zero
%   param(7) - Regional amplification-modulation coefficient, p
%              Not used in the current project; fixed to zero


%% Load seeding data and experimental time points

data_dir = fullfile(getenv('HOME'), 'DataInput');

load(fullfile(data_dir, 'Pathology_Data_Input.mat'), ...
    'Seed_Data', 'tpts');


%% Retrieve the connectome used by the fitting function

global C

C_ = C;


%% Define regional modulation inputs

% Cell-type or regional-feature modulation is not used in the current
% project. One zero-valued feature is retained to preserve the general
% eNDM parameter structure.
U_ = zeros(size(C_, 1), 1);

n_types = size(U_, 2);


%% Extract model parameters

x0_ = param(1) * Seed_Data;

alpha_ = param(2);
beta_ = param(3);
s_ = param(4);

a_ = param(5:(n_types + 4));
b_ = param((n_types + 5):(2 * n_types + 4));
p_ = param((2 * n_types + 5):(3 * n_types + 4));


%% Define the MSA time point

% MSA fitting is restricted to the early 0-3 MPI propagation phase.
ts = tpts.MSA_Average(:, 1);

solvetype_ = 'analytic';

% Volume correction is not used in the current project.
volcorrect_ = 0;


%% Calculate predicted pathology

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


%% Define nonlinear constraints

% fmincon requires c <= 0. Vectorization ensures that the constraint is
% returned as a column vector.
c = y(:) - 1;

% No nonlinear equality constraints are used.
ceq = [];

end