function [c, ceq] = objfun_eNDM_general_dir_nlcon_gene_PFF(param)
% OBJFUN_ENDM_GENERAL_DIR_NLCON_GENE_PFF
% Nonlinear inequality constraint for the PFF gene-modulated model.
%
% The constraint requires predicted pathology at both PFF time points to
% remain less than or equal to 1.
%
% Parameter order:
%   param(1) - Initial seeding rescale factor.
%   param(2) - Amplification/clearance rate.
%   param(3) - Spreading rate.
%   param(4) - Directionality, s.
%   param(5) - Synaptic-locus parameter, lambda.
%   Remaining parameters are a, b, and p and are fixed to zero in the
%   current project.

data_dir = fullfile(getenv('HOME'), 'DataInput');

load(fullfile(data_dir, 'Pathology_Data_Input.mat'), ...
    'Seed_Data', 'tpts');

global C
global Gene_expression_diag

C_ = C;
Gene_expression_diag_ = Gene_expression_diag;

U_ = zeros(size(C_, 1), 1);
n_types = size(U_, 2);

x0_ = param(1) * Seed_Data;
alpha_ = param(2);
beta_ = param(3);
s_ = param(4);
outgoingpara_ = param(5);
a_ = param(6:(n_types + 5));
b_ = param((n_types + 6):(2 * n_types + 5));
p_ = param((2 * n_types + 6):(3 * n_types + 5));

ts = tpts.PFF_Average;
solvetype_ = 'analytic';
volcorrect_ = 0;

y = eNDM_general_dir_geneeffect( ...
    x0_, ts, C_, U_, alpha_, beta_, s_, outgoingpara_, ...
    a_, b_, p_, Gene_expression_diag_, solvetype_, volcorrect_);

% Vectorization applies the constraint to all ROIs at both time points.
c = y(:) - 1;
ceq = [];

end
