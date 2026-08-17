function [y, A] = eNDM_general_dir_geneeffect( ...
    x0, time_stamps, C, U, alpha, beta, s, outgoingpara, ...
    a, b, p, gene_exp_diag, solvetype, volcorrect)
% ENDM_GENERAL_DIR_GENEEFFECT Evaluate the directional gene-modulated eNDM.
%
% The gene-weighted connectome is defined as:
%
%   C_gene = G^lambda * C * G^(1-lambda)
%
% where G is the diagonal regional gene-expression matrix and lambda is
% represented by outgoingpara in this implementation.
%
% Inputs:
%   x0            - Initial seeding vector (n_ROI x 1).
%   time_stamps   - Experimental time points.
%   C             - Directed connectome (n_ROI x n_ROI).
%   U             - Regional feature matrix. Not used in the current
%                   project and supplied as zeros.
%   alpha         - Global amplification/clearance rate.
%   beta          - Global connectome-mediated spreading rate.
%   s             - Directionality parameter in [0,1].
%                   s = 0: fully anterograde.
%                   s = 1: fully retrograde.
%   outgoingpara  - Continuous synaptic-locus parameter lambda in [0,1].
%                   lambda = 1: presynaptic/outgoing modulation.
%                   lambda = 0: postsynaptic/incoming modulation.
%   a             - Source coefficients. Fixed to zero in this project.
%   b             - Transmission coefficients. Retained for compatibility
%                   but not applied in the current implementation.
%   p             - Regional amplification coefficients. Fixed to zero in
%                   this project.
%   gene_exp_diag - Diagonal regional gene-expression matrix.
%   solvetype     - 'analytic' or 'numeric'.
%   volcorrect    - Retained for compatibility; not used in this project.
%
% Outputs:
%   y - Predicted regional pathology at each time point.
%   A - Linear system matrix.

if nargin < 14
    volcorrect = 0;
    if nargin < 13
        solvetype = 'numeric';
    end
end

% volcorrect is retained for interface compatibility but is not applied in
% the current project.
%#ok<NASGU>

%% Source term

if size(a, 2) > size(a, 1)
    a = a.';
end

s_a = U * a;


%% Regional amplification/clearance matrix

if size(p, 2) > size(p, 1)
    p = p.';
end

s_p = U * p;
Gamma = diag(alpha + s_p);


%% Gene-weighted directional graph Laplacian

% The element-wise power formulation is retained from the original model
% implementation to preserve the existing fitting behavior.
C_Gene = (gene_exp_diag.^(outgoingpara)) * ...
    C * (gene_exp_diag.^(1 - outgoingpara));

C_Gene = C_Gene / max(eig(C_Gene));

C_dir = (1 - s) * C_Gene.' + s * C_Gene;
coldegree = sum(C_dir, 1);
L = diag(coldegree) - C_dir;

% b is retained in the function interface but is not used in the current
% project model.
%#ok<NASGU>


%% Define the linear system

A = Gamma - beta * L;
B = s_a;


%% Solve the system

if all(B == 0) && strcmp(solvetype, 'analytic')

    y = zeros(size(C, 1), length(time_stamps));
    y_analytic = @(t, initial_state) expm(A * t) * initial_state;

    for i = 1:length(time_stamps)
        y(:, i) = y_analytic(time_stamps(i), x0);
    end

else

    [~, y] = ode45(@(t, x) A * x + B, time_stamps, x0);
    y = y.';

end

end
