function [y, A] = eNDM_general_dir( ...
    x0, time_stamps, C, U, alpha, beta, s, a, b, p, solvetype, volcorrect)
% eNDM_GENERAL_DIR Evaluate directional eNDM predictions.
%
% The model separates local pathology amplification/clearance from
% connectome-mediated pathology spreading:
%
%       dx/dt = A*x + B
%
% where:
%       A = Gamma - beta*L
%       B = U*a
%
% Inputs:
%   x0          - Initial pathology/seeding vector (n_ROI x 1).
%   time_stamps - Experimental time points.
%   C           - Directed connectivity matrix (n_ROI x n_ROI).
%                 C(i,j) represents connectivity from ROI i to ROI j.
%   U           - Regional feature or cell-type matrix.
%                 Not used in the current project; supplied as zeros.
%   alpha       - Global amplification/clearance rate.
%   beta        - Global connectome-mediated spreading rate.
%   s           - Directionality parameter in [0,1].
%                 s = 0: fully anterograde under the current convention.
%                 s = 1: fully retrograde under the current convention.
%   a           - Source-term coefficients.
%                 Not used in the current project; fixed to zero.
%   b           - Transmission-modulation coefficients.
%                 Not used in the current project; fixed to zero.
%   p           - Regional amplification-modulation coefficients.
%                 Not used in the current project; fixed to zero.
%   solvetype   - 'analytic' or 'numeric'.
%   volcorrect  - Retained for compatibility.
%                 Volume correction is not used in the current project.
%
% Outputs:
%   y - Predicted regional pathology at each time point.
%       Each column corresponds to one time point.
%   A - Linear system matrix.

if nargin < 12
    volcorrect = 0;

    if nargin < 11
        solvetype = 'numeric';
    end
end

% volcorrect is retained as an input for compatibility but is not currently
% applied to the Laplacian.

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


%% Directionally weighted graph Laplacian

C_dir = (1 - s) * C.' + s * C;

coldegree = sum(C_dir, 1);
L_raw = diag(coldegree) - C_dir;


%% Transmission modulation

% This block is retained for compatibility with the general eNDM function.
% b-based regional transmission modulation is not applied in the current
% project because b is fixed to zero and L is set directly to L_raw.

if size(b, 2) > size(b, 1)
    b = b.';
end

s_b = U * b;
S_b = repmat(s_b, 1, length(s_b)) + ones(length(s_b)); 

L = L_raw;


%% Define the system

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