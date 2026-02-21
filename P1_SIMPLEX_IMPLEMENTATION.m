% .-------------------------------------------------.
% |                                                 |
% |     __  __      _ ___              __           |
% |    / / / /___  (_)   |  ____  ____/ /__  _____  |
% |   / / / / __ \/ / /| | / __ \/ __  / _ \/ ___/  |
% |  / /_/ / / / / / ___ |/ / / / /_/ /  __(__  )   |
% |  \____/_/ /_/_/_/  |_/_/ /_/\__,_/\___/____/    |
% |                                                 |
% '-------------------------------------------------'
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ~~~~~~~~~~~~~~~~~~ POINT 1 ~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    Author  : Esteban Rodríguez
%    Date    : 02.2026
%    Purpose : Implementation of the Primal Simplex
%              Method using a finite-state structure.
%              The algorithm:
%              1. Checks optimality
%              2. Selects entering variable
%              3. Performs ratio test
%              4. Updates basis
%              5. Stores iteration trace
%    Notes   : Includes graphical visualization of
%              feasible region and simplex path.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% '-------------------------------------------------'

clear all; clc; clf;

% --------------------------------------------------
% PROBLEM DATA
% --------------------------------------------------
% Minimize f'x
% Subject to Ax <= b , x >= 0

A = [1 -1 1 0;
     1  1 0 1];

b = [2;6];

f = [-2 -1 0 0];

[m_A, n_A] = size(A);

% --------------------------------------------------
% INITIALIZATION
% --------------------------------------------------
basis_idx = [3 4];              % Initial slack-variable basis
STATE = "OPTIMALITY_CHECK";     % Finite state machine
iter = 0;
flag = 0;

% Trace storage
x_trace = zeros(m_A,1);         % Basic variables values history
basis_idx_trace = zeros(m_A,1); % Basis history

fprintf('SIMPLEX METHOD by Estrodz...\n');

% --------------------------------------------------
% SIMPLEX ALGORITHM
% --------------------------------------------------
while STATE ~= "TERMINATE" && flag ~= 1

    switch STATE

        % ------------------------------------------
        % 1. OPTIMALITY CHECK
        % ------------------------------------------
        case "OPTIMALITY_CHECK"

            iter = iter + 1;

            % Identify non-basic indices
            non_basis_idx = setdiff(1:n_A, basis_idx);

            % Construct basis and non-basis matrices
            B = A(:, basis_idx);
            N = A(:, non_basis_idx);

            cB = f(basis_idx);
            cN = f(non_basis_idx);

            % Compute dual variables and reduced costs
            y = cB / B;
            r_n = cN - y * N;

            fprintf('Iter: %d | Basis: %s | Reduced Costs: %s\n', ...
                iter, mat2str(basis_idx), mat2str(r_n));

            % Check optimality condition
            if all(r_n >= 0)
                STATE = "OPTIMAL_FOUND";
            else
                STATE = "SELECT_ENTERING";
            end

        % ------------------------------------------
        % 2. SELECT ENTERING VARIABLE
        % ------------------------------------------
        case "SELECT_ENTERING"

            [~, min_idx] = min(r_n);
            entering_var = non_basis_idx(min_idx);

            x_e = inv(B)*A(:, entering_var);

            if all(x_e <= 0)
                error('Unbounded problem detected.')
            else
                STATE = "RATIO_TEST";
            end

        % ------------------------------------------
        % 3. RATIO TEST AND BASIS UPDATE
        % ------------------------------------------
        case "RATIO_TEST"

            d = B \ A(:, entering_var);   % Direction
            xb = B \ b;                  % Current BFS

            % Store iteration trace
            x_trace = [x_trace, xb];
            basis_idx_trace = [basis_idx_trace, basis_idx.'];

            if all(d <= 0)
                STATE = "UNBOUNDED";
            else
                ratios = xb ./ d;
                ratios(d <= 0) = inf;
                [~, leaving_pos] = min(ratios);

                basis_idx(leaving_pos) = entering_var;
                STATE = "OPTIMALITY_CHECK";
            end

        % ------------------------------------------
        % 4. OPTIMAL SOLUTION FOUND
        % ------------------------------------------
        case "OPTIMAL_FOUND"

            x_final = zeros(n_A,1);
            x_final(basis_idx) = B \ b;

            x_trace = [x_trace, x_final(basis_idx)];
            basis_idx_trace = [basis_idx_trace, basis_idx.'];

            fprintf('\n--- OPTIMAL SOLUTION FOUND ---\n');
            fprintf('Optimal x = %s\n', mat2str(x_final));
            fprintf('Optimal objective value = %f\n\n', f*x_final);

            flag = 1;
            STATE = "TERMINATE";

        % ------------------------------------------
        % 5. UNBOUNDED CASE
        % ------------------------------------------
        case "UNBOUNDED"

            fprintf('\n--- ERROR ---\n');
            disp('The problem is unbounded.');
            STATE = "TERMINATE";
    end
end

% Remove initial dummy column
x_trace = x_trace(:,2:end);
basis_idx_trace = basis_idx_trace(:,2:end);

% --------------------------------------------------
% GRAPHICAL REPRESENTATION
% --------------------------------------------------

x_sample = -2:0.01:10;
y_sample = -2:0.01:10;

[X,Y] = meshgrid(x_sample,y_sample);

% Feasible region definition
const = (X - Y <= 2) & ...
        (X + Y <= 6) & ...
        (X >= 0) & ...
        (Y >= 0);

figure
contourf(X,Y,const,[1 1],'FaceAlpha',0.3)
colormap(summer)
axis equal
xlim([-0.5 7])
ylim([-0.5 7])
xlabel('x_1')
ylabel('x_2')
title('FEASIBLE REGION AND SIMPLEX PATH')
hold on
grid on

% --------------------------------------------------
% RECONSTRUCT (x1,x2) FROM BASIS HISTORY
% --------------------------------------------------

x_plot = zeros(2,iter);

for k = 1:iter

    % Extract x1
    if any(basis_idx_trace(:,k) == 1)
        idx_1 = find(basis_idx_trace(:,k) == 1);
        x_1_plot = x_trace(idx_1,k);
    else
        x_1_plot = 0;
    end

    % Extract x2
    if any(basis_idx_trace(:,k) == 2)
        idx_2 = find(basis_idx_trace(:,k) == 2);
        x_2_plot = x_trace(idx_2,k);
    else
        x_2_plot = 0;
    end

    x_plot(:,k) = [x_1_plot; x_2_plot];
end

% Plot simplex vertices
plot(x_plot(1,:), x_plot(2,:), ...
     'wo','MarkerFaceColor','w')

% Plot direction arrows
for k = 1:size(x_plot,2)-1

    x0 = x_plot(1,k);
    y0 = x_plot(2,k);

    dx = x_plot(1,k+1) - x0;
    dy = x_plot(2,k+1) - y0;

    quiver(x0, y0, dx, dy, 0, ...
           'r','LineWidth',2,'MaxHeadSize',0.5);
end

Z = -2*X -1*Y;
% Draw surface of Z over feasible region
surf(X, Y, Z, 'EdgeColor','none', 'FaceAlpha',0.8)
colormap(winter)
% Overlay contour of objective on feasible region
contour3(X, Y, Z, 10, 'k')

% Highlight optimal point if found
if exist('x_final','var')
    plot3(x_final(1), x_final(2), f*x_final, 'kh', 'MarkerSize',12, ...
          'MarkerFaceColor','y')
end
view(45,30)
colorbar
title('OPTIMAL SOLUTION OBTAIN BY SIMPLEX METHOD')
xlim([-0.5 7])
ylim([-0.5 7])
zlim([-15 1])