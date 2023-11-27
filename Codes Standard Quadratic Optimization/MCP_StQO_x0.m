function [x_end, u_end, time] = MCP_StQO_x0(Q, x_0, time_x_0)
% Mountain Climbing Procedure (MCP) based on Bi_StQO: Improve the initial solution
format long
warning off

%________________Bi_StQP Problem Definition_________________%
% min   0.5 * {x' * (Q+) * x + y' * (Q+) * y} + x' * (Q-) * y %
% s.t.   e' * x = 1, x >= 0, e' * y = 1, y >= 0.               %
%___________________________________________________________%

%====================================================%
% ____________ INPUT:                                %
% Q .......... Matrix of the objective function      %
% x_0 ......... Initial point                        %
% time_x_0 .... Time to reach initial point          %
% ____________ OUTPUT:                               %
% x_end ....... Candidate solution                   %
% u_end ....... Upper bound value                    %
% time ........ Total solver times to reach UB       %
%====================================================%

% Requirements: Yalmip & Mosek

% Setting parameters
[~, n] = size(Q);
I = eye(n);

% Selecting the Representation of Q
% Uncomment the desired representation:
% la = min(eig(Q)); Q1 = Q - ((la - 0.00001) * I); Q2 = (la - 0.00001) * I; % Representation 1: Smallest Eigenvalues of Q
la = max(eig(Q)); Q1 = ((la + 0.00001) * I); Q2 = Q - Q1; % Representation 2: Largest Eigenvalues of Q  

% Initializing variables for MCP
y_matrix = [x_0];
y_value = [x_0' * Q * x_0];
time_y_cal = 0; 

% MCP iteration
while true
    xe = y_matrix(:, end);
    clear('yalmip');
    y = sdpvar(n, 1); % y is an n-dimensional optimization variable

    % Optimization problem
    sol = optimize([sum(y) == 1, y >= 0], 0.5 * y' * Q1 * y + xe' * Q2 * y, sdpsettings('solver', 'mosek', 'verbose', 0));
    
    % Check for optimization problem
    if sol.problem ~= 0
        disp("Error in optimization:");
        yalmiperror(sol.problem)
    else
        y_matrix = [y_matrix, value(y)];
        y_value = [y_value, value(y)' * Q * value(y)];
        time_y_cal = time_y_cal + sol.solvertime;
    end

    % Convergence check
    if abs(y_value(end) - y_value(end - 1)) < 0.0000001
        x_end = y_matrix(:, end);
        u_end = y_value(end);
        Approach = {'Mountain Climbing Procedure: x_0'};
        optimal_value = u_end; 
        time = time_y_cal + time_x_0; % Total time
        MCP = table(Approach, optimal_value, time); % Print the result
        break;
    end
end
end
