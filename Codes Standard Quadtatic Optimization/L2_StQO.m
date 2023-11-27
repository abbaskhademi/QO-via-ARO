function [LB, x_k, Time] = L2_StQO(Q)
% Lower/Upper bound by FME_StQO (L2_StQO)
format long
warning off

%__________StQO Problem Definition__________%
% min    x'Qx                               %
% s.t.   e'x = 1, x >= 0.                   %
%___________________________________________%

%=================================================================%
% ____________ INPUT:                                             %
% Q ........... Matrix of the objective function                  %
% ____________ OUTPUT:                                            %
% x_k .......... Worst-case scenario from (L2_StQO)               %
% LB ........... Lower bound value                                %
% Time ......... Total solver time                                %
%=================================================================%

% Requirements: Yalmip & Mosek

% Setting parameters
[~, n] = size(Q);
I = eye(n);

% Representation of Q
% Uncomment the desired representation:
% la = min(eig(Q)); Q1 = Q - ((la - 0.00001) * I); Q2 = (la - 0.00001) * I; % Representation 1: Smallest Eigenvalues of Q
la = max(eig(Q)); Q1 = ((la + 0.00001) * I); Q2 = Q - Q1; % Representation 2: Largest Eigenvalues of Q  

% Lower bound Fourier–Motzkin Elimination (L2) without dual, direct approach
S = zeros(n, 1);
Time = zeros(n, 1);
Scenario = zeros(n, n);
for k = 1:n
    clear('yalmip');
    x = sdpvar(n, 1);  
    q = Q2(k, :);
    P1 = optimize([sum(x) == 1, x >= 0], 0.5 * (x') * Q1 * x + q * x, sdpsettings('verbose', 0, 'solver', 'mosek'));
    if P1.problem ~= 0
        disp("Error in optimization:");
        yalmiperror(P1.problem)
    else
       x = value(x);
       Scenario(:, k) = x;
       Time(k) = P1.solvertime; % Solver time for each scenario
       S(k) = 0.5 * (x') * Q1 * x + q * x;
    end
end
Time_sena = sum(Time); % Total time to reach scenarios 

% Optimization for lower bound
clear('yalmip');
tau = sdpvar;
u = sdpvar(n, 1);
Cons = [];
for j = 1:n
    Cons = [Cons, S(j) - 0.5 * (u') * Q1 * u + Q1(j, :) * u >= tau];
end
sol = optimize(Cons, -tau, sdpsettings('verbose', 0, 'solver', 'mosek'));
if sol.problem ~= 0
    disp("Error in optimization:");
    yalmiperror(sol.problem)
else
    Time_LB = sol.solvertime;
    Time = sol.solvertime + Time_sena; % Total solver time
    LB = value(tau);
end

% Finding the best scenario and upper bound
upper = zeros(n, 1);
for i = 1:n 
    upper(i) = (Scenario(:, i))' * Q * (Scenario(:, i));
end
[UB0, k] = sort(upper);
UB = UB0(1);
x_k = Scenario(:, k(1)); % Best scenario L2_StQO

% Output display
Approach = {'Lower (L2_StQO)'; 'Best Scenario Value (L2_StQO)'};
Objective_value = [LB; UB];
time = [Time_LB; Time_sena];
On_L2_StQO = table(Approach, Objective_value, time) % Print the result if needed
end
