function [LB, x_u, Time] = L1_StQO(Q)   
% Lower bound with w_x = Zx + z0, v_x = v (L1_StQP); And scenarios (L1_StQP)
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
% x_u .......... Worst-case scenario from (L1_StQP)               %
% LB ........... Lower bound value                                %
% Time ......... Solver time to reach the worst-case scenario     %
%=================================================================%

% Requirements: Yalmip & Mosek

% Setting parameters
[~, n] = size(Q);
e = ones(n, 1);
I = eye(n);

% Selecting the Representation of Q
% Uncomment the desired representation:
% la = min(eig(Q)); Q1 = Q - ((la - 0.00001) * I); Q2 = (la - 0.00001) * I; % Representation 1: Smallest Eigenvalues of Q
la = max(eig(Q)); Q1 = ((la + 0.00001) * I); Q2 = Q - Q1; % Representation 2: Largest Eigenvalues of Q  


% Lower bound with w_x = Zx + z0, v_x = v
clear('yalmip');
tau = sdpvar;      
u = sdpvar(n, 1);
z0 = sdpvar;
z = sdpvar(n, 1);
alpha = sdpvar(n, 1);
beta = sdpvar;
theta = sdpvar(n, 1);

% Setting up the model
Objective = -tau;
Constraints = [];
Constraints = [Constraints,(-0.5*alpha'*Q1*alpha) + beta + z0 - (0.5*(u')*Q1*u)>= tau];
Constraints = [Constraints,(e*beta)-(Q1*alpha)<=z];
Constraints = [Constraints, Q1*u+theta-e*z0>=0];
Constraints = [Constraints, e*theta'<=(-e*z'+Q2)'];

% Optimizing the constraints
sol = optimize(Constraints, Objective, sdpsettings('verbose', 0, 'solver', 'mosek'));

if sol.problem ~= 0
    disp("Error in optimization:");
    yalmiperror(sol.problem)
else
    L1_time = sol.solvertime;  
    LB = value(tau);
    % Optimal solutions for the LB problem
    z = value(z);
    theta = value(theta);
    z0 = value(z0);
    beta = value(beta);
    u = value(u);
    alpha = value(alpha);
end

% Collecting worst-case scenarios
clear('yalmip');
x = sdpvar(n, 1);
X1 = [e' * x == 1, x >= 0];
Time = 0;
OO = [];

ops = sdpsettings('solver', 'mosek', 'verbose', 0);
for i = 1:n
    opt_1 = I(:, i)' * Q * I(:, i);
    OO = [OO, opt_1]; % Compute objective value at e_j's
end

% Finding scenario L1 procedure
g2 = optimize(X1, (0.5 * x' * Q1 * x) + z' * x, ops);
x_2 = value(x);
opt_2 = x_2' * Q * x_2;
OO = [OO, opt_2];
[UB1, index] = sort(OO); % Best objective value all scenarios from L1
scenario = [I, x_2];
x_u = scenario(:, index(1)); % Best scenario from L1

% Total time to reach upper bound with the L1 procedure
Time = Time + g2.solvertime + L1_time;

% Output display
Approach = {'Lower (L1_StQO)'; 'Best Scenario Value (L1_StQO)'};
optimal_value = [LB; UB1(1)]; 
time = [L1_time; Time];
On_L1_StQO = table(Approach, optimal_value, time) % Print the result if needed
end
