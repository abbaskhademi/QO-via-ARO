function [cplex_value,time]=StQO_cplex(Q,time_limit)
% Solves the Standard Quadratic Optimization Problem using CPLEX solver

%__________StQO__________%
% min    x'Qx            %
% s.t.   e'x=1, x>=0.    %
%________________________%

%====================================================%
% ____________ INPUT:                                %
% Q ......... matrix of the objective function       %
% time_limit....... time limite to restriced solver  %
% ____________ OUTPUT:                               %
% cplex_value ......... upper bound value by cplex %
% time ......... solver time to reach to UB          %
%====================================================%

[~,n]=size(Q);                     % Dimension of the problem
e = ones(n,1);
%
yalmip clear;
x = sdpvar(n,1);                   % x is an n-dimensional optimization variable
X1 = [e'*x==1, x>=0 ];             % StQO constraints
% Optimization using CPLEX
ops = sdpsettings('solver','cplex','cplex.TimeLimit',time_limit,'verbose',2); % set solver and max time
sol = optimize(X1,x'*Q*x ,ops);
if sol.problem ~= 0
    disp("Error in optimization:");
    yalmiperror(sol.problem)
else 
    x_ture = value(x);                              % optimal sotution
    cplex_value = x_ture'*Q*x_ture;                % optimal value
    cplex_time = sol.solvertime; 
    % Displaying the result
    Approch = {'Optimal(cplex)'};
    optimal_value = cplex_value;
    time = cplex_time;
    On_cplex = table(Approch,optimal_value,time)  % Print the result
end
%
end
