function [gurobi_value,LB_gurobi,time]=StQO_gurobi(Q,time_limit)
% Solves the Standard Quadratic Optimization Problem using Gurobi solver
 
%__________StQO__________%
% min    x'Qx            %
% s.t.   e'x=1, x>=0.    %
%________________________%

%====================================================%
% ____________ INPUT:                                %
% Q ......... matrix of the objective function       %
% time_limit....... time limite to restriced solver  %
% ____________ OUTPUT:                               %
% gurobi_value ......... upper bound value by Gurobi %
% time ......... solver time to reach to UB          %
%====================================================%

[~,n]=size(Q);                     % Dimension of the problem
e = ones(n,1);
%
yalmip clear;
x = sdpvar(n,1);                   % x is an n-dimensional optimization variable
X1 = [e'*x==1, x>=0 ];             % StQO constraints
%
ops = sdpsettings('solver','gurobi','verbose',0, 'gurobi.TimeLimit',time_limit,'savesolveroutput',1); % set solver and max time
sol = optimize(X1,x'*Q*x ,ops);
if sol.problem ~= 0
    disp("Error in optimization:");
    yalmiperror(sol.problem)
else 
    x_ture = value(x);                              % optimal sotution
    gurobi_value = x_ture'*Q*x_ture;                % optimal value
    gurobi_time = sol.solvertime; 
    LB_gurobi=sol.solveroutput.result.objbound;     % Lower bound
    % Displaying the result
    Approch = {'Optimal(gurobi)'};
    optimal_value = gurobi_value;
    time = gurobi_time;
    On_Gurobi = table(Approch,optimal_value,time)  % Print the result
end
%
end




