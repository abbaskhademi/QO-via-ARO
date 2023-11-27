function [opt_ipopt,time1]=StQO_local_solver(Q)
%  Upper bound by local solver

%__________StQO__________%
% min    x'Qx            %
% s.t.   e'x=1, x>=0.    %
%________________________%

%========================================================%
% ____________ INPUT:                                    %
% Q ......... matrix of the objective function           %
% ____________ OUTPUT:                                   %
% opt_ipopt .......... Upper bound value by IPOPT        %
% time1 ......... solver time to reach to UB by IPOPT    %
%========================================================%

% Requirements: Yalmip & IPOPT

[~,n]=size(Q);                    % dimension
e = ones(n,1);
%%                            ipopt
yalmip clear;
x = sdpvar(n,1);                  % x is n times 1 (optimization variable)
X1 = [e'*x==1, x>=0 ];            % StQO constraints
ops = sdpsettings('solver','ipopt','verbose',1,'ipopt.max_iter',100000);   % set IPOPT as solver
g0 = optimize(X1,x'*Q*x ,ops);
if g0.problem ~= 0
    disp("Error in optimization:")
    g0.info
    yalmiperror(g0.problem)
else
    x_ture = value(x);                                % (sub)optimal sotution
    ipopt_value = x_ture'*Q*x_ture;                   % objective value
    ipopt_time = g0.solvertime; 
    optimal_value1 = ipopt_value;
    time1 = ipopt_time;
    opt_ipopt=optimal_value1;
    Approch = {'UB (IPOPT)'};
    UB =[optimal_value1];
    Time = [time1];
    Local = table(Approch,UB,Time)                    % Print the result
end


%
end

