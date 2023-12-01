function [LB,UB,Total_Time,sol_x]=concaveMIN_gurobi(Q,A,b,c,time_limit)


%===============================%
%	    min    x'Qx + c'x       %
%	    s.t.     Ax >= b,       %
%                 x >= 0.       %
%===============================%

%=================================================================%
% ____________ INPUT:                                             %
% Q ......... the matrix associated with the objective function   %
% c ......... a column vector associated with the linear part of  %
% the objective function                                          %
% A ......... the matrix associated with the constraints          %
% b .........  a column vector associated RHS of the constraints  %
% time_limit ...... Maximum of solver time                        %
%                                                                 %
% ____________ OUTPUT:                                            %
% LB ........ Lower Bound                                         %
% UB ........ Upper Bound                                         %
% sol_x ..... Candidate solution based ARO QO                     %
% Total_Time ...... Total time to reach to Lower and Upper bounds %
%=================================================================%

% Requirements: Yalmip & Gurobi

format long
[~,n_x]=size(A);
%%

clear('yalmip');
x=sdpvar(n_x,1);
X1 = [A*x>=b, x>=0 ];
obj = x'*Q*x+c'*x;
ops = sdpsettings('solver','gurobi','verbose',1, 'gurobi.TimeLimit',time_limit,'savesolveroutput',1);
sol = optimize(X1,obj ,ops);
if sol.problem ~= 0
    disp("error!")
    sol.info
    yalmiperror(sol.problem)
else 
    UB = value(obj);   % optimal value
    LB=sol.solveroutput.result.objbound;
    sol_x = value(x);
    Approch0 = {'Gurobi'};
    optimal_value = UB;                
    Total_Time = sol.solvertime; 
    On_gurobi = table(Approch0,optimal_value,Total_Time)
end
%
end
