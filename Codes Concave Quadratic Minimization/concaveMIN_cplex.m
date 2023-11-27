function [UB,Total_Time,sol_x]=concaveMIN_cplex(Q,A,b,c,time_limit)


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
% UB ........ Upper Bound                                         %
% sol_x ..... Candidate solution based Time limit                 %
% Total_Time ...... Total time to reach to the bound              %
%=================================================================%

format long
[~,n_x]=size(A);
clear('yalmip');
x=sdpvar(n_x,1);
X1 = [A*x>=b, x>=0 ];

obj = x'*Q*x+c'*x;

ops = sdpsettings('solver','cplex','cplex.TimeLimit',time_limit,'verbose',2);
sol = optimize(X1,obj ,ops);
if sol.problem ~= 0
    disp("error!")
    sol.info
    yalmiperror(sol.problem)
else 
    UB = value(obj);                       % optimal value
    cplex_time = sol.solvertime; 
    sol_x = value(x);
    Approch = {'CPLEX'};
    optimal_value = UB;               
    Total_Time = cplex_time;
    On_cplex = table(Approch,optimal_value,Total_Time)
end

end