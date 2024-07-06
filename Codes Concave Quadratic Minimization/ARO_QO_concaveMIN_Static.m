function [LB,UB,Total_Time,sol_x]=ARO_QO_concaveMIN_Static(Q,c,A,b) 

% ARO-QO Algorithm: with full static decision rule for Concave Quadratic Minimization
% This function solves a concave quadratic minimization problem with constraints.
% The problem is defined as:
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
%                                                                 %
% ____________ OUTPUT:                                            %
% LB ........ Lower Bound                                         %
% UB ........ Upper Bound                                         %
% sol_x ..... Candidate solution based ARO QO                     %
% Total_Time ...... Total time to reach to Lower and Upper bounds %
%=================================================================%

% Requirements: Yalmip & Gurobi

format long
[m_x,n_x]=size(A);

% Define variables for the Lower Bound (LB) problem
clear('yalmip');
tau = sdpvar;
z = sdpvar(m_x,1);           
alpha = sdpvar(m_x,1); 
Pi = sdpvar(m_x,n_x,'full');  
T =  sdpvar(m_x,m_x,'full');

% Define constraints
C = [];
C = C + [b'*(alpha+z)>=tau];
C = C + [A'*alpha<=0.5*c];
C = C + [(b'*T)+z'>=0];
C = C + [A'*T<=zeros(m_x,n_x)'];
C = C + [A'*Pi<=(Q)'];
C = C + [b'*Pi+(0.5*c')-(A'*z)'>=0];
C = C + [alpha>=0, Pi(:)>=0, T(:)>=0];

% Solve the LB problem
sol = optimize(C,-tau,sdpsettings('verbose', 0,'solver','gurobi'));
% Analyze error flags
if sol.problem ~= 0
    disp("Error in optimization:");
    yalmiperror(sol.problem)
else 
% Extract LB value 
     LB_Time = sol.solvertime;  
     LB = value(tau);     
     Algorithm = {'(ARO QO): LB with Static Decision Rule'};
     On_QP = table(Algorithm, LB, LB_Time);
     tau = value(tau); z = value(z);    % the optimal solution of LB problem if needed
end
%%
% Upper Bound (UB) Step: collect scenarios based on robust counterpart
clear('yalmip');
x=sdpvar(n_x,1);
sol=optimize([A*x>=b, x>=0],(0.5*c')*x ,sdpsettings('verbose', 0,'solver','gurobi'));
my_time =  sol.solvertime;
UB= value(x'*Q*x+c'*x);
x0=value(x);                   % Collect scenario
S=[];
S=[S,x0];
%
for i=1:n_x
    clear('yalmip');
    x=sdpvar(n_x,1);
    sol=optimize([A*x>=b, x>=0], Q(i,:)*x ,sdpsettings('verbose', 0,'solver','gurobi'));
    my_time = my_time + sol.solvertime;
    x0=value(x);
    UB2= value(x'*Q*x+c'*x);
    if UB2<UB   % Check the quality of scenarios and select the best one
       UB=UB2;
       S=[S,x0];
    end
end
%%
% Mountain Climbing Procedure
x0 =S(:,end);   % Set initial point
while 1
    clear('yalmip');
    y=sdpvar(n_x,1);
    p3=optimize([A*y>=b, y>=0],x0'*Q*y+0.5*c'*y ,sdpsettings('verbose', 0,'solver','gurobi'));
    my_time = my_time + p3.solvertime;
    x0=value(y);
    S=[S,x0];
    if abs(S(:,end)'*Q*S(:,end)+c'*S(:,end)-S(:,end-1)'*Q*S(:,end-1)-c'*S(:,end-1))<=0.00001
        UB= (S(:,end)'*Q*S(:,end)+c'*S(:,end));
        sol_x = S(:,end);    % Solution candidate
        UB_Time=my_time;
        Algorithm = {'(ARO_QO): UB'};
        On_QP = table(Algorithm,UB,UB_Time)
        Total_Time = my_time + LB_Time 
        break
    end
end
    %            
end
