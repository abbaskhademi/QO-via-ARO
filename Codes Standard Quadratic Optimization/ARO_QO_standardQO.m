function [LB1,LB2,x_up,UB,Time] = ARO_QO_standardQO(Q)
% ARO-QO Algorithm for Standard Quadratic Optimization Problem (StQO)   

%__________StQO Problem Definition__________%
% min    x'Qx                               %
% s.t.   e'x = 1, x >= 0.                   %
%___________________________________________%

%====================================================%
% ____________ INPUT:                                %
% Q ........... Matrix of the objective function     %
% ____________ OUTPUT:                               %
% x_up ......... Candidate solution                  %
% UB ........... Upper bound value                   %
% LB1 .......... Lower bound value from L1-StQO      %
% LB2 .......... Lower bound value from L2-StQO      %
% Time ......... Total solver time to reach bounds   %
%====================================================%

[LB2, x_2, T_2] = L2_StQO(Q);              % Worst-case scenarios from EME_StQP (L2-StQO) 
disp("Initial solution (x_2 MC):");
[x2, l2, t2] = MCP_StQO_x0(Q, x_2, T_2);   % Improve the initial solution

[LB1, x_1, T_1] = L1_StQO(Q);              % Worst-case scenarios from L1-StQO
disp("Initial solution (x_1 MC):");
[x1, l1, t1] = MCP_StQO_x0(Q, x_1, T_1);   % Improve the initial solution

Time = t1 + t2;                            % Sum of total solver times 
[UB, ind] = min([l1, l2]);                 % Pick the best objective value as the UB value
X = [x1, x2];
x_up = X(:, ind);                          % Choose x_up as the solution with the best objective value

% Displaying the final results
Approach = {'ARO_QO_StQO'};
objective_value = UB;
time = Time;
On_StQO = table(Approach, objective_value, time) % Print the result

end
