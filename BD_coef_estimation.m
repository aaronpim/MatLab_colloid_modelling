function [C,BDdata,kplus,kminus] = BD_coef_estimation(name,epsilon,tol,k0)
%%
% This code takes positions of colloids over time, computes the change in
% component size over time, and then computes Smoluchoski coefficients to
% fit the data.

% The variable 'name' takes type of string and denotes the name of the folder
% that the files are contained in.

% The variable 'epsilon' takes a data type of double and must be positive.
% This determines the cut off point for the minimisation scheme. This does
% repeatedly runs the fminsearch function to optimise the coefficients,
% using the minimum value from the previous run as an initial value for the
% new one. If the change in error from one run to the next is less than
% epsilon then the code terminates.
if nargin == 0
    name = 'Wall_0/Coll_exp_2'; epsilon = 1.0e-4; tol = 1e-4;
elseif nargin == 1
    name = 'Wall_0/Coll_exp_2'; epsilon = 1.0e-4;
elseif nargin == 2
    name = 'Wall_0/Coll_exp_2';
end
%% 
%This section takes the csv files from the folders and computes the number
%of components of each size for each time. This is the matrix C.
[~,Comp] = Adj_mat(2.5,name);
B = sum(Comp,2);
for n = 0:(length(B)-1)
    if B(length(B)-n)~=0
        break
    end
end
C = Comp(1:(length(B)-n),:);
% C(i,j) = number of simulated aggregates of size i and timestep j.
%%
N = height(C);
tspan = 0:10000:400000;
tdiscrete = 0:10000:400000;
if nargin < 4
    k0 = zeros(2*N-2,1);
elseif length(k0) < 2*N-2
    n = length(k0)/2;
    k0 = [k0(1:n);zeros(N-1-n,1);k0((1:n)+n);zeros(N-1-n,1)];
elseif length(k0) > 2*N-2
    n = length(k0)/2;
    k0 = [k0(1:N-1);k0((1:N-1)+n)];
end
global stiff
stiff = 0;
fun = @(k) BD_error(k(1:N-1),k(N:2*N-2),C,tspan,tdiscrete);
options = optimset('Display','iter','MaxFunEvals',50000,'MaxIter',50000,'TolFun',tol,'TolX',tol);
k_min_vec = abs(fminsearch(fun,k0,options));
er_old = 0;
er = fun(k_min_vec);
global stopping
stopping = 0;
N_max = 100;
n = 0;
while abs(er-er_old)>epsilon && stopping==0 && n<N_max 
    er_old = er;
    k_min_vec = abs(fminsearch(fun,k_min_vec,options));
    er = fun(k_min_vec);
    n = n+1;
end
kplus = k_min_vec(1:N-1); kminus = k_min_vec((1:(N-1))+(N-1));
if stiff
    [~,y] = ode15s(@(t,y) BD_ODE(t,y,kplus,kminus), tspan, C(:,1));
else
    [~,y] = ode45(@(t,y) BD_ODE(t,y,kplus,kminus), tspan, C(:,1));
end
BDdata = interp1(tspan,y,tdiscrete);
end