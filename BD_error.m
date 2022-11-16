function [err] = BD_error(kplus,kminus,C,tspan,tdiscrete)
global stiff stopping
tic
if stiff
    [~,y] = ode15s(@(t,y) BD_ODE(t,y,kplus,kminus), tspan, C(:,1));
else
    [~,y] = ode45(@(t,y) BD_ODE(t,y,kplus,kminus), tspan, C(:,1));
end
if length(tspan)==length(tdiscrete)
    err = norm(y'-C,"fro");
else
    BDdata = interp1(tspan,y,tdiscrete);
    err = norm(BDdata'-C,"fro");
end
    if toc > 1 && stiff == 0
        stiff = 1;
    elseif toc > 1 && stiff == 1
        stopping = 1;
    end
end