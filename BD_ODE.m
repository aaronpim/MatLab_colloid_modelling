function dydt = BD_ODE(t,y,kplus,kminus)
    k = abs(kplus);
    q = abs(kminus);
    N = length(y);
    dydt = zeros(N,1);
    dydt(1) = -k(1).*y(1).*y(1)+q(1)*y(2);
    for i = 2:N-1
        dydt(i) = -k(i).*y(i).*y(1)+q(i)*y(i+1)+k(i-1).*y(i-1).*y(1)-q(i-1)*y(i);
    end
    dydt(N) = k(N-1).*y(N-1).*y(1)-q(N-1)*y(N);
end