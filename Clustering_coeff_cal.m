function [C] = Clustering_coeff_cal(name)
try 
    [A,~] = Adj_mat(2.5,name);
catch
    [A,~] = Adj_mat_no_WD(2.5,name);
end
[N,~,T] = size(A);
C = zeros(N,T);
for t = 1:T
    k = sum(A(:,:,t),1).*(sum(A(:,:,t),1)-1);
    a = reshape(A(:,:,t),N,N);
    a3 = a^3;
    for i=1:N
       C(i,t) = a3(i,i)./k(i);
    end    
end
C(isnan(C))=0;
end