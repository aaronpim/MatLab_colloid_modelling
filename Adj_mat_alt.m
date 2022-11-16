function [Adjacency,Comp] = Adj_mat_alt(error,directory,radius,time,No_coll)
if nargin < 3
    time = 0:10000:400000;
    No_coll = 100;
    radius = 5;
end
Adjacency = zeros(No_coll,No_coll,length(time));
Comp = zeros(No_coll,length(time));
for i = 1:length(time)
    if time(i) == 0
        Fname = strcat(directory,'/col-cds00000000.csv');
    elseif time(i) <100000
        Fname = strcat(directory,'/col-cds000',int2str(time(i)),'.csv');
    else
        Fname = strcat(directory,'/col-cds00',int2str(time(i)),'.csv');
    end
    [X] = importfile(Fname, 2, 468);
    N = height(X);
    A=zeros(N,N);
    for k=1:N
        for j=1:N
            A(k,j) = (sqrt((X(k,1)-X(j,1))^2+(X(k,2)-X(j,2))^2+(X(k,3)-X(j,3))^2)<2*radius + error);
        end
    end
    Adjacency(:,:,i)=A-eye(N,N);
    [~,b] = conncomp(graph(A-eye(N,N)));
    for j = 1:No_coll
        Comp(j,i) = sum(b==j);
    end
end
end

