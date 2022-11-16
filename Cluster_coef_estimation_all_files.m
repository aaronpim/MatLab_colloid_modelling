name_cell = {'Coll_exp_0','Coll_exp_0p5','Coll_exp_1','Coll_exp_1p5','Coll_exp_2','Coll_exp_2p5','Coll_exp_3','Coll_exp_3p5','Coll_exp_4','Coll_exp_4p5','Coll_exp_5'};
W_vec = exp([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]);
Time = 0:10000:400000;
for i = 1:length(name_cell)
    if i==1
        C = Clustering_coeff_cal(strcat('Wall_exp_1/',name_cell{i}));
        Clustering_Tensor_1 = zeros(height(C),width(C),length(name_cell));
        Adjacency_Tensor_1  = zeros(height(C),height(C),width(C),length(name_cell));
        Clustering_Tensor_1(:,:,i) = C;
        [A,~] = Adj_mat_no_WD(2.5,strcat('Wall_exp_1/',name_cell{i}));
        Adjacency_Tensor_1(:,:,:,i)= A;
    else
        C = Clustering_coeff_cal(strcat('Wall_exp_1/',name_cell{i}));
        Clustering_Tensor_1(:,:,i) = C;
        [A,~] = Adj_mat_no_WD(2.5,strcat('Wall_exp_1/',name_cell{i}));
        Adjacency_Tensor_1(:,:,:,i)= A;
    end
    save('cluster_coef_wall_1.mat','name_cell','W_vec','Time','Clustering_Tensor_1','Adjacency_Tensor_1')
end
mean_clustering_tensor_1 = reshape(mean(Clustering_Tensor_1,1),41,11);
save('cluster_coef_wall_1.mat','name_cell','W_vec','Time','Clustering_Tensor_1','Adjacency_Tensor_1','mean_clustering_tensor_1')
[T,Wlog] = meshgrid(Time', log(W_vec'));
s = pcolor(T,Wlog,mean_clustering_tensor_1');
s.EdgeColor = 'none';
shading interp;
ylabel('$\log W$','interpreter','latex')
xlabel('Time','interpreter','latex')
c = colorbar; c.Label.String = 'Clustering Coeficient'; c.Label.Interpreter = 'latex';

rescale = 10^6;
rate  = @(x) x(1).*(1-exp(-(Time')*x(2)/rescale));
options = optimset('Display','iter','MaxFunEvals',100000,'MaxIter',100000,'TolFun',1.0e-8,'TolX',1.0e-8);
a0_limit = zeros(length(name_cell),1);
lambda   = zeros(length(name_cell),1);
for i = 1:length(name_cell)
    er = @(x) norm(rate(x)-mean_clustering_tensor_1(:,i));
    x_min = fminsearch(er,[1,1],options);
    a0_limit(i) = x_min(1);
    lambda(i) = x_min(2)/rescale;
end
save('cluster_coef_wall_1.mat','name_cell','W_vec','Time','Clustering_Tensor_1','Adjacency_Tensor_1','mean_clustering_tensor_1','a0_limit','lambda')
%j = 3;
%plot(Time, rate([a0_limit(j),lambda(j)*rescale]),'b--', Time, mean_clustering_tensor_2(:,j),'k-')
degree = reshape(sum(Adjacency_Tensor_1,1),100,41,11);
deg_freq = zeros(13,41,11);
for i=1:11
    for t=1:41
        for j = 1:13
            deg_freq(j,t,i) = sum(degree(:,t,i)==j-1);
        end
        deg_freq(:,t,i) = deg_freq(:,t,i)/sum(deg_freq(:,t,i));
    end
end
[D,T]=meshgrid(0:12,Time);
surf(D',T',reshape(deg_freq(:,:,11),13,41))
mean_degree = reshape(mean(degree,1),41,11);
save('cluster_coef_wall_1.mat','name_cell','W_vec','Time','Clustering_Tensor_1','Adjacency_Tensor_1','mean_clustering_tensor_1','a0_limit','lambda','degree','deg_freq','mean_degree')