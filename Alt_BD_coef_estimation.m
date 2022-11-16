t = 0:10000:400000;
epsilon = 1.0e-5;
tol = 1.0e-5;
name_cell = {...
    'Coll_exp_0',  'Coll_exp_0p5','Coll_exp_1',  'Coll_exp_1p5',...
    'Coll_exp_1p6','Coll_exp_1p7','Coll_exp_1p8','Coll_exp_1p9',...
    'Coll_exp_2',  'Coll_exp_2p1','Coll_exp_2p2','Coll_exp_2p3',...
    'Coll_exp_2p4','Coll_exp_2p5','Coll_exp_2p6','Coll_exp_2p7',...
    'Coll_exp_2p8','Coll_exp_2p9','Coll_exp_3',  'Coll_exp_3p1',...
    'Coll_exp_3p2','Coll_exp_3p3','Coll_exp_3p4','Coll_exp_3p5',...
    'Coll_exp_4',  'Coll_exp_4p5','Coll_exp_5'};

W_vec = exp([0, 0.5, 1, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4,...
 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4 3.5, 4, 4.5, 5]);

kp_cell = cell(length(name_cell),1);
km_cell = cell(length(name_cell),1);
LB_Simulation = cell(length(name_cell),1);
BD_Simulation = cell(length(name_cell),1);

[C,BDdata,kplus,kminus] = BD_coef_estimation(strcat('Wall_0/',name_cell{4}),epsilon,tol);
kp_cell{4} = kplus;
km_cell{4} = kminus;
LB_Simulation{4} = C';
BD_Simulation{4} = BDdata;
save('Alt_BD_Coef.mat','kp_cell','km_cell','LB_Simulation','BD_Simulation','W_vec')
for i = 4:length(name_cell)
    [C,BDdata,kplus,kminus] = BD_coef_estimation(strcat('Wall_0/',name_cell{i}),epsilon,tol,k0);
    kp_cell{i} = kplus;
    km_cell{i} = kminus;
    LB_Simulation{i} = C';
    BD_Simulation{i} = BDdata;
    k0 = [kplus;kminus];
    save('Alt_BD_Coef.mat','kp_cell','km_cell','LB_Simulation','BD_Simulation','W_vec')
end
k0 = [kp_cell{4};km_cell{4}];
[C,BDdata,kplus,kminus] = BD_coef_estimation(strcat('Wall_0/',name_cell{3}),epsilon,tol,k0);
kp_cell{3} = kplus; km_cell{3} = kminus; LB_Simulation{3} = C'; BD_Simulation{3} = BDdata;
k0 = [kplus;kminus];
save('Alt_BD_Coef.mat','kp_cell','km_cell','LB_Simulation','BD_Simulation','W_vec')
[C,BDdata,kplus,kminus] = BD_coef_estimation(strcat('Wall_0/',name_cell{2}),epsilon,tol,k0);
kp_cell{2} = kplus; km_cell{2} = kminus; LB_Simulation{2} = C'; BD_Simulation{2} = BDdata;
k0 = [kplus;kminus];
save('Alt_BD_Coef.mat','kp_cell','km_cell','LB_Simulation','BD_Simulation','W_vec')
[C,BDdata,kplus,kminus] = BD_coef_estimation(strcat('Wall_0/',name_cell{1}),epsilon,tol,k0);
kp_cell{1} = kplus; km_cell{1} = kminus; LB_Simulation{1} = C'; BD_Simulation{1} = BDdata;
k0 = [kplus;kminus];
save('Alt_BD_Coef.mat','kp_cell','km_cell','LB_Simulation','BD_Simulation','W_vec')