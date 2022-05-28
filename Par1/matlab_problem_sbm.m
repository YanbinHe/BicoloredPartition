% reproducible code

%Part of these functions are based on the well-known GSP_Toolbox
%https://epfl-lts2.github.io/gspbox-html/
%Please cite this paper and theirs if you find this code useful.



warning off
clc
clear
% generate graph with 4 communities
load 'graph1.mat' % load coordinates. 
%This graph is not exactly what we showed in the paper, but the same
%conclusion can be obtained. Or you can generate your own graph using the
%functions provided here. 
%use creat_coords_m.m to generate coordinates.
%%
num_nodes = 40;

alpha = 0;
beta = 0;
intra_clu = 0.9-alpha;
inter_clu_homo = 0.2+alpha;
inter_clu_hetero = 0.1+beta;
P = [intra_clu inter_clu_homo 0 inter_clu_hetero;
    inter_clu_homo intra_clu inter_clu_hetero 0;
    0 inter_clu_hetero intra_clu inter_clu_homo;
    inter_clu_hetero 0 inter_clu_homo intra_clu];


G = gsp_sbm(num_nodes,P,XCoords, YCoords);


b = [zeros(20,1);ones(20,1)];
G.label_N = b;
G.sz = 50;
signal = ones(num_nodes,1);


f1 = figure;
gsp_plot_signal_ox(G,signal);% dipcting the graph
set(gca,'fontname','DejaVuSans')
print(f1,'Grpah','-dpng','-r600')

W = G.A;
L = diag(W*ones(num_nodes,1)) - W;
D = full(diag(W*ones(num_nodes,1)));
Db = diag(b);
Du = diag(1-b);
Wtilde = Db*W*Du;

Wbar = 2*Wtilde;
Ws = (Wbar + Wbar');
Lbar = 1*L;
k = ones(1,num_nodes)*Wtilde*ones(num_nodes,1);
q = (Wtilde + Wtilde')*ones(num_nodes,1);

warning on
num_tol = 1e-6;
[Ub, Eb] = eig(Lbar);
Ub2 = Ub*pinv(SQRT(Eb),num_tol);
[Ua,Ea] = eig(0.5*Ub2'*Ws*Ub2+0.5*(Ub2'*Ws*Ub2)');
V = Ub2*Ua;

Veff = V(:,abs(diag(Ea)) > num_tol);

fncCost = @(x) (x'*L*x)/(x'*Wtilde*x) + ( (1-x)'*L*(1-x) )/( (1-x)'*Wtilde*(1-x) );

Vbin = 1*(Veff > 0);
cst_genEig = zeros(size(Vbin,2),1);
for kk = 1:size(Vbin,2)
    cst_genEig(kk) = (fncCost(Vbin(:,kk)));
end
%%
cst_genEig = fliplr(cst_genEig');
f3 = figure;
plot(cst_genEig,'-r','LineWidth',2)
title('Sum Ratios')
xlabel('Eigenvector Index [k]')
ylabel('Sum of Ratios Value')
set(gca,'fontname','DejaVuSans')
print(f3,'m1_sweep','-dpng','-r600')
%%
Ea = diag(Ea);
Ea = Ea(abs(Ea) > num_tol);
[~, id_gEig] = max(Ea);% largest eigenvalue

v_gEigbin = Vbin(:,id_gEig);% solution

fins2 = figure;
G.label_N = b;
gsp_plot_signal_ox(G, v_gEigbin+1)% binary
set(gca,'fontname','DejaVuSans')
print(fins2,'m1_coloring','-dpng','-r600')

% regular spectral clustering
[Vsc3,egv3] = eig(L);
Vscbi3 = 1*(Vsc3(:,2)>0);

% compare two methods
fncCost(v_gEigbin)
fncCost(Vscbi3)

fins3 = figure;
gsp_plot_signal_ox(G, Vscbi3+1)% partition solution
set(gca,'fontname','DejaVuSans')
print(fins3,'m1_sc','-dpng','-r600')


%% useless
% compute accuracy
% desired clusters
% indicator_desired = [ones(10,1);zeros(20,1);ones(10,1)];
% 
% v_gEigbin = check_label(v_gEigbin);
% Vscbi1 = check_label(Vscbi1);
% Vscbi2 = check_label(Vscbi2);
% Vscbi3 = check_label(Vscbi3);

%     accuracy_prop(i,j) = sum(1*(v_gEigbin == indicator_desired))/num_nodes+accuracy_prop(i,j);
%     accuracy_sc1(i,j) = sum(1*(Vscbi1 == indicator_desired))/num_nodes+accuracy_sc1(i,j);
%     accuracy_sc2(i,j) = sum(1*(Vscbi2 == indicator_desired))/num_nodes+accuracy_sc2(i,j);
%     accuracy_sc3(i,j) = sum(1*(Vscbi3 == indicator_desired))/num_nodes+accuracy_sc3(i,j);

%     end
% end
% end
%% useless
% f1 = figure
% imagesc(beta_list,alpha_list,accuracy_prop/10);%alpha_list,beta_list,
% colorbar
% title('The proposed method')
% xlabel('\beta');
% ylabel('\alpha');
% set(gca,'fontname','DejaVuSans')
% print(f1,'proposed','-dpng','-r600')

% f2 = figure
% imagesc(alpha_list,beta_list,accuracy_sc1/10);
% colorbar
% title('Spectral clustering with normalized Laplacian')
% xlabel('\alpha');
% ylabel('\beta');
% set(gca,'fontname','DejaVuSans')
% print(f2,'sc_nor','-dpng','-r600')
% figure
% imagesc(alpha_list,beta_list,accuracy_sc2/10);
% colorbar
% f3 = figure
% imagesc(beta_list,alpha_list,accuracy_sc3/10);
% colorbar
% title('Spectral clustering with unnormalized Laplacian')
% xlabel('\beta');
% ylabel('\alpha');
% set(gca,'fontname','DejaVuSans')
% print(f3,'sc_unnor','-dpng','-r600')