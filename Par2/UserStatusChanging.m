% clc
% clear
%% 
% input description:
% W is the initial adjacency matrix
% L is the initial laplacian matrix

% output: we can plot graph here in this function, so currently we don't
% have return anything

% new status should be a fat matrix \in R^{#base stations \times #users}
% we keep adding users to the cellular network using users in the newStatus matrix
%% definition of some basic values
epi = 1e-4;
numClu = 4;% #tracking components
numC = 3;% the number of clusters
numN = 37;% the number of base stations
num_nodes = size(Wini, 1);

b = zeros(1,size(Wini,1));
b([1:37]) = 1;
Db = diag(b);
u = 1 - b;
Du = diag(u);
Wtilde = Db*Wini*Du;    
Wbar = 2*Wtilde;
q = (Wtilde + Wtilde')*ones(num_nodes,1);

userUpEach = 20;% the number of users during each update
numSimu = 1200;
time = zeros(2,numSimu/userUpEach);
%% computing different W-related and L-related matrix
W = (Wbar + Wbar');
L = 2*Lini;

num_tol = 1e-6;
[Ub, Eb] = eig(L);
Ub2 = Ub*pinv(SQRT(Eb),num_tol);
[Ua,Ea] = eig(0.5*((Ub2'*W*Ub2).'+(Ub2'*W*Ub2)));
V = Ub2*Ua;
Veff = V(:,abs(diag(Ea)) > num_tol);

W = sparse(W);
L = sparse(L);
% large <- small
GenEigV = Veff(:,end-numClu + 1:end); % tracked components
EigV = diag(Ea(end-numClu + 1:end,end-numClu + 1:end)); % tracked eigenvalues
% %% give initial graph partitions
% labelIni = kmeans(GenEigV(:,end-numC + 1:end),numC,'Replicates',10);
%% give initial graph partitions
% for the sake of running speed, this part can be ignored
labelIni = kmeans(GenEigV(:,end-numC + 1:end),numC,'Replicates',10);
% figure
% for i = 1:num_bs_g
%     X(:,i) = distance0(1,i) + unit*cos(t);
%     Y(:,i) = distance0(2,i) + unit*sin(t);
% end
% for i = 1:num_bs_g
%     plot(X(:,i),Y(:,i),'k')
%     hold on
% end
% for i = 1:1200
%     if labelIni(i+37) == 1
%     scatter(user_coord(1,i),user_coord(2,i),'b.');
%     hold on
%     end
% if labelIni(i+37) == 2
%     scatter(user_coord(1,i),user_coord(2,i),'r.');
%     hold on
% end
%     if labelIni(i+37) == 3
%     scatter(user_coord(1,i),user_coord(2,i),'g.');
%     hold on
%     end
% end
% set(gca,'fontname','DejaVuSans');
%%
GenEigVi = GenEigV;
EigVi = EigV;
%%
% --------------------initialization ends-------------------------------- %
for simu = 1:numSimu/userUpEach
simu
num_nodes = size(Wini, 1);
num_users = num_nodes - numN;

newCol = newStatus(:,(simu-1)*userUpEach+1:simu*userUpEach);% new data
Worinew = Wini;
Worinew(1:37,(simu-1)*userUpEach+38:simu*userUpEach+37) = newCol;
Worinew((simu-1)*userUpEach+38:simu*userUpEach+37,1:37) = newCol';

tic;
Wtildenew = Db*Worinew*Du;    
Wbarnew = 2*Wtildenew;
qnew = (Wtildenew + Wtildenew')*ones(num_nodes,1);
Lorinew = diag(Worinew*ones(num_nodes,1)) - Worinew;

Wnew = sparse(Wbarnew + Wbarnew');
Lnew = sparse(2*Lorinew);    


dW = sparse(Wnew - W);
dL = sparse(Lnew - L);
dx = sparse(zeros(size(Lnew,2),numClu));

dlambda = zeros(numClu,1);

for kkk = 1:numClu % for each component
for iii = 1:2     
    
        dlambda(kkk) = ((GenEigV(:,kkk))'*(dW-(EigV(kkk))*dL)*(GenEigV(:,kkk)+dx(:,kkk)))/((GenEigV(:,kkk))'*(L+dL)*(GenEigV(:,kkk)+dx(:,kkk)));
        
        K = sparse(W+dW-(EigV(kkk)+ dlambda(kkk))*(L+dL));
        h = sparse((dlambda(kkk)*(L)+(EigV(kkk))*dL+dlambda(kkk)*dL-dW)*(GenEigV(:,kkk)));

        qk = sparse(K'*h);
        % conjugate gradient descent
        dx(:,kkk) = CG(K,epi*speye(size(K,1),size(K,1)),qk,sparse(zeros(size(dx(:,1)))));        
end
end

GenEigV = GenEigV + dx;
EigV = EigV + dlambda;
[EigV,order] = sort(EigV,'ascend');
GenEigV = GenEigV(:,order);

label3 = kmeans(GenEigV(:,end-numC + 1:end),numC,'Replicates',10);% using approximate

time(1,simu) = toc;

% exact computation
Wnew = full(Wnew);
Lnew = full(Lnew);

tic;
Wtildenew = Db*Worinew*Du;    
Wbarnew = 2*Wtildenew;
qnew = (Wtildenew + Wtildenew')*ones(num_nodes,1);
Lorinew = diag(Worinew*ones(num_nodes,1)) - Worinew;


% filled with code generated Ws and Lbar

Wnew = (Wbarnew + Wbarnew');
Lnew = (2*Lorinew);    

num_tol = 1e-6;
[Ubext, Ebext] = eig(Lnew);
Ub2ext = Ubext*pinv(SQRT(Ebext),num_tol);
[Uaext,Eaext] = eig(0.5*((Ub2ext'*Wnew*Ub2ext).'+Ub2ext'*Wnew*Ub2ext));
Vext = Ub2ext*Uaext;
Veffext = Vext(:,abs(diag(Eaext)) > num_tol);
Eaext = diag(Eaext);


label2 = kmeans(Veffext(:,end-numC + 1:end),numC,'Replicates',10);% base line

time(2,simu) = toc;

Wnew = sparse(Wnew);
Lnew = sparse(Lnew);

dif(1,simu) = sum(abs(Eaext(end-numC + 1:end) - EigV(end-numC + 1:end)));
dif(4,simu) = sum(abs(Eaext(end-numC + 1:end) - EigVi(end-numC + 1:end)));

[Projext,~] = svd(Veffext(:,end-numC + 1:end));
Projext = Projext(:,1:numC)*Projext(:,1:numC)';

[Basis1,~] = svd(GenEigV(:,end-numC + 1:end));
Basis1Proj = Basis1(:,1:numC)*Basis1(:,1:numC)';

[Basis2,~] = svd(GenEigVi(:,end-numC + 1:end));
Basis2Proj = Basis2(:,1:numC)*Basis2(:,1:numC)';

dif(2,simu) = 1 - 1/3*(abs((Veffext(:,end + 1 -3)'*GenEigV(:,end + 1 -3))/(norm(Veffext(:,end + 1 -3))*norm(GenEigV(:,end + 1 -3)))) + abs((Veffext(:,end + 1 -2)'*GenEigV(:,end + 1 -2))/(norm(Veffext(:,end + 1 -2))*norm(GenEigV(:,end + 1 -2))))+ abs((Veffext(:,end + 1 -1)'*GenEigV(:,end + 1 -1))/(norm(Veffext(:,end + 1 -1))*norm(GenEigV(:,end + 1 -1)))));
dif(3,simu) = 1 - 1/3*(abs((Veffext(:,end + 1 -3)'*GenEigVi(:,end + 1 -3))/(norm(Veffext(:,end + 1 -3))*norm(GenEigVi(:,end + 1 -3))))+abs((Veffext(:,end + 1 -2)'*GenEigVi(:,end + 1 -2))/(norm(Veffext(:,end + 1 -2))*norm(GenEigVi(:,end + 1 -2))))+abs((Veffext(:,end + 1 -1)'*GenEigVi(:,end + 1 -1))/(norm(Veffext(:,end + 1 -1))*norm(GenEigVi(:,end + 1 -1)))));
dif(5,simu) = norm(-Projext/trace(Projext) + Basis1Proj/trace(Basis1Proj),'fro');
dif(6,simu) = norm(-Projext/trace(Projext) + Basis2Proj/trace(Basis2Proj),'fro');

label4 = kmeans(GenEigVi(:,end-numC + 1:end),numC,'Replicates',10);% using initial

G = gsp_graph(Worinew);
cost(1,simu) = MinMaxCut(label2,G,Db,Du);
cost(2,simu) = MinMaxCut(label3,G,Db,Du);
cost(3,simu) = MinMaxCut(label4,G,Db,Du);

W = Wnew;
L = Lnew;
Wini = Worinew;

end
%%
figure
plot(dif(1,:),'-.r','DisplayName','AE with update')
hold on
plot(dif(2,:),'-^r','DisplayName','DC with update')
hold on
plot(dif(5,:),'-xr','DisplayName','SD with update')
hold on
plot(dif(4,:),'-.b','DisplayName','AE without update')
hold on
plot(dif(3,:),'-Vb','DisplayName','DC without update')
hold on
plot(dif(6,:),'o-b','DisplayName','SD without update')
hold on
xlabel('#Perturbations');
ylabel('Estimation Error');
grid on
set(gca,'fontname','DejaVuSans','yscale','log');%
figure
plot(cost(1,:),'-^r','DisplayName','by EVD')
hold on
plot(cost(2,:),'-om','DisplayName','by Update')
hold on
plot(cost(3,:),'.-b','DisplayName','no Update')
hold on
xlabel('#Perturbations');
ylabel('Cost function value');
grid on
set(gca,'fontname','DejaVuSans');%,'yscale','log'
figure
plot(time(1,:),'-.','DisplayName','Time Consumption by Update')
hold on
plot(time(2,:),'-x','DisplayName','Time Consumption by EVD')
hold on
xlabel('#Perturbations');
ylabel('Second');
grid on
set(gca,'fontname','DejaVuSans','yscale','log');
%%
% figure
% for i = 1:num_bs_g
%     X(:,i) = distance0(1,i) + unit*cos(t);
%     Y(:,i) = distance0(2,i) + unit*sin(t);
% end
% for i = 1:num_bs_g
%     plot(X(:,i),Y(:,i),'k')
%     hold on
% end
% for i = 1:1200
%     if label3(i+37) == 1
%     scatter(Newuser_coord(1,i),Newuser_coord(2,i),'r.');
%     hold on
%     end
% if label3(i+37) == 2
%     scatter(Newuser_coord(1,i),Newuser_coord(2,i),'g.');
%     hold on
% end
%     if label3(i+37) == 3
%     scatter(Newuser_coord(1,i),Newuser_coord(2,i),'b.');
%     hold on
%     end
% end
% set(gca,'fontname','DejaVuSans');
%%
% figure
% for i = 1:num_bs_g
%     X(:,i) = distance0(1,i) + unit*cos(t);
%     Y(:,i) = distance0(2,i) + unit*sin(t);
% end
% for i = 1:num_bs_g
%     plot(X(:,i),Y(:,i),'k')
%     hold on
% end
% for i = 1:1200
%     if label(i+37) == 1
%     scatter(Newuser_coord(1,i),Newuser_coord(2,i),'r.');
%     hold on
%     end
% if label(i+37) == 2
%     scatter(Newuser_coord(1,i),Newuser_coord(2,i),'g.');
%     hold on
% end
%     if label(i+37) == 3
%     scatter(Newuser_coord(1,i),Newuser_coord(2,i),'b.');
%     hold on
%     end
% end
% set(gca,'fontname','DejaVuSans');