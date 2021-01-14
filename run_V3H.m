%NOTE THAT our pape "V3H: View Variation and View Heredity for Incomplete Multi-view Clustering"
% performs the experiments in MATLAB R2020a and our codes run on a Windows 10 machine with 3:30 GHz E3-1225 CPU, 64 GB main memory.
clear;
clc
%hyper-parameter setting
alpha  = 1e-3;
beta  = 1e-4;
gamma  = 1e-1;
% load multi-view dataset
Dataname = 'YaleB';
load(Dataname);
% truthF: real label
truthF = y;  
num_cluster = length(unique(truthF));
num_instance=length(truthF);
num_view = length(X);
%missing rate
missing_rate=0.2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ;

W = cell(num_view,1);
ind = ones(num_instance, num_view);
ind = splitDigitData(ind, missing_rate, 1 );
for i = 1:num_view
    item = ind(:,i);
    W{i} = diag(item);
    temp = item == 0; 
    X{i}(temp,:)= 0;
    X{i} = mapminmax(X{i},0,1);
    X{i} = NormalizeFea(X{i},1);
end

    
for iv = 1:num_view
    options = [];
    options.NeighborMode = 'KNN';
    options.k = 5;
    options.WeightMode = 'Binary';
    Z1 = constructW(X{iv}',options);
    Z_ini{iv} = full(Z1);  
    clear Z1;
end

F_ini = solveF(Z_ini,W,num_cluster);

max_iter = 30;
miu = 0.01;
rho = 1.1;

[H,obj] = V3H_R1(X,num_view,W,Z_ini,F_ini,num_cluster,alpha,beta,gamma,miu,rho,max_iter,num_instance);
norm_mat = repmat(sqrt(sum(H.*H,2)),1,size(H,2));
for i = 1:size(norm_mat,1)
    if (norm_mat(i,1)==0)
        norm_mat(i,:) = 1;
    end
end
H = H./norm_mat; 

repeat = 10;
for iter_num = 1:repeat
    % pre_labels: clustering result label
    pre_labels    = kmeans(real(H),num_cluster,'emptyaction','singleton','replicates',20,'display','off');
    result = ClusteringMeasure(truthF, pre_labels);       
    ACC(iter_num)    = result(1)*100;
    NMI(iter_num) = result(2)*100;
    Purity(iter_num)= result(3)*100;
end

%ACC,NMI,and Purity: the clustering evaluation metrics
mean_ACC = mean(ACC)
mean_NMI = mean(NMI)
mean_PUR = mean(Purity)


