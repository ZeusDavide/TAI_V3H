function [H,obj] = V3H(X,num_view,W,Z_ini,F_ini,num_cluster,alpha,beta,gamma,miu,rho,max_iter,num_instance)
% Code for revision "V3H: View Variation and View Heredity for Incomplete Multi-view Clustering"
% writen by Xiang Fang (xfang9508@gmail.com)
%%
%Input
%X: incomplete multi-view data with num_view views
%W:  the incomplete index matrix
%num_cluster: cluster number
%num_view: view number
%alpha,beta,gamma: hyper-parameter
%Z_ini:Initialized subspace matrix
%F_ini:Initialized cluster indicator matrix
%miu,rho,max_iter:Intermediate parameters
%num_instance:instance number
%%
%Output
%H:  the consensus cluster indicator matrix
%obj:objective function value
%%
Z = Z_ini;
M = zeros(num_instance);
p=10;
F = F_ini;
epislon=1e0;
w=1;
seta_E=100;
for i = 1:num_view
    E{i}  = zeros(size(X{i}));
    C1{i} = zeros(size(X{i}));
    C2{i} = zeros(num_instance);
    N{i} = zeros(num_instance);
end
for iter = 1:max_iter
    XXv = 0;
    % ------------Update H -------------- %
    FFF = 0;
    for i = 1:num_view
        FFF = FFF+F{i}*F{i}';
    end
    FFF(isnan(FFF))=0;   
    FFF(isinf(FFF))=1e5;  
    try
        [V,D] = eig(FFF);   
    catch ME
        if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))   
            [V,D] = eig(FFF, eye(size(FFF)));
        else
            rethrow(ME); 
        end
    end  
    [~, ind] = sort(diag(D),'descend');
    H = V(:, ind(1:num_cluster));   
        
    FFF = 0; 
    linshi_obj1 = 0;
    linshi_obj2 = 0;    
    for i = 1:num_view
        % --------------Update Z{i} ------------ % 
        G1 = X{i}'*X{i};
        G2 =W{i}'*W{i};
        G3 = X{i}'*(X{i}-E{i}+C1{i}/miu)+W{i}*(M-p*N{i}-C2{i}/miu)*W{i}';
%         G1 = X{i}-E{i}+C1{i}/miu;
%         G2 = M-N{i}-C2{i}/miu;
%         Z{i} = (X{i}'*X{i}+G{i}'*G{i})\(X{i}'*G1+G{i}*G2*G{i}');
        Z{i}=lyap(G1,G2,G3);
        clear G1 G2 G3
        % ------------Update N{i}-------------- %
        P = F{i};
        Q = L2_distance_1(P',P');
        C = M-W{i}'*Z{i}*W{i}-C2{i}/miu;
        linshi_W = C/p-0.5*alpha*Q/miu;
        linshi_W = linshi_W-diag(diag(linshi_W));
        for ic = 1:size(Z{i},2)
            ind = 1:size(Z{i},2);
            ind(ic) = [];
            N{i}(ic,ind) = EProjSimplex_new(linshi_W(ic,ind));
        end
        clear linshi_W P Q M ind ic
        % ---------------- Update M ---------------- %
        temp = p*N{i}+W{i}'*Z{i}*W{i}+C2{i}/miu;
%         [AU,SU,VU] = svd(temp,'econ'); 
        [AU,SU,VU] = svd(temp,'econ'); 
         AU(isnan(AU)) = 0;
         VU(isnan(VU)) = 0;
         SU(isnan(SU)) = 0;
        T0=zeros(num_instance,1);
%         for t = 1:10
            lambda=0.5/rho;
            S0 = diag(SU);
            grad=(w+epislon)*epislon./(epislon+w*S0).^2;
            T1=max(S0-lambda*grad,0);
         M=AU*diag(T1)*VU';
        clear temp AU SU VU;
        % -------------- Update F{i} --------- %
        WW = p*(abs(N{i})+abs(N{i}'))*0.5;
        LL = diag(sum(WW))-WW;
        C = gamma*(H*H')-alpha*W{i}'*LL*W{i};
        C(isnan(C))=0;
        C(isinf(C))=1e5;
        try
            [V,D] = eig(C);
        catch ME
            if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
                [V,D] = eig(C, eye(size(C)));
            else
                rethrow(ME);
            end
        end  
        [~, ind] = sort(diag(D),'descend');
        F{i} = V(:, ind(1:num_cluster));   
        clear V D WW D_sort ind
        FFF = FFF+F{i}*F{i}';
        % -------------5 E{i}     shrinkage-------- %
        for k = 1:size(E{i},1)
            tmp(k)=norm(E{i}(k,:),2);
            D_E{i}(k,k)=(seta_E+1)*(seta_E*tmp(k)+2)/(seta_E*tmp(k)+1)^2;
        end
        temp1 = X{i}-X{i}*Z{i}+C1{i}/miu;
        invD = diag(1./diag(D_E{i}));
        E{i} = (beta/miu* eye(size(D_E{i}))+invD)\temp1;      
        clear temp1
        % ----------- C1 C2 C3 --------- %
        leq1 = X{i}-X{i}*Z{i}-E{i};
        leq2 = W{i}'*Z{i}*W{i}+p*N{i}-M;
        C1{i} = C1{i} + miu*leq1;
        C2{i} = C2{i} + miu*leq2;
        % ---------- obj --------- %
        linshi_obj1 = linshi_obj1+sum(abs(S0))+alpha*trace(F{i}'*W{i}'*LL*W{i}*F{i})+beta*sum(1./diag(D_E{i}));
        linshi_obj2 = linshi_obj2+norm(leq1,'fro')^2+norm(leq2,'fro')^2;  
        XXv = XXv + norm(X{i},'fro')^2;
    end
    % ----------------- obj ----------- %
    obj(iter) = (linshi_obj1+linshi_obj2+gamma*(num_view*num_cluster-trace(H'*FFF*H)))/XXv;
    clear FFF
    % ---------- miu ------------- %
    miu = min(rho*miu,1e8);
    if iter > 2 && abs(obj(iter)-obj(iter-1))<1e-7
        iter
        break;
    end
end

end
