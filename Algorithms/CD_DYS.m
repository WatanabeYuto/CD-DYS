function [x,y,error,res] = CD_DYS(cliques,W, WW, A, b, lambda, n, d, G, x0, maxiter, x_opt,f_opt)

    %% algorithmic params
    % gamma_i = 1;
    % phi_i = 1;
    L_fi = zeros(n,1);
    for i = 1:n
        L_fi(i,1) = max(eig(A(:,:,i)'*A(:,:,i)));
    end

    alpha = eye(n)*2/max(L_fi)*0.99;

    Qi = sum(W,1);
    
    x = zeros(d,n,maxiter);
    x(:,:,1) = x0;

    y = {};
    u = {};
    for l = 1:length(cliques)
        y{l} = zeros(length(cliques{l})*d, maxiter);

        y{l}(:,1) = kron(WW{l},eye(d)) * reshape(x(:,:,1),[],1); 

        u{l} = zeros(length(cliques{l})*d, maxiter); 
    end


    error = zeros(maxiter,1);
    res = zeros(maxiter,1);

    for kk = 1:maxiter-1
        WTy = zeros(d*n,1);

        for l = 1:length(cliques) 
            WTy = WTy + kron(WW{l}',eye(d)) * y{l}(:,kk);
        end

        for i = 1:n
            tmp = 1/Qi(i) *  WTy(d*(i-1)+1:d*i,:);
            x(:,i,kk+1) = prox_l1(tmp, lambda*alpha(i,i)/Qi(i));
        end
        
        for l = 1:length(cliques) 
            %% y^{k+1/2}
            x_Cl = kron(WW{l},eye(d)) * reshape(x(:,:,kk+1),[],1);

            %% y^{k+1}
            %% projection onto D_l
            tmpp = 2*x_Cl - y{l}(:,kk) - kron(WW{l} * inv(W'*W) , eye(d)) * kron(alpha,eye(d)) * grad_quad(A,b,reshape(x(:,:,kk+1),[],1));
            tmp = 1/length(cliques{l}) * kron(ones(1,length(cliques{l})),eye(d)) * tmpp;

            %% z^{k+1}
            y{l}(:,kk+1) = y{l}(:,kk) - x_Cl + kron(ones(length(cliques{l}),1),tmp);
        end

        res(kk,1) = obj_quad(A,b,reshape(x(:,:,kk),[],1));
        res(kk,1) = abs(res(kk,1)-f_opt)/f_opt;

        for i = 1:n
            error(kk,1) = error(kk,1) + norm(x_opt - x(:,i,kk))^2;
        end
        error(kk,1) = sqrt(error(kk,1))/norm(kron(ones(n,1),x_opt));
        
        tmp=["k:", kk, "Objective residual:", res(kk,1), "Error:", error(kk,1)];
        disp(tmp)
    end
end
