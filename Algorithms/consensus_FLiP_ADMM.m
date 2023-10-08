function [x,y,error,res] = consensus_FLiP_ADMM(cliques,W, WW, A, b, lambda, n, d, G, x0, maxiter, x_opt,f_opt, para_flag)

    %% algorithmic params
    Qi = sum(W,1);
    
    for i = 1:n
        L_fi(i,1) = max(eig(A((i-1)*d+1:i*d,(i-1)*d+1:i*d)'*A((i-1)*d+1:i*d,(i-1)*d+1:i*d)));
    end

    if para_flag
        %% distributed algoritmic parameters
        gamma_i = 1;
        phi_i = 1;
        alpha = zeros(n,1);
        for i = 1:n
            alpha(i,1) = 1/(gamma_i * Qi(i) + L_fi(i,1));
        end
        beta = 1/(gamma_i);
    else
        %% common algoritmic parameters
        gamma_i = 1;
        phi_i = 1;
        alpha = zeros(n,1);
        for i = 1:n
            alpha(i,1) = 1/(gamma_i * Qi(i) + L_fi(i,1));
        end
        alpha = min(alpha) * ones(n,1);
        beta = 1/(gamma_i);
    end

    inv_for_prox_quad = zeros(d,d,n);
    for i = 1:n
        inv_for_prox_quad(:,:,i) = inv( A((i-1)*d+1:i*d,(i-1)*d+1:i*d)'*A((i-1)*d+1:i*d,(i-1)*d+1:i*d) + 1/alpha(i,1) * eye(d));
    end

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
        WTu = zeros(d*n,1);
        for l = 1:length(cliques) 
            WTy = WTy + kron(WW{l}',eye(d)) * y{l}(:,kk);
            WTu = WTu + kron(WW{l}',eye(d)) * u{l}(:,kk);
        end

        tmp = reshape(x(:,:,kk),[],1) - kron(diag(alpha),eye(d)) * (grad_quad(A,b,reshape(x(:,:,kk),[],1)) + WTu + gamma_i * (kron(diag(Qi),eye(d)) * reshape(x(:,:,kk),[],1) - WTy));
        x(:,:,kk+1) = reshape(tmp,[d,n]);

        for l = 1:length(cliques)
            x_Cl = kron(WW{l},eye(d)) * reshape(x(:,:,kk),[],1);

            %% projection
            tmp_yl = 1/length(cliques{l}) * kron(ones(length(cliques{l}),1),eye(d))' *  (y{l}(:,kk) - beta * ( - u{l}(:,kk) - gamma_i * ( x_Cl - y{l}(:,kk)) ) );
            
            sum_prox_const = 0;
            for j = 1:length(cliques{l})
                ii = cliques{l}(j);
                sum_prox_const = sum_prox_const + 1/Qi(ii);
            end

            y{l}(:,kk+1) = kron(ones(length(cliques{l}),1), prox_l1( tmp_yl, lambda / length(cliques{l}) * sum_prox_const * beta ));
            u{l}(:,kk+1) = u{l}(:,kk) + phi_i* gamma_i * (x_Cl - y{l}(:,kk+1)); 
        end

        res(kk,1) = obj_quad(A,b,reshape(x(:,:,kk),[],1))+ lambda * norm(reshape(x(:,:,kk),[],1),1);
        res(kk,1) = abs(res(kk,1)-f_opt)/f_opt;


        for i = 1:n
            error(kk,1) = error(kk,1) + norm(x_opt - x(:,i,kk))^2;
        end

        error(kk,1) = sqrt(error(kk,1))/norm(kron(ones(n,1),x_opt));

        tmp=["k:", kk, "Objective residual:", res(kk,1), "Error:", error(kk,1)];
        disp(tmp)
    end
end
