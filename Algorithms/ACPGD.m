function [x,y,error,fval,constraints_violate] = ACPGD(cliques,W, WW, A, b, lambda, n, d, G, x0, maxiter, x_opt,f_opt, para_flag)

    %% algorithmic params
    % gamma_i = 1;
    % phi_i = 1;
    L_fi = zeros(n,1);
    for i = 1:n
        L_fi(i,1) = max(eig(A(:,:,i)'*A(:,:,i)));
    end
    % alpha = 2/max(L_fi)*0.99;
    alpha = 0.01; 

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
    % A_diag = zeros(size)
    A_diag = zeros(n*d,n*d);
    b_vec = zeros(n*d,1);
    for i = 1:n
        A_diag((i-1)*d+1:i*d,(i-1)*d+1:i*d) = A(:,:,i);
        b_vec((i-1)*d+1:i*d,1) = b(:,i);
    end

    error = zeros(maxiter,1);
    fval = zeros(maxiter,1);
    constraints_violate = zeros(maxiter,1);

    x_k =reshape(x(:,:,1),[],1);
    sigma_k = 0;

    for kk = 1:maxiter-1
        % WTy = zeros(d*n,1);
        % WTu = zeros(d*n,1);
        % for l = 1:length(cliques) 
        %     WTy = WTy + kron(WW{l}',eye(d)) * y{l}(:,kk);
        %     % WTu = WTu + kron(WW{l}',eye(d)) * u{l}(:,kk);
        % end
        if para_flag
            x_plus = reshape(x(:,:,kk),[],1) - alpha *  A_diag'*( A_diag*reshape(x(:,:,kk),[],1) - b_vec );
        else
            x_plus = reshape(x(:,:,kk),[],1) - 1/(kk) *  A_diag'*( A_diag*reshape(x(:,:,kk),[],1) - b_vec );
        end
        % x_plus = reshape(x(:,:,kk),[],1) - alpha *  A_diag'*( A_diag*reshape(x(:,:,kk),[],1) - b_vec );
        
        % for i = 1:n
        %     tmp = 1/Qi(i) *  WTy(d*(i-1)+1:d*i,:);
        %     x(:,i,kk) = prox_l1(tmp, lambda*alpha/Qi(i));
        % end
        x_k_minus = x_k;
        
        for tt = 1:10
        for l = 1:length(cliques) 
            x_Cl = kron(WW{l},eye(d)) * x_plus;

            % tmp = 1/length(cliques{l}) * kron(ones(1,length(cliques{l})),eye(d)) * x_Cl;                            
            q_l = WW{l}*inv(W'*W)*ones(n,1);

            tmp = 1/sum(q_l) * kron(q_l',eye(d)) * x_Cl;

            y{l}(:,kk+1) =  kron(ones(length(cliques{l}),1),tmp);
        end

        WTy = zeros(d*n,1);
        for l = 1:length(cliques) 
            WTy = WTy + kron(WW{l}',eye(d)) * y{l}(:,kk+1);
        end
       
        x_plus = kron(inv(W'*W),eye(d)) * WTy;
        end 

        x_k = x_plus;
        
        sigma_k_minus = sigma_k;
        sigma_k = (1+sqrt(1+4*sigma_k_minus^2))/2;
        tmp = x_k + (sigma_k_minus-1)/sigma_k * (x_k-x_k_minus);

        for i = 1:n
            x(:,i,kk+1) = tmp(d*(i-1)+1:d*i,1);
        end


        constraints_violate(kk,1) = sqrt(constraints_violate(kk,1));

        for i = 1:n
            fval(kk,1) = fval(kk,1) + 0.5 * norm( A(:,:,i) * x(:,i,kk) - b(:,i) )^2 + lambda * norm(x(:,i,kk),1);
        end

        fval(kk,1) = abs(fval(kk,1)-f_opt)/f_opt;

        for i = 1:n
            error(kk,1) = error(kk,1) + norm(x_opt - x(:,i,kk))^2;
        end

        error(kk,1) = sqrt(error(kk,1))/norm(kron(ones(n,1),x_opt));
        [kk, error(kk,1), fval(kk,1),constraints_violate(kk,1)]
    end
end

function vv = prox_l1(ww, ll)
    I = abs(ww) > ll; %% logical 配列を生成
    vv = zeros(size(ww));

    vv(I) = ww(I) - ll* sign(ww(I));
end

function vv = prox_quad(vec, inv, AA, bb, eta)
    vv = inv * ( AA'*bb + 1/eta * vec );
    % vv = ( A'*A + 1/eta *eye(length(b)) ) \ ( A'*b + 1/eta * vec ); 
end