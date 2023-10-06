function [x,y,error,fval,constraints_violate] = consensus_FLiP_ADMM(cliques,W, WW, A, b, lambda, n, d, G, x0, maxiter, x_opt,f_opt, para_flag)

    %% algorithmic params
    % gamma_i = 1;
    % phi_i = 1;
    % alpha = zeros(n,1);
    inv_for_prox_quad = zeros(d,d,n);
    Qi = sum(W,1);
    
    for i = 1:n
        L_fi(i,1) = max(eig(A(:,:,i)'*A(:,:,i)));
    end

    if para_flag
        gamma_i = 1;
        phi_i = 1;
        alpha = zeros(n,1);
        for i = 1:n
            alpha(i,1) = 1/(gamma_i * Qi(i) + L_fi(i,1));
        end
        beta = 1/(gamma_i);
    else
        gamma_i = 1;
        phi_i = 1;
        alpha = zeros(n,1);
        for i = 1:n
            alpha(i,1) = 1/(gamma_i * Qi(i) + L_fi(i,1));
        end
        alpha = min(alpha) * ones(n,1);
        beta = 1/(gamma_i);
    end

    for i = 1:n
        % alpha(i,1) = 1/(gamma_i * Qi(i));
        inv_for_prox_quad(:,:,i) = inv( A(:,:,i)'*A(:,:,i) + 1/alpha(i,1) * eye(d));
    end

    % beta = zeros(n,1);
    % for i = 1:n
    %     beta(i,1) = 1/(gamma_i);
    % end

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
    fval = zeros(maxiter,1);
    constraints_violate = zeros(maxiter,1);

    for kk = 1:maxiter-1
        WTy = zeros(d*n,1);
        WTu = zeros(d*n,1);
        for l = 1:length(cliques) 
            WTy = WTy + kron(WW{l}',eye(d)) * y{l}(:,kk);
            WTu = WTu + kron(WW{l}',eye(d)) * u{l}(:,kk);
        end

        for i = 1:n
            x(:,i,kk+1) = x(:,i,kk) - alpha(i,1) * ( A(:,:,i)'*( A(:,:,i)*x(:,i,kk) - b(:,i) ) + WTu(d*(i-1)+1:d*i,1) + gamma_i * ( Qi(i) *  x(:,i,kk) - WTy(d*(i-1)+1:d*i,1)) );
            % x(:,i,kk+1) = prox_quad( vec, inv_for_prox_quad(:,:,i), A(:,:,i), b(:,i), alpha(i,1) );
        end

        for l = 1:length(cliques)
            x_Cl = kron(WW{l},eye(d)) * reshape(x(:,:,kk+1),[],1);
            % beta_l = WW{l} * beta;
            % beta_l_inv = inv(diag(WW{l} * beta)) * ones(length(cliques{l}),1);

            tmp_yl = 1/length(cliques{l}) * kron(ones(length(cliques{l}),1),eye(d))' *  (y{l}(:,kk) - beta * ( - u{l}(:,kk) - gamma_i * ( x_Cl - y{l}(:,kk)) ) );
            % tmp_yl = 1/sum(beta_l_inv(:,1)) * kron(beta_l_inv,eye(d))' *  (y{l}(:,kk) - kron(diag(beta_l),eye(d)) * ( - u{l}(:,kk) - gamma_i * ( x_Cl - y{l}(:,kk)) ) );
            
            sum_prox_const = 0;
            for j = 1:length(cliques{l})
                ii = cliques{l}(j);
                sum_prox_const = sum_prox_const + 1/Qi(ii);
            end

            y{l}(:,kk+1) = kron(ones(length(cliques{l}),1), prox_l1( tmp_yl, lambda / length(cliques{l}) * sum_prox_const * beta ));
            % y{l}(:,kk+1) = kron(ones(length(cliques{l}),1), prox_l1( tmp_yl, lambda * sum_prox_const / sum(beta_l_inv)));

            u{l}(:,kk+1) = u{l}(:,kk) + phi_i* gamma_i * (x_Cl - y{l}(:,kk+1)); 

            constraints_violate(kk,1) = constraints_violate(kk,1) + norm(x_Cl - y{l}(:,kk+1))^2;
        end

        constraints_violate(kk,1) = sqrt(constraints_violate(kk,1));

        for i = 1:n
            fval(kk,1) = fval(kk,1) + 0.5 * norm( A(:,:,i) * x(:,i,kk) - b(:,i) )^2 + lambda * norm(x(:,i,kk),1);
        end

        fval(kk,1) = abs(fval(kk,1)-f_opt)/f_opt;

        for i = 1:n
            error(kk,1) = error(kk,1) + norm(x_opt - x(:,i,kk+1))^2;
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