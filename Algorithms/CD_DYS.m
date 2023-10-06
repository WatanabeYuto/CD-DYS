function [x,y,error,fval,constraints_violate] = CD_DYS(cliques,W, WW, A, b, lambda, n, d, G, x0, maxiter, x_opt,f_opt, para_flag)

    %% algorithmic params
    % gamma_i = 1;
    % phi_i = 1;
    L_fi = zeros(n,1);
    for i = 1:n
        L_fi(i,1) = max(eig(A(:,:,i)'*A(:,:,i)));
    end
    % alpha = 2/max(L_fi)*0.99;

    % alpha = diag(2/L_fi)*0.99;

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

    for kk = 1:maxiter-1
        WTy = zeros(d*n,1);
        % WTu = zeros(d*n,1);
        for l = 1:length(cliques) 
            WTy = WTy + kron(WW{l}',eye(d)) * y{l}(:,kk);
            % WTu = WTu + kron(WW{l}',eye(d)) * u{l}(:,kk);
        end


        for i = 1:n
            tmp = 1/Qi(i) *  WTy(d*(i-1)+1:d*i,:);
            % x(:,i,kk) = prox_l1(tmp, lambda*alpha/Qi(i));

            x(:,i,kk+1) = prox_l1(tmp, lambda*alpha(i,i)/Qi(i));
        end

        for l = 1:length(cliques) 
            x_Cl = kron(WW{l},eye(d)) * reshape(x(:,:,kk+1),[],1);

            tmpp = 2*x_Cl - y{l}(:,kk) - kron(WW{l} * inv(W'*W) , eye(d)) * kron(alpha,eye(d)) * A_diag'*( A_diag*reshape(x(:,:,kk+1),[],1) - b_vec );
            tmp = 1/length(cliques{l}) * kron(ones(1,length(cliques{l})),eye(d)) * tmpp;
            y{l}(:,kk+1) = y{l}(:,kk) - x_Cl + kron(ones(length(cliques{l}),1),tmp);
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