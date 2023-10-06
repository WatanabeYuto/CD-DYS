function [x,error,fval,constraints_violate] = PG_EXTRA(mixing_matrix,lambda,x0,A,b,n,d,x_opt,f_opt,maxiter,flag)

    x = zeros(n*d,maxiter);
    w = zeros(n*d,maxiter);
    
    x(:,1) = reshape(x0,[],1);

    error = zeros(maxiter,1);
    fval = zeros(maxiter,1);
    constraints_violate = zeros(maxiter,1); 

    AA = zeros(d*n,d*n);
    bb = zeros(d*n,1);

    for i = 1:n
        AA(d*(i-1)+1:d*i,d*(i-1)+1:d*i) = A(:,:,i);
        bb(d*(i-1)+1:d*i,1) = b(:,i);
    end

    alpha = 0.99 * (1 + min(eig(mixing_matrix))) / max(eig(AA'*AA)) ;

    for kk = 1:maxiter-1
        if flag
            x(:,kk+1) = kron(mixing_matrix, eye(d)) * x(:,kk) - 1/sqrt(kk) * AA'*(AA*x(:,kk) - bb); 
        else
            alpha=0.01;
            x(:,kk+1) = kron(mixing_matrix, eye(d)) * x(:,kk) - alpha * AA'*(AA*x(:,kk) - bb);
        end
        % x(:,kk+1) = prox_l1( vec , lambda * alpha );
        % w(:,kk+1) = w(:,kk) + 0.5 * kron((eye(n)-mixing_matrix), eye(d)) * x(:,kk);

        fval(kk,1) =  0.5 * norm( AA(:,:) * x(:,kk) - bb(:,1) )^2 + lambda * norm(x(:,kk),1);

        fval(kk,1) = abs(fval(kk,1)-f_opt)/f_opt;

        for i = 1:n
            error(kk,1) = error(kk,1) + norm(x_opt - x(d*(i-1)+1:d*i,kk))^2;
        end

        % constraints_violate(kk,1) = 1/n * norm(  )
        % error(kk,1) = sqrt(error(kk,1))/n
        error(kk,1) = sqrt(error(kk,1))/norm(kron(ones(n,1),x_opt));
        [kk, error(kk,1), fval(kk,1)]
    end
end

function vv = prox_l1(ww, ll)
    I = abs(ww) > ll; %% logical 配列を生成
    vv = zeros(size(ww));

    vv(I) = ww(I) - ll* sign(ww(I));
end