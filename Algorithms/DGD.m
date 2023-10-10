function [x,error,res] = DGD(mixing_matrix,lambda,x0,obj_params,n,d,x_opt,f_opt,maxiter,stepsize_flag)

    x = zeros(n*d,maxiter);
    w = zeros(n*d,maxiter);
    
    x(:,1) = reshape(x0,[],1);

    error = zeros(maxiter,1);
    res = zeros(maxiter,1);

    % alpha = 0.99 * (1 + min(eig(mixing_matrix))) / max(eig(obj_params.A'*obj_params.A)) ;
    alpha=0.01;

    for kk = 1:maxiter-1

        if stepsize_flag
            x(:,kk+1) = kron(mixing_matrix, eye(d)) * x(:,kk) - alpha * grad_quad(obj_params.A,obj_params.b,x(:,kk));
        else
            x(:,kk+1) = kron(mixing_matrix, eye(d)) * x(:,kk) - 1/sqrt(kk) * grad_quad(obj_params.A,obj_params.b,x(:,kk)); 
        end
        % x(:,kk+1) = prox_l1( vec , lambda * alpha );
        % w(:,kk+1) = w(:,kk) + 0.5 * kron((eye(n)-mixing_matrix), eye(d)) * x(:,kk);

        res(kk,1) =  obj_quad(obj_params.A,obj_params.b,x(:,kk)) + lambda * norm(x(:,kk),1);

        res(kk,1) = abs(res(kk,1)-f_opt)/f_opt;

        for i = 1:n
            error(kk,1) = error(kk,1) + norm(x_opt - x(d*(i-1)+1:d*i,kk))^2;
        end
        error(kk,1) = sqrt(error(kk,1))/norm(kron(ones(n,1),x_opt));

        tmp=["k:", kk, "Objective residual:", res(kk,1), "Error:", error(kk,1)];
        disp(tmp)
        
    end
end