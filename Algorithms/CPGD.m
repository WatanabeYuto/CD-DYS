function [x,y,error,res] = CPGD(cliques,D, DD, obj_params, lambda, n, d, G, x0, maxiter, x_opt,f_opt, stepsize_flag)
    %% algorithmic params
    alpha = 0.01;
    Qi = sum(D,1);    
    x = zeros(d,n,maxiter);
    x(:,:,1) = x0;

    y = {};
    u = {};
    for l = 1:length(cliques)
        y{l} = zeros(length(cliques{l})*d, maxiter);

        y{l}(:,1) = kron(DD{l},eye(d)) * reshape(x(:,:,1),[],1); 

        u{l} = zeros(length(cliques{l})*d, maxiter); 
    end

    error = zeros(maxiter,1);
    res = zeros(maxiter,1);

    p = 10;

    for kk = 1:maxiter-1

        if stepsize_flag
            %% fixed stepsize
            x_plus = reshape(x(:,:,kk),[],1) - alpha *  grad_quad(obj_params.A,obj_params.b,reshape(x(:,:,kk),[],1));
        else
            %% diminishing stepsize
            x_plus = reshape(x(:,:,kk),[],1) - 1/sqrt(kk) * grad_quad(obj_params.A,obj_params.b,reshape(x(:,:,kk),[],1));
        end

        %% clique-based projection
        for tt = 1:p
            for l = 1:length(cliques)
                x_Cl = kron(DD{l},eye(d)) * x_plus;
                
                q_l = DD{l}*inv(D'*D)*ones(n,1);
            
                tmp = 1/sum(q_l) * kron(q_l',eye(d)) * x_Cl;
                y{l}(:,kk+1) =  kron(ones(length(cliques{l}),1),tmp);
            end

            DTy = zeros(d*n,1);
            for l = 1:length(cliques) 
                DTy = DTy + kron(DD{l}',eye(d)) * y{l}(:,kk+1);
            end
            
            x_plus = kron(inv(D'*D),eye(d)) * DTy;
        end

        for i = 1:n
            x(:,i,kk+1) = x_plus(d*(i-1)+1:d*i,1);
        end

        res(kk,1) = obj_quad(obj_params.A,obj_params.b,reshape(x(:,:,kk),[],1)) + lambda * norm(reshape(x(:,:,kk),[],1),1);
        res(kk,1) = abs(res(kk,1)-f_opt)/f_opt;

        for i = 1:n
            error(kk,1) = error(kk,1) + norm(x_opt - x(:,i,kk))^2;
        end

        error(kk,1) = sqrt(error(kk,1))/norm(kron(ones(n,1),x_opt));

        tmp=["k:", kk, "Objective residual:", res(kk,1), "Error:", error(kk,1)];
        disp(tmp)

    end
end
