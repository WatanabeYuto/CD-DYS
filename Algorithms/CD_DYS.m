function [x,y,error,res] = CD_DYS(cliques, D, DD, obj_params, lambda, n, d, G, x0, maxiter, x_opt, f_opt)

    %% algorithmic params
    % gamma_i = 1;
    % phi_i = 1;
    % L_fi = zeros(n,1);
    % for i = 1:n
    %     L_fi(i,1) = max(eig(A((i-1)*d+1:i*d,(i-1)*d+1:i*d)'*A((i-1)*d+1:i*d,(i-1)*d+1:i*d)'));
    % end

    alpha = 2/max(eig(obj_params.A'*obj_params.A))*0.99;

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

    for kk = 1:maxiter-1
        DTy = zeros(d*n,1);

        for l = 1:length(cliques) 
            DTy = DTy + kron(DD{l}',eye(d)) * y{l}(:,kk);
        end

        for i = 1:n
            tmp = 1/Qi(i) *  DTy(d*(i-1)+1:d*i,:);
            x(:,i,kk) = prox_l1(tmp, lambda*alpha/Qi(i));
        end
        
        for l = 1:length(cliques) 
            %% y^{k+1/2}
            x_Cl = kron(DD{l},eye(d)) * reshape(x(:,:,kk),[],1);

            %% y^{k+1}
            %% projection onto D_l
            tmp1 = 2*x_Cl - y{l}(:,kk) - kron(DD{l} * inv(D'*D), eye(d)) * alpha * grad_quad(obj_params.A,obj_params.b,reshape(x(:,:,kk),[],1));
            tmp2 = 1/length(cliques{l}) * kron(ones(1,length(cliques{l})),eye(d)) * tmp1;

            %% z^{k+1}
            y{l}(:,kk+1) = y{l}(:,kk) - x_Cl + kron(ones(length(cliques{l}),1),tmp2);
        end

        res(kk,1) = obj_quad(obj_params.A,obj_params.b,reshape(x(:,:,kk),[],1))+ lambda * norm(reshape(x(:,:,kk),[],1),1);
        res(kk,1) = abs(res(kk,1)-f_opt)/f_opt;

        for i = 1:n
            error(kk,1) = error(kk,1) + norm(x_opt - x(:,i,kk))^2;
        end
        error(kk,1) = sqrt(error(kk,1))/norm(kron(ones(n,1),x_opt));
        
        tmp=["k:", kk, "Objective residual:", res(kk,1), "Error:", error(kk,1)];
        disp(tmp)
    end
end
