% run

params_ex2;

[x,y,error,res] = CD_DYS(cliques, D, DD, obj_params, lambda, n, d, G, x0, maxiter, x_opt, f_opt);
[x_e,y_e,error_e,res_e] = CD_DYS(edges, D_e, DD_e, obj_params, lambda, n, d, G, x0, maxiter, x_opt, f_opt);
[x_extra,error_extra,res_extra] = PG_EXTRA(mixing_matrix,lambda,x0,obj_params,n,d,x_opt,f_opt,maxiter);

[x_fl,y_fl,error_fl,res_fl] = consensus_FLiP_ADMM(cliques, D, DD, obj_params, lambda, n, d, G, x0, maxiter, x_opt, f_opt,true);

plot_residual_ex2
plot_error_ex2