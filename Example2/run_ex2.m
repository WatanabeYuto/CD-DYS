% params_ex2;

[x,y,error,fval,constraints_violate] = CD_DYS(cliques, D, DD, A, b, lambda, n, d, G, x0, maxiter, x_opt, f_opt, true);
[x_e,y_e,error_e,fval_e,constraints_violate_e] = CD_DYS(edges, D_e, DD_e, A, b, lambda, n, d, G, x0, maxiter, x_opt, f_opt, true);
[x_extra,error_extra,fval_extra,constraints_violate_extra] = PG_EXTRA(mixing_matrix,lambda,x0,A,b,n,d,x_opt,f_opt,maxiter);

[x_fl,y_fl,error_fl,fval_fl,constraints_violate_fl_para] = consensus_FLiP_ADMM(cliques, D, DD, A, b, lambda, n, d, G, x0, maxiter, x_opt, f_opt,true);

plot_residual_ex2
plot_error