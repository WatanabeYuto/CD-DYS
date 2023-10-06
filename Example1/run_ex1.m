%

% params_ex1

[x,y,error,fval,constraints_violate] = CD_DYS(cliques,D, DD, A, b, lambda, n, d, G, x0, maxiter, x_opt, f_opt, true);

[x_cpgd_dim,y,error_cpgd_dim,fval_cpgd_dim,constraints_violate_cpgd_dim] = CPGD(cliques,D, DD, A, b, 0, n, d, G, x0, maxiter, x_opt, f_opt, false);
[x_cpgd,y,error_cpgd,fval_cpgd,constraints_violate_cpgd] = CPGD(cliques, D, DD, A, b, 0, n, d, G, x0, maxiter, x_opt, f_opt, true);
[x_acpgd,y,error_acpgd,fval_acpgd,constraints_violate_acpgd] = ACPGD(cliques, D, DD, A, b, 0, n, d, G, x0, maxiter, x_opt, f_opt, true);

[x_extra,error_extra,fval_extra,constraints_violate_extra] = PG_EXTRA(mixing_matrix,lambda,x0,A,b,n,d,x_opt,f_opt,maxiter);
[x_dgd,error_dgd,fval_dgd,constraints_violate_dgd] = DGD(mixing_matrix,lambda,x0,A,b,n,d,x_opt,f_opt,maxiter,false);


plot_residual_ex1
plot_error_ex1