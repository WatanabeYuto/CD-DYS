%

params_ex1

[x,y,error,fval,constraints_violate] = CD_DYS(cliques,D, DD, A, b, lambda, n, d, G, x0, maxiter, x_opt, f_opt);

%% {true,false} = {fixed,diminishing} step size

[x_cpgd_dim,y,error_cpgd_dim,res_cpgd_dim] = CPGD(cliques,D, DD, A, b, 0, n, d, G, x0, maxiter, x_opt, f_opt, false);
[x_cpgd,y,error_cpgd,res_cpgd] = CPGD(cliques, D, DD, A, b, 0, n, d, G, x0, maxiter, x_opt, f_opt, true);
[x_acpgd,y,error_acpgd,res_acpgd] = ACPGD(cliques, D, DD, A, b, 0, n, d, G, x0, maxiter, x_opt, f_opt, true);

[x_extra,error_extra,res_extra] = PG_EXTRA(mixing_matrix,lambda,x0,A,b,n,d,x_opt,f_opt,maxiter);
[x_dgd,error_dgd,res_dgd] = DGD(mixing_matrix,lambda,x0,A,b,n,d,x_opt,f_opt,maxiter,false);


plot_residual_ex1
plot_error_ex1