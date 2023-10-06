%% plot
figure;

len = length(fval(1:end));

loglog([1:len],fval(1:len))
hold on

loglog([1:len],fval_cpgd_dim(1:len))

loglog([1:len],fval_cpgd(1:len))
loglog([1:len],fval_acpgd(1:len))

loglog([1:len],fval_extra(1:len))
loglog([1:len],fval_dgd(1:len))

% hold off

xlim([1 len])
ylim([1/10^10 10^3])
xlabel('Iterations')
ylabel('Relative objective residual')

legend({"(a) CD-DYS","(b) CPGD ($\lambda^k = 1/\sqrt{k+1}$)","(c) CPGD ($\lambda^k = \alpha$)","(d) ACPGD ($\lambda^k = \alpha$)","(e) EXTRA","(f) DGD"},'Interpreter', 'latex')
