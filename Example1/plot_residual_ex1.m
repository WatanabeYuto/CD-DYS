%% plot
figure;

len = length(res(1:end));

loglog([1:len],res(1:len))
hold on

loglog([1:len],res_cpgd_dim(1:len))

loglog([1:len],res_cpgd(1:len))
loglog([1:len],res_acpgd(1:len))

loglog([1:len],res_extra(1:len))
loglog([1:len],res_dgd(1:len))

% hold off

xlim([1 len])
ylim([1/10^10 10^3])
xlabel('Iterations')
ylabel('Relative objective residual')

legend({"(a) CD-DYS","(b) CPGD ($\lambda^k = 1/\sqrt{k+1}$)","(c) CPGD ($\lambda^k = \alpha$)","(d) ACPGD ($\lambda^k = \alpha$)","(e) EXTRA","(f) DGD"},'Interpreter', 'latex')

pbaspect([2 1 1])
fontsize(15,"points")
