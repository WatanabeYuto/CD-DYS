%% plot
figure;

len = length(error(1:end));

loglog([1:len],error(1:len))
hold on

loglog([1:len],error_cpgd_dim(1:len))

loglog([1:len],error_cpgd(1:len))
loglog([1:len],error_acpgd(1:len))

loglog([1:len],error_extra(1:len))
loglog([1:len],error_dgd(1:len))

% hold off
xlim([1 len])
ylim([1/10^5 10^2])
xlabel('Iterations')
ylabel({'$\|\mathbf{x}^k-\mathbf{x}^*\|/\|\mathbf{x}^*\|$'},'Interpreter', 'latex')

legend({"(a) CD-DYS","(b) CPGD ($\lambda^k = 1/\sqrt{k+1}$)","(c) CPGD ($\lambda^k = \alpha$)","(d) ACPGD ($\lambda^k = \alpha$)","(e) EXTRA","(f) DGD"},'Interpreter', 'latex')
