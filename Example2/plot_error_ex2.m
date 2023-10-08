%% plot
figure;

len = length(error(1:end));

semilogy([1:len],error(1:len))
hold on
semilogy([1:len],error_e(1:len))
semilogy([1:len],error_extra(1:len))
semilogy([1:len],error_fl(1:len))

% hold off
xlim([1 len])
ylim([1/10^5 10^4])
xlabel('Iterations')
ylabel({'$\|\mathbf{x}^k-\mathbf{x}^*\|/\|\mathbf{x}^*\|$'},'Interpreter', 'latex')

legend({"(a) CD-DYS ($\mathcal{Q}_\mathcal{G}=\mathcal{Q}_\mathcal{G}^\mathrm{max}$)","(b) CD-DYS ($\mathcal{Q}_\mathcal{G}=\mathcal{Q}_\mathcal{G}^\mathrm{edge}$)","(c) PG-EXTRA","(d) CL-FLiP-ADMM ($\mathcal{Q}_\mathcal{G}=\mathcal{Q}_\mathcal{G}^\mathrm{max}$)"},'Interpreter', 'latex')


pbaspect([2 1 1])
fontsize(15,"points")
