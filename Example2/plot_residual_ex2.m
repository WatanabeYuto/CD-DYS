%% plot
figure;

len = length(fval(1:end));

semilogy([1:len],fval(1:len))
hold on
semilogy([1:len],fval_e(1:len))
semilogy([1:len],fval_extra(1:len))
semilogy([1:len],fval_fl(1:len))

% hold off
xlim([1 len])
ylim([1/10^11 10^7])
xlabel('Iterations')
ylabel('Relative objective residual')

legend({"(a) CD-DYS ($\mathcal{Q}_\mathcal{G}=\mathcal{Q}_\mathcal{G}^\mathrm{max}$)","(b) CD-DYS ($\mathcal{Q}_\mathcal{G}=\mathcal{G}^\mathrm{edge}$)","(c) PG-EXTRA","(d) CL-FLiP-ADMM ($\mathcal{Q}_\mathcal{G}=\mathcal{Q}_\mathcal{G}^\mathrm{max}$)"},'Interpreter', 'latex')