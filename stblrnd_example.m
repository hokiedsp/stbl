clear; close all; drawnow;

alpha = [.5 .75 1.5 2];
beta = 0;
gam = 1;
delta = 0;

X = stbl.rnd(alpha,beta,gam,delta,10000,1);

hg_args = {201,'BinLimits',[-10,10],'Normalization','probability','DisplayStyle','stairs'};
histogram(X(:,1),hg_args{:},'DisplayName',sprintf('\\alpha = %g',alpha(1)))
hold on
for n = 2:size(X,2)
   histogram(X(:,n),hg_args{:},'DisplayName',sprintf('\\alpha = %g',alpha(n)));
end
hold off
title(sprintf('\\alpha-stable random numbers with \\beta = %g, \\gamma = %g, and \\delta = %g',beta,gam,delta))
legend show