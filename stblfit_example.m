clear; close all; drawnow;

N = 300;
sampsize = 100;
s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',50); % For reproducibility
X = zeros(N,1);
for i = 1:N
   % Generate a normalized sum of Pareto-type random variables
   Samp = 1./rand(s,sampsize,1).^(4/3);
   X(i) = sum(Samp)/sampsize^(4/3); % Normalize sum
end
% estimate parameters
p = stbl.fit(X,'ecf',statset('Display','iter'));
% plot data with fit parameters
xmax = 15;

H = figure(1);
set(H,'Position', [517 626 939 410]);
clf;
title('Stable fit to sums of Pareto random variables');
subplot(1,2,1)
hold on
stem(X(X < xmax),stbl.pdf(X(X<xmax),p(1),p(2),p(3),p(4),'quick'));
x = 0:.1:xmax;
plot(x,stbl.pdf(x,p(1),p(2),p(3),p(4),'quick'),'r-')
hold off
xlabel(['\alpha_0 = ',num2str(p(1)),'  \beta_0 = ',num2str(p(2)),'  \gamma_0 = ',num2str(p(3)),'  \delta_0 = ',num2str(p(4))]);
legend('Data','Fit stable density')
subplot(1,2,2)

CDF = prctile(X,[1:75]);
cmin = CDF(1);
cmax = CDF(end);
x = cmin:.1:cmax;
estCDF = stbl.cdf(x,p(1),p(2),p(3),p(4));

plot(CDF,[.01:.01:.75],'b.',x,estCDF,'r-')
legend('Empirical CDF','Estimated CDF','Location','northwest')
