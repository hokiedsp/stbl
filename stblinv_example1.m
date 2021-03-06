clear; close all; drawnow;

figure(4)
u = .01:.01:.99;
beta = 0;
gam = 1;
delta = 0;
plot(u , stbl.inv(u,.5,beta,gam,delta),...
   u , stbl.inv(u,.75,beta,gam,delta),...
   u , stbl.inv(u,1,beta,gam,delta),...
   u , stbl.inv(u,1.25,beta,gam,delta),...
   u , stbl.inv(u,1.5,beta,gam,delta) )
axis([0 1 -2 2]);
title('Symmetric inverse \alpha-stable CDFs, \beta = 0, \gamma = 1, \delta = 0');
legend('\alpha = 0.5',...
   '\alpha = 0.75',...
   '\alpha = 1.0',...
   '\alpha = 1.25',...
   '\alpha = 1.5',...
   'Location','northwest');


figure(5)
u = .01:.01:.99;
beta = .5;
gam = 1;
delta = 0;
plot( u , stbl.inv(u,.5,beta,gam,delta),...
   u , stbl.inv(u,.75,beta,gam,delta),...
   u , stbl.inv(u,1,beta,gam,delta),...
   u , stbl.inv(u,1.25,beta,gam,delta),...
   u , stbl.inv(u,1.5,beta,gam,delta) )
axis([0 1 -2 2]);
title('Skewed inverse \alpha-stable CDFs, \beta = .5, \gamma = 1, \delta = 0');
legend( '\alpha = 0.5',...
   '\alpha = 0.75',...
   '\alpha = 1.0',...
   '\alpha = 1.25',...
   '\alpha = 1.5',...
   'Location','northwest');
