clear; close all; drawnow

x = 10.^(-2:.25:7);
alpha = 1.5;
beta = 0;
gam = 1;
delta = 0;
C = (1 - alpha)/(gamma(2 - alpha)*cos(pi*alpha/2));
loglog( x , stbl.pdf( x,alpha,beta,gam,delta,'quick'), 'g-*', ...
   x , stbl.pdf( x,alpha,beta,gam,delta),         'b-o', ...
   x , C/2 * alpha * x.^(-alpha - 1),            'r--')
legend('quick','full', 'asymptotic formula')
title('Comparison of "quick" option to full algorithm for large x' )
