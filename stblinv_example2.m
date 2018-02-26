clear; close all; drawnow

u = [.01:.01:.99];
alpha = 1.5;
beta = .5;
gam = 1;
delta = 0;
w = stbl.cdf( stbl.inv( u,alpha,beta,gam,delta),alpha,beta,gam,delta) % This should return u.
check = all( abs(w-u) < 1e-6 )
