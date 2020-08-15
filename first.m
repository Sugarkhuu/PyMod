% add IRIS to the path
% addpath /home/sugarkhuu/Documents/IRIS-Toolbox-Release-20200119
%% edited from github
addpath /home/sugarkhuu/Documents/IRIS-Toolbox-Release-20180319

% configure IRIS
irisstartup


p = struct();
p.alpha=0.2;
p.beta=0.4;
p.theta=0.9;
p.mpc=0.7;

m = model('first.mod','linear',true,'assign',p);


m = solve(m);
m = sstate(m);



d=struct();
d.obs_x = tseries();
% d.obs_x = tseries(qq(2018,1):qq(2019,4),[1.5 1.2 0.9 0.8 1.5 2.1 1.3 1.2]);
d.obs_x(qq(2018,1):qq(2019,4)) = [1.5 1.2 0.9 0.8 1.5 2.1 1.3 1.2];
d.e = tseries();
d.e(qq(2018,1)) = 5;
df = d;

[~, f, v, ~, pe, co] = filter(m, d, qq(2018,1):qq(2019,4)+10);


s = zerodb(m, qq(2018,1):qq(2019,4));
% s.e(qq(2018,1)) = 5;
s.e(qq(2018,4)) = 5;
% res = simulate(m,s,qq(2018,1):qq(2019,4),'deviation',true,'anticipate',false);
res = simulate(m,s,qq(2018,1):qq(2019,4),'deviation',true,'anticipate',true);
