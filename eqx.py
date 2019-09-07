var
pi
e_pi 
y 
c 
tune_shock_pi
tune_shock_e1q_pi
varexo 
shock_pi
shock_y
shock_e1q_pi
varobs
obs_y
obs_shock_pi
obs_tune_shock_e1q_pi
parameters
alpha
beta
theta
mpc
paramval
alpha=0.2
beta=0.4
theta=0.9
mpc=0.7

model
pi = alpha*e_pi + (1-alpha)*pi(-5) + beta*y + shock_pi + tune_shock_e1q_pi(-1);
e_pi = pi(+1);
y = theta*y(-1) + shock_y;
c = mpc*y;
tune_shock_pi = shock_pi;
tune_shock_e1q_pi = shock_e1q_pi;

model_obs
obs_y = y;
obs_shock_pi = tune_shock_pi;
obs_tune_shock_e1q_pi = tune_shock_e1q_pi;