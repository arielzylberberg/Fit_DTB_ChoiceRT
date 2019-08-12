addpath('functions')

load test_data

%% one run
kappa = 15;
ndt_m = 0.2;
ndt_s = 0.01;
B0    = 1;
a     = 0;
d     = 0;
coh0  = 0;
y0    = 0;
plot_flag =1;
pars = [];
theta = [kappa,ndt_m,ndt_s,B0,a,d,coh0,y0];
[err,P] = wrapper_dtb_parametricbound_rt(theta,rt,coh,choice,c,pars,plot_flag);

%% fitting


% kappa, ndt_mu, ndt_sigma, B0, a, d, coh0, y0
tl = [5,  0.1, .01 ,0.5  , -1, -3,0,0];
th = [40, 0.7, .15 ,4    , 4 ,4,0,0];
tg = [15, 0.2, .02 ,1    , 0.1 ,1,0,0];

plot_flag = true;
pars = [];

fn_fit = @(theta) (wrapper_dtb_parametricbound_rt(theta,rt,coh,choice,c,pars,plot_flag));

options = optimset('Display','final','TolFun',.1,'FunValCheck','on');


ptl = tl;
pth = th;
[theta, fval, exitflag, output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);




