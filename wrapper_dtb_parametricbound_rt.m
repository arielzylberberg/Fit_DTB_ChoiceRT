function [err,P] = wrapper_dtb_parametricbound_rt(theta,rt,coh,choice,c,pars,plot_flag)
% function [err,P] = wrapper_dtb_parametricbound_rt(theta,rt,coh,choice,c,pars,plot_flag)
% written by ariel zylberberg (ariel.zylberberg@gmail.com)

%%
kappa  = theta(1);
ndt_m  = theta(2);
ndt_s  = theta(3);
B0     = theta(4);
a      = theta(5);
d      = theta(6);
coh0   = theta(7);
y0a    = theta(8);

%%
if ~isempty(pars) && isfield(pars,'notabs_flag')
    notabs_flag = pars.notabs_flag;
else
    notabs_flag = false;
end

%%

if ~isempty(pars) && isfield(pars,'t')
    t = pars.t;
    %dt = t(2)-t(1);
else
    dt = 0.0005;
    t  = 0:dt:10;
end

%% bounds
if ~isempty(pars) && isfield(pars,'USfunc')
    USfunc = pars.USfunc;
else
    USfunc = 'Exponential';
end
[Bup,Blo] = expand_bounds(t,B0,a,d,USfunc);

%%

% y  = linspace(min(Blo)-0.3,max(Bup)+0.3,512)';%esto puede ser problematico
y  = symmetric_scale(max(Bup)*1.2,0.004);

y0 = zeros(size(y));
y0(findclose(y,0)) = 1;
y0 = y0/sum(y0);


%%
% prior = Rtable(coh)/sum(Rtable(coh));

%%
drift = kappa * unique(coh + coh0);

if ~isempty(pars) && isfield(pars,'methodforward_flag')
    methodforward_flag = pars.methodforward_flag;
else
    methodforward_flag = 1;%chang cooper
end
switch methodforward_flag
    case 1 % chang cooper
        P = dtb_fp_cc_vec(drift,t,Bup,Blo,y,y0,notabs_flag);
    case 2 % fft
        P = spectral_dtbAA(drift(:)',t,Bup,Blo,y,y0,notabs_flag);
end


%% likelihood
err = logl_choiceRT_1d(P,choice,rt,coh,ndt_m,ndt_s);

%% print
fprintf('err=%.3f kappa=%.2f ndt_mu=%.2f ndt_s=%.2f B0=%.2f a=%.2f d=%.2f coh0=%.2f y0=%.2f \n',...
    err,kappa,ndt_m,ndt_s,B0,a,d,coh0,y0a);

%%
if plot_flag
    
    figure(11);clf
    set(gcf,'Position',[241   249   582   202])
    
    subplot(1,2,1);
    curva_media(choice,coh,[],1);
    hold all
    ucoh = unique(coh);
    plot(ucoh,P.up.p,'.-');
    xlabel('Motion coherence (%)')
    ylabel('P rightward choice')
    xlim([-0.6,0.6])
    hl = legend('data','model');
    set(hl,'location','best');
    
    subplot(1,2,2);
    rt_model = (P.up.mean_t.*P.up.p+P.lo.mean_t.*P.lo.p)./(P.up.p+P.lo.p) + ndt_m; %Just for plotting; is not exact because it
    % doesn't take into account the curtailing of the non-decision time
    % distribution
    curva_media(rt,coh,[],1);
    hold all
    plot(ucoh,rt_model,'.-');
    xlabel('Motion coherence (%)')
    ylabel('Response time (s)')
    xlim([-0.6,0.6])
    format_figure(gcf,'FontSize',14);
    
    drawnow
    
end