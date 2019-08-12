function [err,P] = wrapper_dtb_parametricbound_rt(theta,rt,coh,choice,c,pars,plot_flag)
% function [err,P] = wrapper_dtb_parametricbound_rt(theta,rt,coh,choice,c,pars,plot_flag)

% INPUTS
% ~~~~~~~~~~~
% theta. Model parameters: [kappa,ndt_m, ndt_s, B0, a, d, coh0, y0]
% 'kappa': signal-to-noise ratio.
% 'ndt_m' and 'ndt_s': mean and standard deviation of the non-decision time (assumed Gaussian). 
% 'B0', 'a', 'd': parameters that specify the shape of the bounds.
% 'coh0' and'y0': bias terms, the first as a change in the drift rate, and the second one as an offset in the starting point of the accumulation. 
% 'pars': struct with optional flags and parameters. If 'pars.notabs_flag = True', the output includes the distribution of non-absorved particles.
% Most useful for variable-duration experiments. 
% 'pars.USfunc' specifies the shape of the bound (see the 'expand_bounds' function. 
% 'pars.t' specifies the times at which to solve the Fokker-Planck
% equation (otherwise, it uses a default).
% 'plot_flag': boolean. If true, a plot the data & fit is generated and
% updated during the optimization.

% OUPUTS
% ~~~~~~~~~~~
% 'err': negative log likelihood of the params given the data
% 'P': struct containing the model's distribution of decision times per difficulty
% level, time course of the bounds, expected accuracy, ...

% Written by ariel zylberberg (ariel.zylberberg@gmail.com)


%%
kappa  = theta(1);
ndt_m  = theta(2);
ndt_s  = theta(3);
B0     = theta(4);
a      = theta(5);
d      = theta(6);
coh0   = theta(7);
y0a    = theta(8);

%% save the unabsorbed
if ~isempty(pars) && isfield(pars,'notabs_flag')
    notabs_flag = pars.notabs_flag;
else
    notabs_flag = false;
end

%% set the time vector

if ~isempty(pars) && isfield(pars,'t')
    t = pars.t;
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

%% discretize space and initialize the decision variable

y  = linspace(min(Blo)-0.3,max(Bup)+0.3,750)';

y0a = clip(y0a,Blo(1),Bup(1));

y0 = zeros(size(y));
y0(findclose(y,y0a)) = 1;
y0 = y0/sum(y0);


%% numerical solution using Chang-Cooper's method
drift = kappa * unique(coh + coh0);
P = dtb_fp_cc_vec(drift,t,Bup,Blo,y,y0,notabs_flag);

%% likelihood

err = logl_choiceRT_1d(P,choice,rt,coh,ndt_m,ndt_s);


%% print
fprintf('err=%.3f kappa=%.2f ndt_mu=%.2f ndt_s=%.2f B0=%.2f a=%.2f d=%.2f coh0=%.2f y0=%.2f \n',...
    err,kappa,ndt_m,ndt_s,B0,a,d,coh0,y0a);


%%

if plot_flag
    m = prctile(rt,99.5);
    
    figure(1);clf
    set(gcf,'Position',[293  388  828  560])
    
    subplot(2,2,1);
    plot(t,P.Bup,'k');
    hold all
    plot(t,P.Blo,'k');
    title('Bounds (t)')
    if ~isnan(m)
        xlim([0,m])
    end
    ylabel('Accumulated evidence')
    xlabel('Time (s)')
    
    subplot(2,2,2);
    [tt,xx,ss] = curva_media(choice,coh,[],0);
    terrorbar(tt,xx,ss,'color','k','LineStyle','none','Marker','.');
    hold all
    ucoh = unique(coh);
    plot(ucoh,P.up.p,'k-');
    xlim([min(ucoh),max(ucoh)])
    xlabel('Motion coherence');
    ylabel('P rightward');
    
    subplot(2,2,3);
    [out,~,rt_prctiles] = rt_quant_plot(coh(choice==1),rt(choice==1),c(choice==1), 0);
    hold all
    colors = [0.8,0,0; 0,0.8,0]; % colors for error and correct responses
    for j=1:2
        plot(out(j).uc,out(j).rt_prc','linestyle','-','color','k','color',colors(j,:));
        hold all
        plot(out(j).uc,out(j).rt_prc','linestyle','none','marker','x','color',colors(j,:));
    end
    
    
    rt_prc_up = nan(length(ucoh),length(rt_prctiles));
    for i=1:length(ucoh)
        if sum(P.up.cdf_t(i,:))>0
            pup_cum_norm = P.up.cdf_t(i,:)/P.up.cdf_t(i,end);
            idx = findclose(pup_cum_norm,rt_prctiles/100);
            rt_prc_up(i,:) = P.t(idx) + ndt_m; %not exact, should truncate
        end
    end
    plot(P.up.p,rt_prc_up,'k')
    
%     P.up.rt_prc_up = rt_prc_up; 
%     P.up.rt_prctiles = rt_prctiles;
    
    xlabel('Response proportion')
    ylabel('Response time (s)')
    
    title('Rightward choices only')
    
    
    subplot(2,2,4);
    
    ind = P.drift>=0;
    rt_model_c(ind) = P.up.mean_t(ind) + ndt_m;
    rt_model_nc(ind) = P.lo.mean_t(ind) + ndt_m;
    
    ind = P.drift<0;
    rt_model_c(ind) = P.lo.mean_t(ind) + ndt_m;
    rt_model_nc(ind) = P.up.mean_t(ind) + ndt_m;
    
    % don't plot errors if there probability is too low
    ind = P.up.p>1e-5 &  P.up.p<(1-1e-5);
    rt_model_nc(~ind) = nan;
    
    [tt,xx,ss] = curva_media(rt,coh,c~=0,0);
    terrorbar(tt,xx,ss,'color','k','LineStyle','none','Marker','.');
    hold all
    plot(ucoh,rt_model_c,'k-');
    
    [tt,xx,ss] = curva_media(rt,coh,c==0,0);
    terrorbar(tt,xx,ss,'color','r','LineStyle','none','Marker','.');
    hold all
    plot(ucoh,rt_model_nc,'r-');
    
    xlim([min(ucoh),max(ucoh)])
    xlabel('Motion coherence');
    ylabel('Response time');
    
    format_figure(gcf);
    
    drawnow
    
end


