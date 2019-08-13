function [out, p, rt_prctiles] = rt_quant_plot(coh,rt,correct, do_plot)
% function rt_quant_plot(coh,rt,choice,req_choice)
% Does the rt-quantile plot ala Ratcliff

if nargin<4 || isempty(do_plot)
    do_plot = 1;
end

rt_prctiles = [10 30 50 70 90];
nprc   = length(rt_prctiles);
colors = [0.8,0,0; 0,0.8,0]; % colors for error and correct responses

min_trials_threshold = 5; % min number of trials per condition to include the data in the plot

%%
% correct = choice==req_choice;
[~,~,idx] = unique(abs(coh));

for j=1:2 % incorrect, correct
    if j==1
        depvar = 1 - correct;
    else
        depvar = correct;
    end
    [~,c] = curva_media(depvar,abs(coh),[],0); % average prop of correct or incorrect per coherence level
    uc = unique(c(c>0)); % unique value of non-zero correct/incorrect
    val = index_to_val(idx,c); % replaces the coh by the av correct or incorrect
    
    % fills with nans the values that we don't care about
    if j==1
        to_nan = correct>0;
    else
        to_nan = correct==0;
    end
    val(to_nan) = nan;
    
    out(j).rt_prc = nan(nprc,length(uc));
    out(j).uc = uc;
    for i=1:length(uc)
        I = val==uc(i);
        if sum(I)>min_trials_threshold
            out(j).rt_prc(:,i) = prctile(rt(I),rt_prctiles);
        end
    end
    
end


if (do_plot)
    p = publish_plot(1,1);
    for j=1:2
        plot(out(j).uc,out(j).rt_prc','linestyle','-','color','k','color',colors(j,:));
        hold all
        plot(out(j).uc,out(j).rt_prc','linestyle','none','marker','x','color',colors(j,:));
        
        hold all
    end
    
    xlabel('Response proportion')
    ylabel('Response time (s)')
    
    p.format();
else
    p = [];
end



end