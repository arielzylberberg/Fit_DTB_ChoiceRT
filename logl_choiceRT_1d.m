function [nlogl,pPred] = logl_choiceRT_1d(P,choice,rt,coh,ndt_m,ndt_s)
% Computes the negative log-likelihood of the parameters given the choice
% and reaction time data. 
% Comments to: ariel.zylberberg@gmail.com

t = P.t;

%convolve for non-decision times
%sanity check
if t(1)~=0
    error('for conv to work, t(1) has to be zero');
end
nt = length(t);
dt = t(2)-t(1);
ntr = length(P.drift);

ndt = normpdf(t,ndt_m,ndt_s)*dt;
upRT = conv2(1,ndt(:),P.up.pdf_t);
loRT = conv2(1,ndt(:),P.lo.pdf_t);
upRT = upRT(:,1:nt);
loRT = loRT(:,1:nt);

rt_step = ceil(rt/dt);
ucoh = unique(coh);
ncoh = length(ucoh);
p_up = nan(ntr,1);
p_lo = nan(ntr,1);
for i=1:ncoh
    inds = coh == ucoh(i);
    %if RT too long, clamp to last
    J = min(rt_step(inds),nt);
    p_up(inds) = upRT(i,J);
    p_lo(inds) = loRT(i,J);
end

%clip
p_up(p_up<eps) = eps;
p_lo(p_lo<eps) = eps;
pPred = p_up.*(choice==1) + p_lo.*(choice==0);
nlogl = -nansum(log(pPred));