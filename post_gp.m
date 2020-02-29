%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gaussian Process posterior mean and variance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mu,v] = post_gp(x,u,is_det,gp_data)

s= x.*(1-u); % need to be changed for different types of control


snorm = (s-gp.mIn)./(gp.sdIn);

for i=1:gp.n_dim
    [~,~,out]=gp.gp_handle{i}(snorm);
    muTmp(:,i)=out.pred*gp_data.sdOut(i)+gp.mOut(i);
    if gp.isLog==1
        muMed(:,i)  = exp( muTmp(:,i)).*s(:,i); %log normal median
        
        
        mu(:,i) = muMed(:,i).*exp(0.5*out.var).*s(:,i); %log normal mean
    else
        mu(:,i) =  muTmp(:,i);
    end
    v(:,i) = out.var*gp.sdOut(i)^2;
    
    
    
end
%     end

% add the uncertainties
if ~is_det
    if gp.isLog==1
        mu = s.*exp(normrnd(muTmp,sqrt(v)));
    else
        mu = normrnd(muTmp,sqrt(v));
    end
end

mu(mu<0)=0;
end







