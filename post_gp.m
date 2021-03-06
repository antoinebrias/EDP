function [mu,v] = post_gp(x,u,is_det,gp)
% POST_GP  Gaussian Process posterior mean and variance.
%   [MU,V] = POST_GP(X,U,IS_DET,GP) return the posterior mean MU and
%   variance V, given the data X, the control U, and the Gaussian Process
%   structure GP.
%   If is_det is 0, generates a sample MU, according the posterior
%   distribution.
%   
%   !!! Need to handle more control rules.



s= x.*(1-u); % !!! need to be changed for different types of control


snorm = (s-gp.mIn)./(gp.sdIn);

for i=1:length(gp.ind_available_var) % !!!! 16/04/21 change gp.n_dim
    [~,~,out]=gp.gp_handle{i}(snorm);
    muTmp(:,i)=out.pred*gp.sdOut(i)+gp.mOut(i);
    if gp.is_log==1
        muMed(:,i)  = exp( muTmp(:,i)).*s(:,i); %log normal median
        
        
        mu(:,i) = muMed(:,i).*exp(0.5*out.var); %log normal mean
    else
        mu(:,i) =  muTmp(:,i);
    end
    v(:,i) = out.var*gp.sdOut(i)^2;
    
    
    
end
%     end

% add the uncertainties
if ~is_det
    if gp.is_log==1
        mu = s.*exp(normrnd(muTmp,sqrt(v)));
    else
        mu = normrnd(muTmp,sqrt(v));
    end
end

% if is_hist
%    % return [X_{t+1} X_{t} X_{t-1} ...] (GP regression)
% else
%    % return only X_{t+1} (value fct GP)
%     
% end


end







