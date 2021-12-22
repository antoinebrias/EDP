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

 %%%%%% CHANGE 18/11/21
 if ~gp.is_value_function
s = log(s+eps);
 end
 
snorm = (s-gp.mIn)./(gp.sdIn);

for i=1:length(gp.ind_available_var) % !!!! 16/04/21 change gp.n_dim
    [~,~,out]=gp.gp_handle{i}(snorm);
    muTmp(:,i)=out.pred*gp.sdOut(i)+gp.mOut(i);
        v(:,i) = out.var(:,i)*gp.sdOut(i)^2;
        
    if gp.is_log==1
        muMed(:,i)  = exp( muTmp(:,i)).*s(:,i); %log normal median
        
        
        mu(:,i) = muMed(:,i).*exp(0.5*out.var(:,i)); %log normal mean
        
        % add the uncertainties
        if ~is_det
                mu(:,i) = s(:,i).*exp(normrnd(muTmp(:,i),sqrt(v(:,i))));
        end
        
    else
        mu(:,i) =  muTmp(:,i);
        
        % add the uncertainties
        if ~is_det
                mu(:,i) = normrnd(muTmp(:,i),sqrt(v(:,i)));
    
        end
        
    end

    

    
end
%     end



     %%%%%% CHANGE 18/11/21
     if ~gp.is_value_function && ~gp.is_log
mu = exp(mu+eps);
     end

% if is_hist
%    % return [X_{t+1} X_{t} X_{t-1} ...] (GP regression)
% else
%    % return only X_{t+1} (value fct GP)
%
% end


end







