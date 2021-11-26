function [mu,v] = post_nn(x,u,nn)
% POST_NN  Nearest neaighboor normalized prediction.
%   [MU,V] = POST_NN(X,U,GP) return the posterior mean MU and
%   variance V, given the data X, and the Nearest Neighbor
%   structure GP.
% work for 1 species only
s= x.*(1-u); % !!! need to be changed for different types of control


snorm = (s-nn.mIn)./(nn.sdIn);

for i=1:length(nn.ind_available_var) % !!!! 16/04/21 change gp.n_dim
    if  size(snorm,2)> size(nn.nn_handle{i}.Y,2)
         [out.pred]=predict(nn.nn_handle{i},snorm(:,1));
             out.var = out.pred.*0;
    else
    [out.pred]=predict(nn.nn_handle{i},snorm);
    out.var = out.pred.*0;
    end
    muTmp(:,i)=out.pred*nn.sdOut(i)+nn.mOut(i);
    
    if nn.is_log==1
        muMed(:,i)  = exp( muTmp(:,i)).*s(:,i); %log normal median
        
        
        mu(:,i) = muMed(:,i).*exp(0.5*out.var); %log normal mean
    else
        mu(:,i) =  muTmp(:,i);
    end
    v(:,i) = out.var*nn.sdOut(i)^2;
    
    
    
end


end







