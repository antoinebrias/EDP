%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GP regression vizualisation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function []=td_show(mdlstruct,optstruct,dpstruct,indX,indY)
% name = mdlstruct.name;

if mdlstruct.n_dim>1
    
    if nargin == 3
        indX = 1;
        indY = 2;
    end
    
    nxgrid=10;
    x=mdlstruct.data;
    x1d = [linspace(0*min(x(:,indX)),max(x(:,indX)),nxgrid);linspace(0*min(x(:,indY)),max(x(:,indY)),nxgrid)]';
    
    [x1d1,x1d2] =  meshgrid(x1d(:,1),x1d(:,2));
    
    
    xgrid=repmat(mean(x),nxgrid^2,1);
    xgrid(:,indX)=x1d1(:);
    xgrid(:,indY)=x1d2(:);
    
    M = xgrid;
    
    
    
    [opt_control,unweighted_opt_value,weighted_opt_value]=td_policy(M,optstruct,dpstruct);
    
    
    figure;
    
    
    
    name = mdlstruct.name;
    n_dim = mdlstruct.n_dim;
    data = mdlstruct.data;
    
    

    [pplot,nplot]=numSubplots(n_dim);
    for i=1:n_dim
            subplot(pplot(1),pplot(2),i);
        
            contourf(x1d1,x1d2,reshape(opt_control(:,i),nxgrid,nxgrid))

        
        
        grid on;  box on
    end
    
else
    
    % to do
end



% if isTD
%     
%     [upred,~,~]=tdPolicy(xgrid,gamma,Mdl,MdlV,ugrid,paramReward);
% else
%     oneSpMdl=MdlV;
%     [upred,vres,vUnweightedTmp]=oneSpPolicy(xgrid,oneSpMdl,paramReward);
%     
%     
% end
% 
% 
% for i=1%:size(ugrid,2)
%     ctrlMap{i} = reshape(upred(:,i),nxgrid,nxgrid);
%     % subplot(1,size(ugrid,2),i)
%     contourf(x1d1,x1d2,ctrlMap{i})
% end

end