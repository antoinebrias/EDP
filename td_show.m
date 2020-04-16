function []=td_show(mdlstruct,optstruct,dpstruct,indX,indY)
% TD_SHOW  displays the control map obtained by the EDP policy
%   []=TD_SHOW(MDLSTRUCT,OPTSTRUCT,DPSTRUCT,INDX,INDY)
%    returns the slice of control map according the axis INDX and INDY
%
% !!! Need to handle more cases

% name = mdlstruct.name;
n_lags = mdlstruct.n_lags;
n_dim= mdlstruct.n_dim;
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
    
    
    
    [opt_control,unweighted_opt_value,weighted_opt_value]=td_policy(M,optstruct,dpstruct,mdlstruct);
    
    
    figure;
    
    
    
    name = mdlstruct.name;
    n_dim = mdlstruct.n_dim;
    data = mdlstruct.data;
    
    
    
    [pplot,nplot]=numSubplots(n_dim);
    for i=1:n_dim
        subplot(pplot(1),pplot(2),i);
        
        contourf(x1d1,x1d2,reshape(opt_control(:,i),nxgrid,nxgrid))
          xlabel('lag 0')
            ylabel('lag 1')
        
        
        grid on;  box on
    end
    
else
    
    % to do
    
    if  mdlstruct.n_lags == 1
        
    else
        % by default the two first lags are displayed
        indX = 1;
        indY = 2;
        nxgrid=10;
        x=mdlstruct.data;
        x_tmp=[];
        for i=1:1%n_lags(1) %!!!
            x_tmp = [x_tmp linspace(0*min(x(:,indX)),max(x(:,indX)),nxgrid)'];
            
        end
        x1d = x_tmp;
        [x1d1,x1d2] =  meshgrid(x1d,x1d);
        
        x_grid_tmp = [];
        for i=1:n_dim
            x_grid_tmp = [x_grid_tmp repmat(mean(x(:,i)),nxgrid^2,n_lags(i))];
            
        end
        
        x_grid=x_grid_tmp;
        x_grid(:,indX)=x1d1(:);
        x_grid(:,indY)=x1d2(:);
        
        
        
        [opt_control,unweighted_opt_value,weighted_opt_value]=td_policy(x_grid,optstruct,dpstruct,mdlstruct);
        
        
        figure;
        
        
        
        name = mdlstruct.name;
        n_dim = mdlstruct.n_dim;
        data = mdlstruct.data;
        
        
        
        [pplot,nplot]=numSubplots(n_dim);
        for i=1:n_dim
            subplot(pplot(1),pplot(2),i);
            
            contourf(x1d1,x1d2,reshape(opt_control(:,i),nxgrid,nxgrid))
            
            xlabel('lag 0')
            ylabel('lag 1')
            
            grid on;  box on
        end
        
    end
    
    
    
    
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