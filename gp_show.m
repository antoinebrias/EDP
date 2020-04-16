% GP regression vizualisation

function []=gp_show(mdlstruct,indX,indY)
% name = mdlstruct.name;

if mdlstruct.n_dim>1
    
    if nargin == 1
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
    
    figure;
    
    
    
    name = mdlstruct.name;
    n_dim = mdlstruct.n_dim;
    data = mdlstruct.data;
    
    
    if isfield(mdlstruct,'gp_model')
        [mu_gp,v_gp] = mdlstruct.gp_model(M,M*0,1);
        
    end
    
    hold on;
    
    if  isfield(mdlstruct,'model')
        [mu,v] = mdlstruct.model(M,M*0,1);
    end
    
    [pplot,nplot]=numSubplots(n_dim);
    for i=1:n_dim
        subplot(pplot(1),pplot(2),i);
        
        if isfield(mdlstruct,'gp_model')
            surf(x1d1,x1d2,reshape(mu(:,i),nxgrid,nxgrid))
        end
        
        hold on;
        
        if  isfield(mdlstruct,'model')
            surf(x1d1,x1d2,reshape(mu_gp(:,i),nxgrid,nxgrid))
        end
        
        
        
        
        
        switch mdlstruct.control_type
            case 'rate'
                plot3(mdlstruct.data(1:end-1,indX).*(1-mdlstruct.control_data(1:end-1,indX)),mdlstruct.data(1:end-1,indY).*(1-mdlstruct.control_data(1:end-1,indY)),mdlstruct.data(2:end,i),'rx')
                
            case 'single'
                % to do
            case 'global'
                % to do
                
        end
        
        %         legend(name{i},leg_control);
        
        grid on;  box on
    end
    
else
    ind_available_var= mdlstruct.ind_available_var;
    
    if mdlstruct.n_dim==1
        if mdlstruct.n_lags(ind_available_var(1))>1
            
            
            
        else
            nxgrid=100;
            x=mdlstruct.data;
            x1d = [linspace(0*min(x(:,ind_available_var(1))),max(x(:,ind_available_var(1))),nxgrid)]';
            
            [x1d1] = x1d;
            
            
%             xgrid(:,1)=x1d1(:);
            xgrid=repmat(mean(x),nxgrid,1);
    xgrid(:,ind_available_var(1))=x1d1(:);

    
            
            M = xgrid;
           
            
            figure;
            
            
            
            name = mdlstruct.name;
            n_dim = mdlstruct.n_dim;
            data = mdlstruct.data;
            
            if isfield(mdlstruct,'gp_model')
                [mu_gp,v_gp] = mdlstruct.gp_model(M,M*0,1);
                
            end
            
            hold on;
            
            if  ~isfield(mdlstruct,'var_availability') & isfield(mdlstruct,'model')
                [mu,v] = mdlstruct.model(M,M*0,1);
            end
            
            [pplot,nplot]=numSubplots(n_dim);
            for i=1:n_dim
                subplot(pplot(1),pplot(2),i);
                
                if isfield(mdlstruct,'gp_model')
                    plot(x1d1,mu_gp(:,i),'g','Linewidth',2)
                end
                
                hold on;
                
                if  ~isfield(mdlstruct,'var_availability') & isfield(mdlstruct,'model')
                    plot(x1d1,mu(:,i),'b','Linewidth',2)
                end
                
                
                
                
                
                switch mdlstruct.control_type
                    case 'rate'
                        plot(mdlstruct.data(1:end-1,1).*(1-mdlstruct.control_data(1:end-1,1)),mdlstruct.data(2:end,i),'rx','Linewidth',2)
                        
                    case 'single'
                        % to do
                    case 'global'
                        % to do
                        
                end
                
                %         legend(name{i},leg_control);
                
                grid on;  box on
            end
            
        end
    end
    
    
    
end