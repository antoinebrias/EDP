function []=data_show(mdlstruct)
% DATA_SHOW  Display the available data.
%   []=DATA_SHOW(MDLSTRUCT) displays data in MDLSTRUCT.
%
%   See also GP_SHOW, TD_SHOW.

name = mdlstruct.name;
n_dim = mdlstruct.n_dim;
data = mdlstruct.data;

[pplot,nplot]=numSubplots(n_dim);
for i=1:n_dim
    subplot(pplot(1),pplot(2),i);
    hold on
    plot(data(:,i),'b','Linewidth',2);
    
    switch mdlstruct.control_type
        case 'rate'
            plot(data(:,i).*(mdlstruct.control_data(:,i)),'c','Linewidth',2);
            leg_control = [name{i} ' catch'];
        case 'single'
            % to do
        case 'global'
            % to do
            
    end
    
    legend(name{i},leg_control);
    
    grid on;  box on
end

end