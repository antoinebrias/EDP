% clear all
% close all

% load('save_gp_11072021.mat')
% load('compGP-290821.mat')

res = resTmp;
% parset harvestfun model vset  Tset iter


%% init
res_chaos  = [];res_limit_cycle  = [];res_stable  = [];res_noise_level_0  = [];res_noise_level_1  = [];res_noise_level_2  = [];
res_training_interval_30  = [];res_training_interval_50  = [];res_training_interval_100  = [];res_two_locations  = [];res_mat_effect  = [];res_time_varying_r  = [];res_competition  = [];
iter_res = 0;
for parset=1:3
for harvestfun=3%1:3    
for model=[2 3 4 10]
for vset=1:3
for Tset=1:3
for iter=1:100
    iter_res = iter_res +1;
    %% dynamic regime
    if parset == 1 & harvestfun == 3
    res_chaos = [res_chaos; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    if parset == 2 & harvestfun == 3
    res_limit_cycle = [res_limit_cycle; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    if parset == 3 & harvestfun == 3
    res_stable = [res_stable; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    
    %% noise level
    if vset == 1 & harvestfun == 3
    res_noise_level_0 = [res_noise_level_0; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    if vset == 2 & harvestfun == 3
    res_noise_level_1 = [res_noise_level_1; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    if vset == 3 & harvestfun == 3
    res_noise_level_2 = [res_noise_level_2; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    
       %% training interval
    if Tset == 1 & harvestfun == 3
    res_training_interval_30 = [res_training_interval_30; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    if Tset == 2 & harvestfun == 3
    res_training_interval_50 = [res_training_interval_50; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    if Tset == 3 & harvestfun == 3
    res_training_interval_100 = [res_training_interval_100; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    
    
        %% scenario
    if model == 2 & harvestfun == 3
    res_two_locations = [res_two_locations; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    if model == 3 & harvestfun == 3
    res_mat_effect = [res_mat_effect; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    if model == 4 & harvestfun == 3
    res_time_varying_r= [res_time_varying_r; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
      if model == 10 & harvestfun == 3
    res_competition = [res_competition; resTmp(iter_res,3)./resTmp(iter_res,1) resTmp(iter_res,4)./resTmp(iter_res,1) resTmp(iter_res,3)./resTmp(iter_res,2) resTmp(iter_res,4)./resTmp(iter_res,2)];
    end
    
    
end;end;end;end;end;end

% disp("dynamic regime")
a  = quantile(res_chaos,[0.25, 0.5, 0.75]);
fprintf( 'Dynamic & Chaos & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))

a = quantile(res_limit_cycle,[0.25, 0.5, 0.75]);
fprintf( 'regime & Limit cycle & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))

a=quantile(res_stable,[0.25, 0.5, 0.75]);
fprintf( ' & Stable & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n\\hline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))


% disp("noise level")
a=quantile(res_noise_level_0,[0.25, 0.5, 0.75]);
fprintf( 'Noise & 0.0 & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))

a=quantile(res_noise_level_1,[0.25, 0.5, 0.75]);
fprintf( 'level & 0.1 & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))

a=quantile(res_noise_level_2,[0.25, 0.5, 0.75]);
fprintf( ' & 0.2 & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n\\hline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))


% disp("training interval")
a=quantile(res_training_interval_30,[0.25, 0.5, 0.75]);
fprintf( 'Training & 30 & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))

a=quantile(res_training_interval_50,[0.25, 0.5, 0.75]);
fprintf( 'interval & 50 & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))

a=quantile(res_training_interval_100,[0.25, 0.5, 0.75]);
fprintf( ' & 100 & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n\\hline\n ',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))


% disp("scenario")
a=quantile(res_time_varying_r,[0.25, 0.5, 0.75]);
fprintf( 'Scenario & Time varying r & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))

a=quantile(res_mat_effect,[0.25, 0.5, 0.75]);
fprintf( ' & Mat. effect & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))

a=quantile(res_two_locations,[0.25, 0.5, 0.75]);
fprintf( ' & Two locations & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))

a=quantile(res_competition,[0.25, 0.5, 0.75]);
fprintf( ' & Competition & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)  & %.2f (%.2f , %.2f)  & %.2f (%.2f, %.2f)\\tabularnewline\n\\hline\n',  a(2,1),a(1,1),a(3,1), a(2,3),a(1,3) ,a(3,3),a(2,2),a(1,2),a(3,2), a(2,4),a(1,4) ,a(3,4))
