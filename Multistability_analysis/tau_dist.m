clc
clear all

mouselist = ["M19" "M20" "M29" "M34" "M35" "M39"];

for iMouse = 1:length(mouselist)
    
    animal = mouselist(iMouse); % mouse number
    pathall = strcat('G:\My Drive\Fenton_Lab\CellPress_CircleSquare\MNS_Analysis\Additional_codes\Raw_data','\','Dist_',mouselist(iMouse),'.mat');
    load(pathall);
    
    
    for isession = 1:size(Dist,2)
        
        distance{iMouse,isession} = Dist{2,isession};
        tau{iMouse,isession} = Dist{1,isession}(:,1);
        pvalues{iMouse,isession} = Dist{1,isession}(:,2);
        
    end
    
end

%% figure
tiledlayout(2,3);

for iMouse = 1:length(mouselist)
    
    cumulative_distance = [reshape(distance{iMouse,1},[],1)];
    cumulative_tau = [reshape(tau{iMouse,1},[],1)];
    cumulative_pvalue= [reshape(pvalues{iMouse,1},[],1)];
    
    
    for isession = 2:size(Dist,2)

        cumulative_distance = [cumulative_distance; reshape(distance{iMouse,isession},[],1)];
        cumulative_tau= [cumulative_tau; reshape(tau{iMouse,isession},[],1)];
        cumulative_pvalue= [cumulative_pvalue; reshape(pvalues{iMouse,isession},[],1)];
       
    end
    
        remove_zeros = find(cumulative_distance>0);
        cumulative_distance = cumulative_distance(remove_zeros);
        cumulative_tau = cumulative_tau(remove_zeros);
        cumulative_pvalue = cumulative_pvalue(remove_zeros);
    
%         figure()
%         scatter(cumulative_distance,cumulative_tau,20,'.','MarkerFaceColor','black','MarkerEdgeColor','black')
%         grid on
%         xlabel('distance')
%         ylabel('tau')
        
        nexttile
        [fitresult{iMouse}, gof{iMouse}] = power_fit(cumulative_distance, cumulative_tau, cumulative_pvalue)
        title(mouselist(iMouse),'FontSize',20);
        
        xlabel( 'Interneuron distance', 'Interpreter', 'latex' );
        ylabel( '$\tau$', 'Interpreter', 'latex' );        
        xticks([0 150 300 450 600])
        ylim([min(cumulative_tau) 1])
        
        a_parameter(iMouse) = fitresult{1,iMouse}{1,1}.a;
        b_parameter(iMouse) = fitresult{1,iMouse}{1,1}.b;
        r2_parameter(iMouse) = gof{1,iMouse}.rsquare;
        rmse_parameter(iMouse) = gof{1,iMouse}.rmse;
        
end

mean_a = mean(a_parameter);
mean_b = mean(b_parameter);
std_a = std(a_parameter);
std_b = std(b_parameter);
range_a = [min(a_parameter) max(a_parameter)];
range_b = [min(b_parameter) max(b_parameter)];

%% figure
figure()
scatter(distance(1:200),tau(1:200))
grid on
xlabel('distance')
ylabel('tau')


%%
figure()
histfit(tau)

%%

[prob,edges] = histcounts(tau,'Normalization','probability');


log_distance = log(distance);

for iI = 1:length(tau)
    
    log_tau(iI) = log(abs(tau(iI)))*sign(tau(iI));
    
end


%% figure
figure()
scatter(log_distance,log_tau)
grid on
xlabel('log-distance','FontSize',20)
ylabel('log-tau','FontSize',20)

%% figure
figure()
scatter(log_distance(1:90),log_tau(1:90))
grid on
xlabel('log-distance','FontSize',20)
ylabel('log-tau','FontSize',20)


%%

rangex = min(tau)
rangex2 = max(tau)

x = linspace(rangex,rangex2,100)
p = cdf(pd,x);


%%

indx = find(distance>40);

newdistance = distance(indx);
newtau = tau(indx);


%% single recording

%% remove negative

remove_zeros = find(distance{6,26}>0);
distance_single = distance{6,26}(remove_zeros);
tau_single = tau{6,26}(remove_zeros);
pvalues_single = pvalues{6,26}(remove_zeros);

figure()
[fitresult_single, gof_single] = power_fit(distance_single, tau_single, pvalues_single)
xlabel( 'Interneuron distance', 'Interpreter', 'latex' );
ylabel( '$\tau$', 'Interpreter', 'latex' );        
xticks([0 150 300 450 600])
ylim([tau_single 1])
texttittle = strcat(mouselist(6),'Day5','CYL1');
title(texttittle,'FontSize',20);
