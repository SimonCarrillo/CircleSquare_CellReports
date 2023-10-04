clc
clear all

set(0,'DefaultFigureWindowStyle','docked')

path_all = '/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/MNS_Analysis/Data/CircleSquare/';

save_figures = 0; % 1 for saving

%% experiment to analyze

mouselist = ["M29" "M34" "M39" "M35" "M20" "M19"];
ratiolist = ["8" "8" "8" "8" "8" "8"]; % "4" "4" "4" "4" "4" "4"
timestep = 60;
ratio_disc = 8;

for iMouse = 1:length(mouselist)
    
    animal = mouselist(iMouse); % mouse number
    
    time_bin = '1000'; % in ms
    expt = 'CircleSquare'; % experiment setup
    loadStr = 'TempCorr'; % 'Zlinear' or 'TempCorr'
    
    save_output = 0; % 1 for saving angles and proportion of overlap
    
    discard_tech= ["kc" "neg" "randPos"];
    
    excel_path = strcat(path_all, 'FilenameMatch.xlsx');
    T = readtable(excel_path,'Sheet',animal);
    
    for idiscard = 1:3
        
        discardPop = discard_tech(idiscard);
        
        ratio = ratiolist(iMouse);
        remove_outliers = 0; %1 for removing outliers
        alpharadius = Inf; % for alphashape
        ISO_output = 1; % 1 for plot
        
        path_ISO = strcat('/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/Response_Jan2022/Analysis_codes/Final_figures/Raw_data/spiking/','ISO_','3_','CellDiscarded_',discardPop,'Ratio',ratio,loadStr,time_bin,expt,'_',animal,'_','.mat');
        load(path_ISO);
        
        %% retrieve daily data
        
        days{iMouse} = linspace(1,size(ISO,1),size(ISO,1));
        
        for iDay = days{iMouse}
            
            ISOmap{iMouse,iDay}{idiscard} = ISO{iDay,2}; % extract dim reduced ISOMAP
            spikes{iMouse,iDay}{idiscard}= ISO{iDay,1}; % extract dim reduced ISOMAP
            sessList{iMouse,iDay}{idiscard}= ISO{iDay,4}; % extract session information
            sessStrList{iMouse,iDay}{idiscard} = ISO{iDay,5}; % extract session information
            tau{iMouse,iDay}{idiscard} = ISO{iDay,8}(:,1);
            pairs{iMouse,iDay}{idiscard} = ISO{iDay,9};
            
            index_nan{iMouse,iDay}{idiscard} = ~isnan(tau{iMouse,iDay}{idiscard});
            tau{iMouse,iDay}{idiscard} = tau{iMouse,iDay}{idiscard}(index_nan{iMouse,iDay}{idiscard});
            pairs{iMouse,iDay}{idiscard} = pairs{iMouse,iDay}{idiscard}(index_nan{iMouse,iDay}{idiscard},:);
            cell_list{iMouse,iDay}{idiscard} = ISO{iDay,11};
            
            if idiscard == 2 || idiscard == 3
                
                cells_removed{iMouse,iDay}{idiscard} = double(ISO{iDay,12}');
            else
                
                cells_removed{iMouse,iDay}{idiscard} = [0];
            end
            
            
            ncells{iMouse,iDay}{idiscard} = size(spikes{iMouse,iDay}{idiscard},2);
            
        end
        
        for iDay = days{iMouse}
            
            ISO_data = ISOmap{iMouse,iDay}{idiscard};
            sess = sessList{iMouse,iDay}{idiscard};
            sess_Str = sessStrList{iMouse,iDay}{idiscard};
            spikes_data = spikes{iMouse,iDay}{idiscard};
            
            %% figures ISO
            
            if ISO_output == 1
                
                index1 = 0;
                index2 = 0;
                index3 = 0;
                
                for Ii = 1:length(ISO_data(:,1))
                    
                    if sess_Str(Ii,1) == 'H' && (sess_Str(Ii,4) == '1' || sess_Str(Ii,4) == '2')
                        
                        index1 = index1 + 1;
                        
                        HMC_store{iMouse,iDay}{idiscard}(index1,:,:,:) = ISO_data(Ii,:);
                        HMC_spikes{iMouse,iDay}{idiscard}(index1,:,:,:) = spikes_data(Ii,:);
                        
                        
                    elseif sess_Str(Ii,1) == 'R' && (sess_Str(Ii,4) == '1' || sess_Str(Ii,4) == '2')
                        
                        index2 = index2 + 1;
                        
                        RCT_store{iMouse,iDay}{idiscard}(index2,:,:,:) = ISO_data(Ii,:);
                        RCT_spikes{iMouse,iDay}{idiscard}(index2,:,:,:) = spikes_data(Ii,:);
                        
                        
                    elseif sess_Str(Ii,1) == 'C' && (sess_Str(Ii,4) == '1' || sess_Str(Ii,4) == '2')
                        
                        index3 = index3 + 1;
                        
                        CYL_store{iMouse,iDay}{idiscard}(index3,:,:,:) = ISO_data(Ii,:);
                        CYL_spikes{iMouse,iDay}{idiscard}(index3,:,:,:) = spikes_data(Ii,:);
                        
                    end
                    
                end
                
            end
            
            % remove outliers or not
            if remove_outliers == 1
                
                HMC = rmoutliers(HMC_store{iMouse,iDay},'quartiles');
                CYL = rmoutliers(CYL_store{iMouse,iDay},'quartiles');
                RCT = rmoutliers(RCT_store{iMouse,iDay},'quartiles');
                
            else
                
                HMC{iMouse,iDay}{idiscard} = HMC_store{iMouse,iDay}{idiscard};
                CYL{iMouse,iDay}{idiscard} = CYL_store{iMouse,iDay}{idiscard};
                RCT{iMouse,iDay}{idiscard} = RCT_store{iMouse,iDay}{idiscard};
                
            end
            
              %% compute most anticorrelated cells per day per mouse
            
            idx_anti{iMouse,iDay}{idiscard} = find(tau{iMouse,iDay}{idiscard}<-0.05);
            cellpairs_anti_tmp{iMouse,iDay}{idiscard} = reshape(pairs{iMouse,iDay}{idiscard}(idx_anti{iMouse,iDay}{idiscard},:),[],1);
            cellpairs_anti{iMouse,iDay}{idiscard} = cellpairs_anti_tmp{iMouse,iDay}{idiscard}(~ismember(cellpairs_anti_tmp{iMouse,iDay}{idiscard},cells_removed{iMouse,iDay}{idiscard})==1);
            
            [cell_count{iMouse,iDay}{idiscard},anti_cell_idx{iMouse,iDay}{idiscard},~] = groupcounts(reshape(cellpairs_anti{iMouse,iDay}{idiscard},[],1));
            [sorted_anti{iMouse,iDay}{idiscard},idx_sorted{iMouse,iDay}{idiscard}] = sort(cell_count{iMouse,iDay}{idiscard},'descend');
            anticofiring_fields{iMouse,iDay}{idiscard} = (anti_cell_idx{iMouse,iDay}{idiscard}(idx_sorted{iMouse,iDay}{idiscard}));
            
            %% remove the cells we actually remove
                                    
            [~,~,real_id_anti{iMouse,iDay}{idiscard}] = intersect(anticofiring_fields{iMouse,iDay}{idiscard},cell_list{iMouse,iDay}{idiscard}','stable');

            anti_power{iMouse,iDay}{idiscard} = cell_count{iMouse,iDay}{idiscard}/ncells{iMouse,iDay}{1};
            
            max_anti_power{iMouse,iDay}{idiscard} = max(anti_power{iMouse,iDay}{idiscard});
            min_anti_power{iMouse,iDay}{idiscard} = min(anti_power{iMouse,iDay}{idiscard});
            
        end
        
    end
    
    
    
    for iDay = days{iMouse}
        
         anti_powerkc{iMouse,iDay}= cell_count{iMouse,iDay}{1}/ncells{iMouse,iDay}{1};
         anti_powerneg{iMouse,iDay}= cell_count{iMouse,iDay}{2}/ncells{iMouse,iDay}{1};
         anti_powerrand{iMouse,iDay}= cell_count{iMouse,iDay}{3}/ncells{iMouse,iDay}{1};

         if isempty(anti_powerkc{iMouse,iDay}) 
             continue
         else
         max_anti_powerkc(iMouse,iDay) = max(anti_powerkc{iMouse,iDay});
         max_anti_powerneg(iMouse,iDay) = max(anti_powerneg{iMouse,iDay});
         max_anti_powerrand(iMouse,iDay) = max(anti_powerrand{iMouse,iDay});
         end
         
    end
    
end


%% compute anti-cofiring power - v2

for iMouse=1:length(mouselist)
    
        for iDay = days{iMouse}
            
            n_removed{iMouse,iDay} = length(cells_removed{iMouse,iDay}{2});


            all_cells_anti{iMouse,iDay} = anticofiring_fields{iMouse,iDay}{1};
            all_cells_anti_id{iMouse,iDay} = idx_sorted{iMouse,iDay}{1};

            
            neg_anti_idx{iMouse,iDay} = all_cells_anti_id{iMouse,iDay}(n_removed{iMouse,iDay}:end);
            
            [~,~,real_id_anti_v2{iMouse,iDay}] = intersect(all_cells_anti{iMouse,iDay},cell_list{iMouse,iDay}{idiscard}','stable');

             non_anti{iMouse,iDay} = 0.0001*ones(ncells{iMouse,iDay}{1} - length(anti_power{iMouse,iDay}{1}),1);

            
        end
        
end


% %% plotting anti-cofiring power
% 
% for iMouse=1:length(mouselist)
%     
%     figure()
%     
%     t1{iMouse} = tiledlayout(3,3);
%     
%     for iDay = days
%         
%         non_anti{iMouse,iDay}{idiscard} = zeros(ncells{iMouse,iDay}(1) - length(anti_power{iMouse,iDay}{1}),1);
%         anti_power_withall{iMouse,iDay}{idiscard} = [anti_power{iMouse,iDay}{1}; non_anti{iMouse,iDay}{1}];
%         
%         nexttile
%         title(iDay)
%         
%         hold on
%         histfit(anti_power_withall{iMouse,iDay}{1},40,'kernel')
%         ax = gca
%         ax.LineWidth = 4
%         set(gca,'fontname','times');
%         set(gca,'Fontsize',20);
%         ax.FontSize = 20;
%         xlabel('Anti-cofiring power [a.u.]')
%         ylabel('Count')
%         hold off
%         xlim([0 inf])
%             
%     end
% end
% 
% 
% %% save all figures
% 
% if save_figures == 1
%     
%     for iMouse=1:length(mouselist)
%         
%         titlestr = strcat(mouselist(iMouse),'.pdf')
%         exportgraphics(t1{iMouse},titlestr)
%         
%     end
%     
% end

%% all the cells

Apower_allcells = cell2mat(cellflat([anti_powerkc non_anti])');
Apower_allcells = sort(Apower_allcells);

Apower_neg = cell2mat(cellflat([anti_powerneg non_anti])');
Apower_neg = sort(Apower_neg);

Apower_randPos = cell2mat(cellflat([anti_powerrand non_anti])');
Apower_randPos = sort(Apower_randPos);


figure()
histogram(log(Apower_allcells),40,'Normalization','probability')


figure()
pd = fitdist(log(Apower_allcells),'kernel')
y = pdf(pd,log(Apower_allcells));
yyaxis left
hold on
plot(log(Apower_allcells),y,'LineWidth',4,'Color','red')
yyaxis right
histogram(log(Apower_allcells),100,'Normalization','probability')
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',20);
ax.FontSize = 20;
xlabel('Anti-cofiring power [a.u.]')
ylabel('Probability')
hold off
xlim([0 inf])


%% looking for  binormal distribution 

x3 = Apower_allcells;

[f_allcells, x_allcells] = ksdensity(Apower_allcells);

pdf_normmixture = @(x3,p,mu1,mu2,sigma1,sigma2) ...
    0.5*p*gampdf(x3,mu1,sigma1) + normpdf(x3,mu2,sigma2);

pStart = .5;
muStart = [0.3 0.1]; %quantile(x3,[.25 .75])

sigmaStart = [0.12 0.1]; %sqrt(var(x3) - .25*diff(muStart).^2)

start = [pStart muStart sigmaStart];

lb = [0 -Inf -Inf 0 0];
ub = [1 Inf Inf Inf Inf];
paramEsts = mle(x3,'pdf',pdf_normmixture,'Start',start, ...
    'LowerBound',lb,'UpperBound',ub)

statset('mlecustom')

options = statset('MaxIter',3000,'MaxFunEvals',6000);
paramEsts = mle(x3,'pdf',pdf_normmixture,'Start',start, ...
    'LowerBound',lb,'UpperBound',ub,'Options',options)

figure()
histogram(x3,'Normalization','pdf')
hold on
xgrid = linspace(1.1*min(x3),1.1*max(x3),200);
pdfgrid = pdf_normmixture(xgrid, ...
    paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
plot(xgrid,pdfgrid,'-','LineWidth',5)
hold off
xlabel('x3')
ylabel('Probability')
legend('Sample Data','Fitted pdf','Location','best')
ylim([0 40])


%%

x3_2 = Apower_neg;

pdf_normmixture = @(x3_2,p,mu1,mu2,sigma1,sigma2) ...
    gampdf(x3_2,mu1,sigma1) + normpdf(x3_2,mu2,sigma2);

pStart = .1;
muStart = quantile(x3_2,[.25 .75])

sigmaStart = sqrt(var(x3_2) - .25*diff(muStart).^2)

start = [pStart muStart sigmaStart sigmaStart];

lb = [0 -Inf -Inf 0 0];
ub = [1 Inf Inf Inf Inf];
paramEsts = mle(x3_2,'pdf',pdf_normmixture,'Start',start, ...
    'LowerBound',lb,'UpperBound',ub)

statset('mlecustom')

options = statset('MaxIter',3000,'MaxFunEvals',6000);
paramEsts = mle(x3_2,'pdf',pdf_normmixture,'Start',start, ...
    'LowerBound',lb,'UpperBound',ub,'Options',options)

figure()
histogram(x3_2,'Normalization','pdf')
hold on
xgrid = linspace(1.1*min(x3_2),1.1*max(x3_2),200);
pdfgrid = pdf_normmixture(xgrid, ...
    paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
plot(xgrid,pdfgrid,'-','LineWidth',5)
hold off
xlabel('x3_2')
ylabel('Probability')
legend('Sample Data','Fitted pdf','Location','best')
ylim([0 40])

%%

x3_3 = Apower_randPos;

pdf_normmixture = @(x3_3,p,mu1,mu2,sigma1,sigma2) ...
    p*normpdf(x3_3,mu1,sigma1) + (1-p)*normpdf(x3_3,mu2,sigma2);

pStart = .5;
muStart = quantile(x3_3,[.25 .75])

sigmaStart = sqrt(var(x3_3) - .25*diff(muStart).^2)

start = [pStart muStart sigmaStart sigmaStart];

lb = [0 -Inf -Inf 0 0];
ub = [1 Inf Inf Inf Inf];
paramEsts = mle(x3_3,'pdf',pdf_normmixture,'Start',start, ...
    'LowerBound',lb,'UpperBound',ub)

statset('mlecustom')

options = statset('MaxIter',3000,'MaxFunEvals',6000);
paramEsts = mle(x3_3,'pdf',pdf_normmixture,'Start',start, ...
    'LowerBound',lb,'UpperBound',ub,'Options',options)

figure()
histogram(x3_3,'Normalization','pdf')
hold on
xgrid = linspace(1.1*min(x3_3),1.1*max(x3_3),200);
pdfgrid = pdf_normmixture(xgrid, ...
    paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
plot(xgrid,pdfgrid,'-','LineWidth',5)
hold off
xlabel('x3_3')
ylabel('Probability')
legend('Sample Data','Fitted pdf','Location','best')
ylim([0 40])


% %%
% x1 = Apower_allcells;
% 
% pf_truncpoiss = @(x1,lambda) poisspdf(x1,lambda)./(1-poisscdf(0,lambda));
% start = mean(x1)
% 
% [lambdaHat,lambdaCI] = mle(x1,'pdf',pf_truncpoiss,'Start',start, ...
%     'LowerBound',0)
% 
% [lambdaHat2,lambdaCI2] = mle(x1,'Distribution','Poisson', ...
%     'TruncationBounds',[0 Inf])
% 
% avar = mlecov(lambdaHat,x1,'pdf',pf_truncpoiss);
% stderr = sqrt(avar)
% 
% histogram(x1,'Normalization','pdf')
% xgrid = min(x1):max(x1);
% pmfgrid = pf_truncpoiss(xgrid,lambdaHat);
% hold on
% plot(xgrid,pmfgrid,'-')
% xlabel('x1')
% ylabel('Probability')
% legend('Sample Data','Fitted pmf','Location','best')
% hold off