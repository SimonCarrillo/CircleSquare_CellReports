clc
clear all

set(0,'DefaultFigureWindowStyle','docked')

%% experiment to analyze

mouselist = ["M29" "M34" "M39" "M35" "M20" "M19"];
ratiolist = ["8" "8" "8" "8" "8" "8"];

save_output = 0; % 1 for saving files

for iMouse = 1:length(mouselist)
    
    animal = mouselist(iMouse); % mouse number
    
    time_bin = '1000'; % in ms
    expt = 'CircleSquare'; % experiment setup
    loadStr = 'TempCorr'; % 'Zlinear' or 'TempCorr'
    
    discard_tech= ["kc" "neg" "randPos"];
    
    for idiscard = 1:3
        
        discardPop = discard_tech(idiscard);
        
        ratio = ratiolist(iMouse);
        remove_outliers = 0; %1 for removing outliers
        alpharadius = Inf; % for alphashape
        ISO_output = 1; % 1 for plot
        
        path_ISO = strcat('/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/Response_Jan2022/Analysis_codes/Final_figures/Raw_data/Manifold_60D/','ISO_','60_','CellDiscarded_',discardPop,'Ratio',ratio,loadStr,time_bin,expt,'_',animal,'_','.mat');
        load(path_ISO);
        
        %% retrieve daily data
        
        days{iMouse} = linspace(1,size(ISO,1),size(ISO,1));
        
        for iDay = days{iMouse}
            
            ISOmap{iMouse,iDay}{idiscard} = ISO{iDay,2}; % extract dim reduced ISOMAP
            sessList{iMouse,iDay}{idiscard} = ISO{iDay,4}; % extract session information
            sessStrList{iMouse,iDay}{idiscard} = ISO{iDay,5}; % extract session information
            
        end
        
        for iDay = days{iMouse}
            
            display(animal)
            display(iDay)
            
            ISO_data{iMouse,iDay}{idiscard} = ISOmap{iMouse,iDay}{idiscard};
            sess{iMouse,iDay}{idiscard} = sessList{iMouse,iDay}{idiscard};
            sess_Str{iMouse,iDay}{idiscard} = sessStrList{iMouse,iDay}{idiscard};
            
            %% figures ISO
            
            if ISO_output == 1
                
                index1 = 0;
                index2 = 0;
                index3 = 0;
                
                for Ii = 1:length(ISO_data{iMouse,iDay}{idiscard}(:,1))
                    
                    if sess_Str{iMouse,iDay}{idiscard}(Ii,1) == 'H' && (sess_Str{iMouse,iDay}{idiscard}(Ii,4) == '1' || sess_Str{iMouse,iDay}{idiscard}(Ii,4) == '2')
                        
                        index1 = index1 + 1;
                        
                        HMC_store{iMouse,iDay}{idiscard}(index1,:,:,:) = ISO_data{iMouse,iDay}{idiscard}(Ii,:);
                        
                        
                    elseif sess_Str{iMouse,iDay}{idiscard}(Ii,1) == 'R' && (sess_Str{iMouse,iDay}{idiscard}(Ii,4) == '1' || sess_Str{iMouse,iDay}{idiscard}(Ii,4) == '2')
                        
                        index2 = index2 + 1;
                        
                        RCT_store{iMouse,iDay}{idiscard}(index2,:,:,:) = ISO_data{iMouse,iDay}{idiscard}(Ii,:);
                        
                        
                    elseif sess_Str{iMouse,iDay}{idiscard}(Ii,1) == 'C' && (sess_Str{iMouse,iDay}{idiscard}(Ii,4) == '1' || sess_Str{iMouse,iDay}{idiscard}(Ii,4) == '2')
                        
                        index3 = index3 + 1;
                        
                        CYL_store{iMouse,iDay}{idiscard}(index3,:,:,:) = ISO_data{iMouse,iDay}{idiscard}(Ii,:);
                        
                        
                    end
                    
                end
                
            end
            
            %% least square fit and plot normal vector to best fit plane
            
            %%% for plotting
            
            title_expt = strcat('DiscardedCells',discardPop,'Ratio',ratio,animal,'Day',int2str(iDay));
            
            %% best fit plane
            
            % remove outliers or not
            if remove_outliers == 1
                
                HMC{iMouse,iDay} = rmoutliers(HMC_store{iMouse,iDay},'quartiles');
                CYL{iMouse,iDay} = rmoutliers(CYL_store{iMouse,iDay},'quartiles');
                RCT{iMouse,iDay} = rmoutliers(RCT_store{iMouse,iDay},'quartiles');
                
            else
                
                HMC{iMouse,iDay}{idiscard} = HMC_store{iMouse,iDay}{idiscard};
                CYL{iMouse,iDay}{idiscard} = CYL_store{iMouse,iDay}{idiscard};
                RCT{iMouse,iDay}{idiscard} = RCT_store{iMouse,iDay}{idiscard};
                
            end
            
            
            CYL1_tmp{iMouse,iDay}{idiscard} = CYL{iMouse,iDay}{idiscard}(1:round(size(CYL{iMouse,iDay}{idiscard},1)/2),:);
            RCT1_tmp{iMouse,iDay}{idiscard} = RCT{iMouse,iDay}{idiscard}(1:round(size(RCT{iMouse,iDay}{idiscard},1)/2),:);
            CYL2_tmp{iMouse,iDay}{idiscard} = CYL{iMouse,iDay}{idiscard}(round(size(CYL{iMouse,iDay}{idiscard},1)/2)+1:end,:);
            RCT2_tmp{iMouse,iDay}{idiscard} = RCT{iMouse,iDay}{idiscard}(round(size(RCT{iMouse,iDay}{idiscard},1)/2)+1:end,:);
            HMC1_tmp{iMouse,iDay}{idiscard} = HMC{iMouse,iDay}{idiscard}(1:round(size(HMC{iMouse,iDay}{idiscard},1)/2),:);
            HMC2_tmp{iMouse,iDay}{idiscard} = HMC{iMouse,iDay}{idiscard}(round(size(HMC{iMouse,iDay}{idiscard},1)/2)+1:end,:);
            
            
            %% computation of radius of influence for 80 percent of data per environment
            %%% centroid of each environment
            
            centroid_CYL1{iMouse,iDay}{idiscard} = mean(CYL1_tmp{iMouse,iDay}{idiscard});
            centroid_CYL2{iMouse,iDay}{idiscard} = mean(CYL2_tmp{iMouse,iDay}{idiscard});
            centroid_RCT1{iMouse,iDay}{idiscard} = mean(RCT1_tmp{iMouse,iDay}{idiscard});
            centroid_RCT2{iMouse,iDay}{idiscard} = mean(RCT2_tmp{iMouse,iDay}{idiscard});
            centroid_HMC1{iMouse,iDay}{idiscard} = mean(HMC1_tmp{iMouse,iDay}{idiscard});
            centroid_HMC2{iMouse,iDay}{idiscard} = mean(HMC2_tmp{iMouse,iDay}{idiscard});
            
            for itime = 1:size(CYL1_tmp{iMouse,iDay}{idiscard},1)
                
                distance_tocentroid_CYL1{iMouse,iDay}{idiscard}(itime) = norm(centroid_CYL1{iMouse,iDay}{idiscard} - CYL1_tmp{iMouse,iDay}{idiscard}(itime,:));
                
            end
            
            for itime = 1:size(CYL2_tmp{iMouse,iDay}{idiscard},1)
                
                distance_tocentroid_CYL2{iMouse,iDay}{idiscard}(itime) = norm(centroid_CYL2{iMouse,iDay}{idiscard} - CYL2_tmp{iMouse,iDay}{idiscard}(itime,:));
                
            end
            
            for itime = 1:size(RCT1_tmp{iMouse,iDay}{idiscard},1)
                
                distance_tocentroid_RCT1{iMouse,iDay}{idiscard}(itime) = norm(centroid_RCT1{iMouse,iDay}{idiscard} - RCT1_tmp{iMouse,iDay}{idiscard}(itime,:));
                
            end
            
            for itime = 1:size(RCT2_tmp{iMouse,iDay}{idiscard},1)
                
                distance_tocentroid_RCT2{iMouse,iDay}{idiscard}(itime) = norm(centroid_RCT2{iMouse,iDay}{idiscard} - RCT2_tmp{iMouse,iDay}{idiscard}(itime,:));
            end
            
            
%            pts_CYL1{iMouse,iDay} = linspace(0,max([max(distance_tocentroid_CYL1{iMouse,iDay}{1}) max(distance_tocentroid_CYL2{iMouse,iDay}{1}) max(distance_tocentroid_RCT1{iMouse,iDay}{1}) max(distance_tocentroid_RCT2{iMouse,iDay}{1})])+10,100);
%            pts_CYL2{iMouse,iDay} = pts_CYL1{iMouse,iDay};
%            pts_RCT1{iMouse,iDay} = pts_CYL1{iMouse,iDay};
%            pts_RCT2{iMouse,iDay} = pts_CYL1{iMouse,iDay};
            pts_CYL1{iMouse,iDay} = linspace(0,max(distance_tocentroid_CYL1{iMouse,iDay}{1})+10,100);
            pts_CYL2{iMouse,iDay} = linspace(0,max(distance_tocentroid_CYL2{iMouse,iDay}{1})+10,100);
            pts_RCT1{iMouse,iDay} = linspace(0,max(distance_tocentroid_RCT1{iMouse,iDay}{1})+10,100);
            pts_RCT2{iMouse,iDay} = linspace(0,max(distance_tocentroid_RCT2{iMouse,iDay}{1})+10,100);
            
            [f_CYL1{iMouse,iDay}{idiscard},xi_CYL1{iMouse,iDay}{idiscard}] = ksdensity(distance_tocentroid_CYL1{iMouse,iDay}{idiscard},pts_CYL1{iMouse,iDay});
            [f_CYL2{iMouse,iDay}{idiscard},xi_CYL2{iMouse,iDay}{idiscard}] = ksdensity(distance_tocentroid_CYL2{iMouse,iDay}{idiscard},pts_CYL2{iMouse,iDay});
            [f_RCT1{iMouse,iDay}{idiscard},xi_RCT1{iMouse,iDay}{idiscard}] = ksdensity(distance_tocentroid_RCT1{iMouse,iDay}{idiscard},pts_RCT1{iMouse,iDay});
            [f_RCT2{iMouse,iDay}{idiscard},xi_RCT2{iMouse,iDay}{idiscard}] = ksdensity(distance_tocentroid_RCT2{iMouse,iDay}{idiscard},pts_RCT2{iMouse,iDay});
            
            
            [N_CYL1{iMouse,iDay}{idiscard}, pts_CYL1{iMouse,iDay}] = histcounts(distance_tocentroid_CYL1{iMouse,iDay}{idiscard},pts_CYL1{iMouse,iDay}, 'Normalization', 'probability');
            [N_CYL2{iMouse,iDay}{idiscard}, pts_CYL2{iMouse,iDay}] = histcounts(distance_tocentroid_CYL2{iMouse,iDay}{idiscard},pts_CYL2{iMouse,iDay}, 'Normalization', 'probability');
            [N_RCT1{iMouse,iDay}{idiscard}, pts_RCT1{iMouse,iDay}] = histcounts(distance_tocentroid_RCT1{iMouse,iDay}{idiscard},pts_RCT1{iMouse,iDay}, 'Normalization', 'probability');
            [N_RCT2{iMouse,iDay}{idiscard}, pts_RCT2{iMouse,iDay}] = histcounts(distance_tocentroid_RCT2{iMouse,iDay}{idiscard},pts_RCT2{iMouse,iDay}, 'Normalization', 'probability');

        end
        
        clear sess_Str
        
        
    end
    
    for iDay = days{iMouse}
        
        sk_allcells{iMouse,iDay} = [skewness(distance_tocentroid_CYL1{iMouse,iDay}{1}) skewness(distance_tocentroid_CYL2{iMouse,iDay}{1}) skewness(distance_tocentroid_RCT1{iMouse,iDay}{1}) skewness(distance_tocentroid_RCT2{iMouse,iDay}{1})];
        sk_neg{iMouse,iDay} = [skewness(distance_tocentroid_CYL1{iMouse,iDay}{2}) skewness(distance_tocentroid_CYL2{iMouse,iDay}{2}) skewness(distance_tocentroid_RCT1{iMouse,iDay}{2}) skewness(distance_tocentroid_RCT2{iMouse,iDay}{2})];
        sk_rand{iMouse,iDay} = [skewness(distance_tocentroid_CYL1{iMouse,iDay}{3}) skewness(distance_tocentroid_CYL2{iMouse,iDay}{3}) skewness(distance_tocentroid_RCT1{iMouse,iDay}{3}) skewness(distance_tocentroid_RCT2{iMouse,iDay}{3})];

        %%% replace zeros
        
        CYL1_KLD_mix(iMouse,iDay) = getKullbackLeibler(N_CYL1{iMouse,iDay}{2},N_CYL1{iMouse,iDay}{3});
        CYL2_KLD_mix(iMouse,iDay) = getKullbackLeibler(N_CYL2{iMouse,iDay}{2},N_CYL2{iMouse,iDay}{3});
        RCT1_KLD_mix(iMouse,iDay) = getKullbackLeibler(N_RCT1{iMouse,iDay}{2},N_RCT1{iMouse,iDay}{3});
        RCT2_KLD_mix(iMouse,iDay) = getKullbackLeibler(N_RCT2{iMouse,iDay}{2},N_RCT2{iMouse,iDay}{3});        
        
        CYL1_KLD_anti(iMouse,iDay) = getKullbackLeibler(N_CYL1{iMouse,iDay}{1},N_CYL1{iMouse,iDay}{2});
        CYL2_KLD_anti(iMouse,iDay) = getKullbackLeibler(N_CYL2{iMouse,iDay}{1},N_CYL2{iMouse,iDay}{2});
        RCT1_KLD_anti(iMouse,iDay) = getKullbackLeibler(N_RCT1{iMouse,iDay}{1},N_RCT1{iMouse,iDay}{2});
        RCT2_KLD_anti(iMouse,iDay) = getKullbackLeibler(N_RCT2{iMouse,iDay}{1},N_RCT2{iMouse,iDay}{2});

        CYL1_KLD_rand(iMouse,iDay) = getKullbackLeibler(N_CYL1{iMouse,iDay}{1},N_CYL1{iMouse,iDay}{3});
        CYL2_KLD_rand(iMouse,iDay) = getKullbackLeibler(N_CYL2{iMouse,iDay}{1},N_CYL2{iMouse,iDay}{3});
        RCT1_KLD_rand(iMouse,iDay) = getKullbackLeibler(N_RCT1{iMouse,iDay}{1},N_RCT1{iMouse,iDay}{3});
        RCT2_KLD_rand(iMouse,iDay) = getKullbackLeibler(N_RCT2{iMouse,iDay}{1},N_RCT2{iMouse,iDay}{3});

    end
    
end

sk_kc = cell2mat(cellflat(sk_allcells));
sk_neg = cell2mat(cellflat(sk_neg));
sk_rand = cell2mat(cellflat(sk_rand));

KLD_anti = [reshape(CYL1_KLD_anti,[],1); reshape(CYL2_KLD_anti,[],1); reshape(RCT1_KLD_anti,[],1); reshape(RCT2_KLD_anti,[],1)];
KLD_rand = [reshape(CYL1_KLD_rand,[],1); reshape(CYL2_KLD_rand,[],1); reshape(RCT1_KLD_rand,[],1); reshape(RCT2_KLD_rand,[],1)];
KLD_mix = [reshape(CYL1_KLD_mix,[],1); reshape(CYL2_KLD_mix,[],1); reshape(RCT1_KLD_mix,[],1); reshape(RCT2_KLD_mix,[],1)];

KLD_anti = KLD_anti(KLD_anti>0);
KLD_rand = KLD_rand(KLD_rand>0);
KLD_mix = KLD_mix(KLD_mix>0);

%% example figure
figure()

day_toplot = 7;
mouse_toplot = 3;

t1 =plot(xi_RCT2{mouse_toplot,day_toplot}{1},f_RCT2{mouse_toplot,day_toplot}{1},'-','Color','[0.5 0.5 0.5]','LineWidth',4);
hold on
t2 =plot(xi_RCT2{mouse_toplot,day_toplot}{2},f_RCT2{mouse_toplot,day_toplot}{2},'--','Color','[0.5 0.5 0.5]','LineWidth',4);
t3 =plot(xi_RCT2{mouse_toplot,day_toplot}{3},f_RCT2{mouse_toplot,day_toplot}{3},':','Color','[0.5 0.5 0.5]','LineWidth',4);
hold off
ylabel('Probabilty density')
xlabel('Distance [a.u.]')
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',40);
ax.FontSize = 40;
legend([t1 t2 t3],'All cells','Anti-cofiring cells removed','Random cells removed')
xlim([0 Inf])

display(skewness(distance_tocentroid_RCT2{mouse_toplot,day_toplot}{1}))
display(skewness(distance_tocentroid_RCT2{mouse_toplot,day_toplot}{2}))
display(skewness(distance_tocentroid_RCT2{mouse_toplot,day_toplot}{3}))

%% save files

if save_output == 1
    
    path_overlap = strcat('overlap','_60D',animal,'_ratio',ratio,'.mat');
    
    save(path_overlap,'proportion_cylrct','proportion_rctcyl')
    
end
