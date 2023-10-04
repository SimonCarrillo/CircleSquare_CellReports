clc
clear all

set(0,'DefaultFigureWindowStyle','docked')

mouselist = ["M29" "M34" "M39" "M35" "M20" "M19"];
ratiolist = ["8" "8" "8" "8" "8" "8"]; % "4" "4" "4" "4" "4" "4"

save_output = 1; % 1 for saving

%% introduction of variables

pvalkc = nan(6,9); kkc= nan(6,9); Kkc= nan(6,9);
pvalneg = nan(6,9); kneg= nan(6,9); Kneg= nan(6,9);
pvalrand = nan(6,9); krand= nan(6,9); Krand= nan(6,9);

pvalkc2 = nan(6,9); kkc2= nan(6,9); Kkc2= nan(6,9);
pvalneg2 = nan(6,9); kneg2= nan(6,9); Kneg2= nan(6,9);
pvalrand2 = nan(6,9); krand2= nan(6,9); Krand2= nan(6,9);

pvalkcsame = nan(6,9); kkcsame= nan(6,9); Kkcsame= nan(6,9);
pvalnegsame = nan(6,9); knegsame= nan(6,9); Knegsame= nan(6,9);
pvalrandsame = nan(6,9); krandsame= nan(6,9); Krandsame= nan(6,9);

pvalkcsame2 = nan(6,9); kkcsame2= nan(6,9); Kkcsame2= nan(6,9);
pvalnegsame2 = nan(6,9); knegsame2= nan(6,9); Knegsame2= nan(6,9);
pvalrandsame2 = nan(6,9); krandsame2= nan(6,9); Krandsame2= nan(6,9);


for iMouse = 1:length(mouselist)
    
    animal = mouselist(iMouse); % mouse number
    %% experiment to analyze
    
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
        
        path_ISO = strcat('/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/Response_Jan2022/Analysis_codes/',animal,'/','ISO_','CellDiscarded_',discardPop,'Ratio',ratio,loadStr,time_bin,expt,'_',animal,'_','.mat');
        load(path_ISO);
        
        %% retrieve daily data
        
        %%% remove first week of data
        daystoremove = 0;
        
        %days = linspace(1,size(ISO,1),size(ISO,1));
        
        days = linspace(1+daystoremove,size(ISO,1),size(ISO,1)-daystoremove);
        
        for iDay = days
            
            ISOmap{iMouse,iDay} = ISO{iDay,2}; % extract dim reduced ISOMAP
            spikes{iMouse,iDay} = ISO{iDay,1}; % extract dim reduced ISOMAP
            sessList{iMouse,iDay} = ISO{iDay,4}; % extract session information
            sessStrList{iMouse,iDay} = ISO{iDay,5}; % extract session information
            
        end
        
        for iDay = days
            
            ISO_data = ISOmap{iMouse,iDay};
            sess = sessList{iMouse,iDay};
            sess_Str = sessStrList{iMouse,iDay};
            
            %% figures ISO
            
            if ISO_output == 1
                
                index1 = 0;
                index2 = 0;
                index3 = 0;
                
                for Ii = 1:length(ISO_data(:,1))
                    
                    if sess_Str(Ii,1) == 'H' && (sess_Str(Ii,4) == '1' || sess_Str(Ii,4) == '2')
                        
                        index1 = index1 + 1;
                        
                        HMC_store{iMouse,iDay}(index1,:,:,:) = ISO_data(Ii,:);
                        
                        
                    elseif sess_Str(Ii,1) == 'R' && (sess_Str(Ii,4) == '1' || sess_Str(Ii,4) == '2')
                        
                        index2 = index2 + 1;
                        
                        RCT_store{iMouse,iDay}(index2,:,:,:) = ISO_data(Ii,:);
                        
                        
                    elseif sess_Str(Ii,1) == 'C' && (sess_Str(Ii,4) == '1' || sess_Str(Ii,4) == '2')
                        
                        index3 = index3 + 1;
                        
                        CYL_store{iMouse,iDay}(index3,:,:,:) = ISO_data(Ii,:);
                        
                        
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
                
                HMC{iMouse,iDay} = HMC_store{iMouse,iDay};
                CYL{iMouse,iDay} = CYL_store{iMouse,iDay};
                RCT{iMouse,iDay} = RCT_store{iMouse,iDay};
                
            end
            
            
            CYL1_tmp{iMouse,iDay} = CYL{iMouse,iDay}(1:round(size(CYL{iMouse,iDay},1)/2),:);
            RCT1_tmp{iMouse,iDay} = RCT{iMouse,iDay}(1:round(size(RCT{iMouse,iDay},1)/2),:);
            CYL2_tmp{iMouse,iDay} = CYL{iMouse,iDay}(round(size(CYL{iMouse,iDay},1)/2)+1:end,:);
            RCT2_tmp{iMouse,iDay} = RCT{iMouse,iDay}(round(size(RCT{iMouse,iDay},1)/2)+1:end,:);
            HMC1_tmp{iMouse,iDay} = HMC{iMouse,iDay}(1:round(size(HMC{iMouse,iDay},1)/2),:);
            HMC2_tmp{iMouse,iDay} = HMC{iMouse,iDay}(round(size(HMC{iMouse,iDay},1)/2)+1:end,:);
            
            
            %%%%%%%%%
            time_windows = 60;
            timestep = time_windows;
            
            %%%%%%%%%%
            %% CYL1
            
            index = 0;
            if size(CYL1_tmp{iMouse,iDay},1) < 240%max(time_windows)
                
                continue
                
            else
                
                for iter = timestep/2:4:size(CYL1_tmp{iMouse,iDay},1)-timestep/2
                    index = index+1;
                    
                    [n2{iMouse,iDay}{index},~,~] = fitNormal(CYL1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                    
                    if index > 1
                        
                        angle_cyl1{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n2{iMouse,iDay}{index},n2{iMouse,iDay}{1})));
                        
                    end
                    
                    
                    
                end
                
            end
            
            %% CYL2
            
            index = 0;
            
            if size(CYL2_tmp{iMouse,iDay},1) < 240%max(time_windows)
                
                continue
                
            else
                
                for iter = timestep/2:4:size(CYL2_tmp{iMouse,iDay},1)-timestep/2
                    index = index+1;
                    
                    [n2_2{iMouse,iDay}{index},~,~] = fitNormal(CYL2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                    
                    if index > 1
                        
                        angle_cyl2{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n2_2{iMouse,iDay}{index},n2_2{iMouse,iDay}{1})));
                        
                    end
                    
                    angle_cyl1cyl2{iMouse,iDay}{idiscard}(index) = rad2deg(acos(dot(n2_2{iMouse,iDay}{index},n2{iMouse,iDay}{1})));
                    
                    
                end
                
            end
            
            %% RCT1
            
            index = 0;
            if size(RCT1_tmp{iMouse,iDay},1) < 240%max(time_windows)
                
                continue
                
            else
                
                for iter = timestep/2:4:size(RCT1_tmp{iMouse,iDay},1)-timestep/2
                    index = index+1;
                    
                    [n3{iMouse,iDay}{index},~,~] = fitNormal(RCT1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                    
                    if index > 1
                        
                        angle_rct1{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n3{iMouse,iDay}{index},n3{iMouse,iDay}{1})));
                        
                    end
                    
                    angle_rct1cyl1{iMouse,iDay}{idiscard}(index) = rad2deg(acos(dot(n3{iMouse,iDay}{index},n2{iMouse,iDay}{1})));
                    
                end
                
            end
            
            %% RCT2
            
            index = 0;
            
            if size(RCT2_tmp{iMouse,iDay},1) < 240%max(time_windows)
                
                continue
                
            else
                
                for iter = timestep/2:4:size(RCT2_tmp{iMouse,iDay},1)-timestep/2
                    index = index+1;
                    
                    [n3_2{iMouse,iDay}{index},~,~] = fitNormal(RCT2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                    
                    if index > 1
                        
                        angle_rct2{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n3_2{iMouse,iDay}{index},n3_2{iMouse,iDay}{1})));
                        
                    end
                    
                    angle_rct2cyl2{iMouse,iDay}{idiscard}(index) = rad2deg(acos(dot(n3_2{iMouse,iDay}{index},n2_2{iMouse,iDay}{1})));
                    angle_rct1rct2{iMouse,iDay}{idiscard}(index) = rad2deg(acos(dot(n3_2{iMouse,iDay}{index},n3{iMouse,iDay}{1})));
                    angle_rct2cyl1{iMouse,iDay}{idiscard}(index) = rad2deg(acos(dot(n3_2{iMouse,iDay}{index},n2{iMouse,iDay}{1})));
                    
                    
                end
                
            end
            
            %%% summary angles
            
            RCT1RCT2{iMouse,iDay}{idiscard} = [angle_rct1rct2{iMouse,iDay}{idiscard} angle_rct1{iMouse,iDay}{idiscard}];
            CYL1CYL2{iMouse,iDay}{idiscard} = [angle_cyl1cyl2{iMouse,iDay}{idiscard} angle_cyl1{iMouse,iDay}{idiscard}];
            RCT1CYL1{iMouse,iDay}{idiscard} = [angle_rct1cyl1{iMouse,iDay}{idiscard} angle_cyl1{iMouse,iDay}{idiscard}];
            RCT2CYL2{iMouse,iDay}{idiscard} = [angle_rct2cyl2{iMouse,iDay}{idiscard} angle_cyl2{iMouse,iDay}{idiscard}];
            
        end
        
    end
    
    for iDay = days
        
        if isempty(angle_cyl1{iMouse,iDay}) || isempty(angle_rct1cyl1{iMouse,iDay}) || isempty(angle_rct2cyl2{iMouse,iDay}) || isempty(angle_cyl2{iMouse,iDay}) || isempty(angle_rct1{iMouse,iDay})
            continue
        else
            %% kuipers test
            
            [pvalkc(iMouse,iDay), kkc(iMouse,iDay), Kkc(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_rct1cyl1{iMouse,iDay}{1}), deg2rad(angle_cyl1{iMouse,iDay}{1}),100); % deg2rad(angle_cyl1{2,5}{1}), 100, 1)
            [pvalneg(iMouse,iDay), kneg(iMouse,iDay), Kneg(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_rct1cyl1{iMouse,iDay}{2}), deg2rad(angle_cyl1{iMouse,iDay}{2}), 100);
            [pvalrand(iMouse,iDay), krand(iMouse,iDay), Krand(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_rct1cyl1{iMouse,iDay}{3}), deg2rad(angle_cyl1{iMouse,iDay}{3}), 100);
            [pvalkc2(iMouse,iDay), kkc2(iMouse,iDay), Kkc2(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_rct2cyl2{iMouse,iDay}{1}), deg2rad(angle_cyl2{iMouse,iDay}{1}),100); % deg2rad(angle_cyl1{2,5}{1}), 100, 1)
            [pvalneg2(iMouse,iDay), kneg2(iMouse,iDay), Kneg2(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_rct2cyl2{iMouse,iDay}{2}), deg2rad(angle_cyl2{iMouse,iDay}{2}), 100);
            [pvalrand2(iMouse,iDay), krand2(iMouse,iDay), Krand2(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_rct2cyl2{iMouse,iDay}{3}), deg2rad(angle_cyl2{iMouse,iDay}{3}), 100);
            
            
            [pvalkcsame(iMouse,iDay), kkcsame(iMouse,iDay), Kkcsame(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_cyl1cyl2{iMouse,iDay}{1}), deg2rad(angle_cyl1{iMouse,iDay}{1}),100); % deg2rad(angle_cyl1{2,5}{1}), 100, 1)
            [pvalkcsame2(iMouse,iDay), kkcsame2(iMouse,iDay), Kkcsame2(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_rct1rct2{iMouse,iDay}{1}), deg2rad(angle_rct1{iMouse,iDay}{1}),100); % deg2rad(angle_cyl1{2,5}{1}), 100, 1)
            [pvalnegsame(iMouse,iDay), knegsame(iMouse,iDay), Knegsame(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_cyl1cyl2{iMouse,iDay}{2}), deg2rad(angle_cyl1{iMouse,iDay}{2}),100); % deg2rad(angle_cyl1{2,5}{1}), 100, 1)
            [pvalnegsame2(iMouse,iDay), knegsame2(iMouse,iDay), Knegsame2(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_rct1rct2{iMouse,iDay}{2}), deg2rad(angle_rct1{iMouse,iDay}{2}),100); % deg2rad(angle_cyl1{2,5}{1}), 100, 1)
            [pvalrandsame(iMouse,iDay), krandsame(iMouse,iDay), Krandsame(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_cyl1cyl2{iMouse,iDay}{3}), deg2rad(angle_cyl1{iMouse,iDay}{3}),100); % deg2rad(angle_cyl1{2,5}{1}), 100, 1)
            [pvalrandsame2(iMouse,iDay), krandsame2(iMouse,iDay), Krandsame2(iMouse,iDay)] = circ_kuipertest(deg2rad(angle_rct1rct2{iMouse,iDay}{3}), deg2rad(angle_rct1{iMouse,iDay}{3}),100); % deg2rad(angle_cyl1{2,5}{1}), 100, 1)
        end
    end
    
end


%% summmary statistics

statistics_neg = reshape(kneg,[],1);
statistics_kc = reshape(kkc,[],1);
statistics_rand = reshape(krand,[],1);


statistics_neg2 = reshape(kneg2,[],1);
statistics_kc2 = reshape(kkc2,[],1);
statistics_rand2 = reshape(krand2,[],1);

% remove zeros
ind_zeros = find(statistics_neg>0);
statistics_neg = statistics_neg(ind_zeros);
statistics_kc = statistics_kc(ind_zeros);
statistics_rand = statistics_rand(ind_zeros);

statistics_neg2 = statistics_neg2(ind_zeros);
statistics_kc2 = statistics_kc2(ind_zeros);
statistics_rand2 = statistics_rand2(ind_zeros);

%% combine diff environments

kc_all = [statistics_kc; statistics_kc2];
neg_all = [statistics_neg; statistics_neg2];
rand_all = [statistics_rand; statistics_rand2];

%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine same
statistics_sameneg = reshape(knegsame,[],1);
statistics_samekc = reshape(kkcsame,[],1);
statistics_samerand = reshape(krandsame,[],1);


statistics_sameneg2 = reshape(knegsame2,[],1);
statistics_samekc2 = reshape(kkcsame2,[],1);
statistics_samerand2 = reshape(krandsame2,[],1);

% remove zeros
ind_zeros = find(statistics_sameneg>0);
statistics_sameneg = statistics_sameneg(ind_zeros);
statistics_samekc = statistics_samekc(ind_zeros);
statistics_samerand = statistics_samerand(ind_zeros);

statistics_sameneg2 = statistics_sameneg2(ind_zeros);
statistics_samekc2 = statistics_samekc2(ind_zeros);
statistics_samerand2 = statistics_samerand2(ind_zeros);


kcsame_all = [statistics_samekc; statistics_samekc2];
negsame_all = [statistics_sameneg; statistics_sameneg2];
randsame_all = [statistics_samerand; statistics_samerand2];

%% save angle information

if save_output == 1
    
    path_angle = strcat('angles_all','.mat');
    save(path_angle,'angle_rct1cyl1','angle_rct2cyl2', 'angle_rct2cyl1', 'angle_rct1', 'angle_rct2', 'angle_cyl1', 'angle_cyl2');
    
end

%% range of K critical

K_tmp = [Kkc Kkc2 Kneg Kneg2 Krand Krand2 Kkcsame Kkcsame2 Knegsame Knegsame2 Krandsame Krandsame2];

min_K = min(min(K_tmp(K_tmp>0)));
max_K = max(max(K_tmp(K_tmp>0)));