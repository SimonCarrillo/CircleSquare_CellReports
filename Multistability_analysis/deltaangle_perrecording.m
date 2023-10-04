clc
clear all

set(0,'DefaultFigureWindowStyle','docked')

mouselist = ["M29" "M34" "M39" "M35" "M20" "M19"];
ratiolist = ["8" "8" "8" "8" "8" "8"]; % "4" "4" "4" "4" "4" "4"

save_output = 1; % 1 for saving

angle_threshold = [2 5 10 15 20 30 40 50 60 90]; % in degrees

for iThreshold = 1:length(angle_threshold)
    
    for iMouse = 1:length(mouselist)
        
        animal = mouselist(iMouse); % mouse number
        %% experiment to analyze
        
        time_bin = '1000'; % in ms
        expt = 'CircleSquare'; % experiment setup
        loadStr = 'TempCorr'; % 'Zlinear' or 'TempCorr'
        
        discard_tech= ["kc" "neg" "randPos"];
        
        for idiscard = 1:1 %3
            
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
                if size(CYL1_tmp{iMouse,iDay},1) < max(time_windows)
                    
                    continue
                    
                else
                    
                    for iter = timestep/2:4:size(CYL1_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [n2{iMouse,iDay}{index},~,~] = fitNormal(CYL1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        if index > 1
                            
                            angle_cyl1{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n2{iMouse,iDay}{index},n2{iMouse,iDay}{index-1})));
                            
                        end
                        
                        
                        
                    end
                    
                end
                
                shifts_cyl1{iMouse,iDay}{idiscard} = sum(angle_cyl1{iMouse,iDay}{idiscard}>angle_threshold(iThreshold));
                
                
                %% CYL2
                
                index = 0;
                
                if size(CYL2_tmp{iMouse,iDay},1) < max(time_windows)
                    
                    continue
                    
                else
                    
                    for iter = timestep/2:4:size(CYL2_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [n2_2{iMouse,iDay}{index},~,~] = fitNormal(CYL2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        if index > 1
                            
                            angle_cyl2{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n2_2{iMouse,iDay}{index},n2_2{iMouse,iDay}{index-1})));
                            
                        end
                        
                        
                    end
                    
                end
                
                
                shifts_cyl2{iMouse,iDay}{idiscard} = sum(angle_cyl2{iMouse,iDay}{idiscard}>angle_threshold(iThreshold));
                
                
                %% RCT1
                
                index = 0;
                if size(RCT1_tmp{iMouse,iDay},1) < max(time_windows)
                    
                    continue
                    
                else
                    
                    for iter = timestep/2:4:size(RCT1_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [n3{iMouse,iDay}{index},~,~] = fitNormal(RCT1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        if index > 1
                            
                            angle_rct1{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n3{iMouse,iDay}{index},n3{iMouse,iDay}{index-1})));
                            
                        end
                        
                        
                        
                    end
                    
                end
                
                shifts_rct1{iMouse,iDay}{idiscard} = sum(angle_rct1{iMouse,iDay}{idiscard}>angle_threshold(iThreshold));
                
                
                %% RCT2
                
                index = 0;
                
                if size(RCT2_tmp{iMouse,iDay},1) < max(time_windows)
                    
                    continue
                    
                else
                    
                    for iter = timestep/2:4:size(RCT2_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [n3_2{iMouse,iDay}{index},~,~] = fitNormal(RCT2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        if index > 1
                            
                            angle_rct2{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n3_2{iMouse,iDay}{index},n3_2{iMouse,iDay}{index-1})));
                            
                        end
                        
                        
                    end
                    
                end
                
                
                shifts_rct2{iMouse,iDay}{idiscard} = sum(angle_rct2{iMouse,iDay}{idiscard}>angle_threshold(iThreshold));
                
                %% HMC1
                
                index = 0;
                if size(HMC1_tmp{iMouse,iDay},1) < max(time_windows)
                    
                    continue
                    
                else
                    
                    for iter = timestep/2:4:size(HMC1_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [n1{iMouse,iDay}{index},~,~] = fitNormal(HMC1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        if index > 1
                            
                            angle_hmc1{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n1{iMouse,iDay}{index},n1{iMouse,iDay}{index-1})));
                            
                        end
                        
                        
                        
                    end
                    
                end
                
                shifts_hmc1{iMouse,iDay}{idiscard} = sum(angle_hmc1{iMouse,iDay}{idiscard}>angle_threshold(iThreshold));
                
                
                %% HMC2
                
                index = 0;
                
                if size(HMC2_tmp{iMouse,iDay},1) < max(time_windows)
                    
                    continue
                    
                else
                    
                    for iter = timestep/2:4:size(HMC2_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [n1_2{iMouse,iDay}{index},~,~] = fitNormal(HMC2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        if index > 1
                            
                            angle_hmc2{iMouse,iDay}{idiscard}(index-1) = rad2deg(acos(dot(n1_2{iMouse,iDay}{index},n1_2{iMouse,iDay}{index-1})));
                            
                        end
                        
                        
                    end
                    
                end
                
                
                shifts_hmc2{iMouse,iDay}{idiscard} = sum(angle_hmc2{iMouse,iDay}{idiscard}>angle_threshold(iThreshold));
                
            end
            
            for iDay = days
                
                shifts_cyl1kc(iMouse,iDay) = shifts_cyl1{iMouse,iDay}{1};
                shifts_cyl2kc(iMouse,iDay) = shifts_cyl2{iMouse,iDay}{1};
                shifts_rct1kc(iMouse,iDay) = shifts_rct1{iMouse,iDay}{1};
                shifts_rct2kc(iMouse,iDay) = shifts_rct2{iMouse,iDay}{1};
                shifts_hmc1kc(iMouse,iDay) = shifts_hmc1{iMouse,iDay}{1};
                shifts_hmc2kc(iMouse,iDay) = shifts_hmc2{iMouse,iDay}{1};
                
            end
        end
    end
    
         shifts = [shifts_hmc1kc shifts_cyl1kc shifts_rct1kc shifts_cyl2kc shifts_rct2kc shifts_hmc2kc]';
    
    %% save angle information
    
    if save_output == 1
        
        path_angle = strcat('angles_all',num2str(angle_threshold(iThreshold)),'.mat');
        path_shifts = strcat('shifts_all',num2str(angle_threshold(iThreshold)),'.mat');
        save(path_angle,'angle_rct1', 'angle_rct2', 'angle_cyl1', 'angle_cyl2');
        save(path_shifts,'shifts_rct1kc', 'shifts_rct2kc', 'shifts_cyl1kc', 'shifts_cyl2kc', 'shifts_hmc1kc','shifts_hmc2kc', 'shifts');
        
        display('saving')
    end
    
    %%
    
end


