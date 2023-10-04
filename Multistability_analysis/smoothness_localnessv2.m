clc
clear all

set(0,'DefaultFigureWindowStyle','docked')

%% experiment to analyze
mouselist = ["M29" "M34" "M39" "M35" "M20" "M19"];
ratiolist = ["8" "8" "8" "8" "8" "8"]; % "4" "4" "4" "4" "4" "4"

for iMouse = 1:length(mouselist)
    
    animal = mouselist(iMouse); % mouse number
    
    time_bin = '1000'; % in ms
    expt = 'CircleSquare'; % experiment setup
    loadStr = 'TempCorr'; % 'Zlinear' or 'TempCorr'
    
    discard_tech= ["kc" "neg" "randPos"];
    
    for idiscard = 1:1
        
        discardPop = discard_tech(idiscard);
        
        ratio = ratiolist(iMouse);
        remove_outliers = 0; %1 for removing outliers
        alpharadius = Inf; % for alphashape
        ISO_output = 1; % 1 for plot
        
        path_ISO = strcat('/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/Response_Jan2022/Analysis_codes/',animal,'/','ISO_','CellDiscarded_',discardPop,'Ratio',ratio,loadStr,time_bin,expt,'_',animal,'_','.mat');
        load(path_ISO);
        
        %% retrieve daily data
        
        days = linspace(1,size(ISO,1),size(ISO,1));
        
        for iDay = days
            
            ISOmap{iMouse,iDay} = ISO{iDay,2}; % extract dim reduced ISOMAP
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
            time_windows = [4 6 8 10 16 20 30 40 60 90 120 240];
            time= linspace(1,max(time_windows)-6,max(time_windows)-6);
            max_size = max(time_windows);
            
            %% compute localness for CYL1
            
            if size(CYL1_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                for iTime = 1:length(time)
                    
                    iwant = time(randperm(length(time),1));
                    
                    cyl_localness_v2{iMouse,iDay}(iTime) = norm(CYL1_tmp{iMouse,iDay}(iTime+1,:)-CYL1_tmp{iMouse,iDay}(iTime,:));
                    cyl_localness_v2_rand{iMouse,iDay}(iTime) = norm(CYL1_tmp{iMouse,iDay}(iTime+1,:)-CYL1_tmp{iMouse,iDay}(iwant,:));
                    
                    
                end
            end
            %% compute localness for CYL2
            
            
            if size(CYL2_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                for iTime = 1:length(time)
                    
                    iwant = time(randperm(length(time),1));
                    
                    cyl2_localness_v2{iMouse,iDay}(iTime) = norm(CYL2_tmp{iMouse,iDay}(iTime+1,:)-CYL2_tmp{iMouse,iDay}(iTime,:));
                    cyl2_localness_v2_rand{iMouse,iDay}(iTime) = norm(CYL2_tmp{iMouse,iDay}(iTime+1,:)-CYL2_tmp{iMouse,iDay}(iwant,:));
                    
                    
                end
                
            end
            
            %% compute localness for RCT1
            
            
            if size(RCT1_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                for iTime = 1:length(time)
                    
                    iwant = time(randperm(length(time),1));
                    
                    
                    rct_localness_v2{iMouse,iDay}(iTime) = norm(RCT1_tmp{iMouse,iDay}(iTime+1,:)-RCT1_tmp{iMouse,iDay}(iTime,:));
                    rct_localness_v2_rand{iMouse,iDay}(iTime) = norm(RCT1_tmp{iMouse,iDay}(iTime+1,:)-RCT1_tmp{iMouse,iDay}(iwant,:));
                    
                    
                end
                
            end
            
            %% compute localness for RCT2
            
            if size(RCT2_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                for iTime = 1:length(time)
                    
                    iwant = time(randperm(length(time),1));
                    
                    
                    rct2_localness_v2{iMouse,iDay}(iTime) = norm(RCT2_tmp{iMouse,iDay}(iTime+1,:)-RCT2_tmp{iMouse,iDay}(iTime,:));
                    rct2_localness_v2_rand{iMouse,iDay}(iTime) = norm(RCT2_tmp{iMouse,iDay}(iTime+1,:)-RCT2_tmp{iMouse,iDay}(iwant,:));
                    
                end
                
            end
            
            %% compute smoothness per environment and mouse
            
            timestep = 60;
            index = 0;
            
            %% compute smoothness for CYL1
            
            if size(CYL1_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                index = 0;
                for iter = timestep/2:4:max_size-timestep/2
                    
                    index = index + 1;
                    
                    CYL1_random{iMouse,iDay}(:,1) = [min(CYL1_tmp{iMouse,iDay}(:,1)) + (max(CYL1_tmp{iMouse,iDay}(:,1))-min(CYL1_tmp{iMouse,iDay}(:,1))).*rand(size(CYL1_tmp{iMouse,iDay}(:,1),1),1)];
                    CYL1_random{iMouse,iDay}(:,2) = [min(CYL1_tmp{iMouse,iDay}(:,2)) + (max(CYL1_tmp{iMouse,iDay}(:,2))-min(CYL1_tmp{iMouse,iDay}(:,2))).*rand(size(CYL1_tmp{iMouse,iDay}(:,1),1),1)];
                    CYL1_random{iMouse,iDay}(:,3) = [min(CYL1_tmp{iMouse,iDay}(:,3)) + (max(CYL1_tmp{iMouse,iDay}(:,3))-min(CYL1_tmp{iMouse,iDay}(:,3))).*rand(size(CYL1_tmp{iMouse,iDay}(:,1),1),1)];
                    
                    
                    [n2{iMouse,iDay}{index},~,total_residual2{iMouse,iDay}(index)] = fitNormal(CYL1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+timestep/2,:),0,'g');
                    
                    [n2_rand{iMouse,iDay}{index},~,total_residual2_rand{iMouse,iDay}(index)] = fitNormal(CYL1_random{iMouse,iDay}(iter-(timestep/2-1):iter+timestep/2,:),0,'g');
                    
                    if index > 1
                        
                        angle_cyl1{iMouse,iDay}(index-1) = rad2deg(acos(dot(n2{iMouse,iDay}{index},n2{iMouse,iDay}{index-1})));
                        angle_cyl1_rand{iMouse,iDay}(index-1) = rad2deg(acos(dot(n2_rand{iMouse,iDay}{index},n2_rand{iMouse,iDay}{index-1})));
                        
                        angleabs_cyl1{iMouse,iDay}(index-1) = rad2deg(acos(dot(n2{iMouse,iDay}{index},n2{iMouse,iDay}{1})));
                        angleabs_cyl1_rand{iMouse,iDay}(index-1) = rad2deg(acos(dot(n2_rand{iMouse,iDay}{index},n2_rand{iMouse,iDay}{1})));
                        
                    end
                    
                end
                
            end
            
            %% compute smoothness for CYL2
            index = 0;
            if size(CYL2_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                index = 0;
                for iter = timestep/2:4:max_size-timestep/2
                    
                    index = index + 1;
                    
                    
                    CYL2_random{iMouse,iDay}(:,1) = [min(CYL2_tmp{iMouse,iDay}(:,1)) + (max(CYL2_tmp{iMouse,iDay}(:,1))-min(CYL2_tmp{iMouse,iDay}(:,1))).*rand(size(CYL2_tmp{iMouse,iDay}(:,1),1),1)];
                    CYL2_random{iMouse,iDay}(:,2) = [min(CYL2_tmp{iMouse,iDay}(:,2)) + (max(CYL2_tmp{iMouse,iDay}(:,2))-min(CYL2_tmp{iMouse,iDay}(:,2))).*rand(size(CYL2_tmp{iMouse,iDay}(:,1),1),1)];
                    CYL2_random{iMouse,iDay}(:,3) = [min(CYL2_tmp{iMouse,iDay}(:,3)) + (max(CYL2_tmp{iMouse,iDay}(:,3))-min(CYL2_tmp{iMouse,iDay}(:,3))).*rand(size(CYL2_tmp{iMouse,iDay}(:,1),1),1)];
                    
                    [n2_2{iMouse,iDay}{index},~,total_residual2_2{iMouse,iDay}(index)] = fitNormal(CYL2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+timestep/2,:),0,'g');
                    
                    [n2_2_rand{iMouse,iDay}{index},~,total_residual2_2_rand{iMouse,iDay}(index)] = fitNormal(CYL2_random{iMouse,iDay}(iter-(timestep/2-1):iter+timestep/2,:),0,'g');
                    
                    if index > 1
                        
                        angle_cyl2{iMouse,iDay}(index-1) = rad2deg(acos(dot(n2_2{iMouse,iDay}{index},n2_2{iMouse,iDay}{index-1})));
                        angle_cyl2_rand{iMouse,iDay}(index-1) = rad2deg(acos(dot(n2_2_rand{iMouse,iDay}{index},n2_2_rand{iMouse,iDay}{index-1})));
                        
                        angleabs_cyl2{iMouse,iDay}(index-1) = rad2deg(acos(dot(n2_2{iMouse,iDay}{index},n2_2{iMouse,iDay}{1})));
                        angleabs_cyl2_rand{iMouse,iDay}(index-1) = rad2deg(acos(dot(n2_2_rand{iMouse,iDay}{index},n2_2_rand{iMouse,iDay}{1})));
                                                
                    end
                    
                end
            end
            
            
            %% compute smoothness for RCT1
            index = 0;
            if size(RCT1_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                index = 0;
                for iter = timestep/2:4:max_size-timestep/2
                    
                    index = index + 1;
                    
                    RCT1_random{iMouse,iDay}(:,1) = [min(RCT1_tmp{iMouse,iDay}(:,1)) + (max(RCT1_tmp{iMouse,iDay}(:,1))-min(RCT1_tmp{iMouse,iDay}(:,1))).*rand(size(RCT1_tmp{iMouse,iDay}(:,1),1),1)];
                    RCT1_random{iMouse,iDay}(:,2) = [min(RCT1_tmp{iMouse,iDay}(:,2)) + (max(RCT1_tmp{iMouse,iDay}(:,2))-min(RCT1_tmp{iMouse,iDay}(:,2))).*rand(size(RCT1_tmp{iMouse,iDay}(:,1),1),1)];
                    RCT1_random{iMouse,iDay}(:,3) = [min(RCT1_tmp{iMouse,iDay}(:,3)) + (max(RCT1_tmp{iMouse,iDay}(:,3))-min(RCT1_tmp{iMouse,iDay}(:,3))).*rand(size(RCT1_tmp{iMouse,iDay}(:,1),1),1)];
                    
                    [n3{iMouse,iDay}{index},~,total_residual3{iMouse,iDay}(index)] = fitNormal(RCT1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+timestep/2,:),0,'b');
                    
                    [n3_rand{iMouse,iDay}{index},~,total_residual3_rand{iMouse,iDay}(index)] = fitNormal(RCT1_random{iMouse,iDay}(iter-(timestep/2-1):iter+timestep/2,:),0,'b');
                    
                    
                    if index > 1
                        
                        angle_rct1{iMouse,iDay}(index-1) = rad2deg(acos(dot(n3{iMouse,iDay}{index},n3{iMouse,iDay}{index-1})));
                        
                        angle_rct1_rand{iMouse,iDay}(index-1) = rad2deg(acos(dot(n3_rand{iMouse,iDay}{index},n3_rand{iMouse,iDay}{index-1})));
                        
                        angleabs_rct1{iMouse,iDay}(index-1) = rad2deg(acos(dot(n3{iMouse,iDay}{index},n2{iMouse,iDay}{1})));
                        angleabs_rct1_rand{iMouse,iDay}(index-1) = rad2deg(acos(dot(n3_rand{iMouse,iDay}{index},n2_rand{iMouse,iDay}{1})));
                        
                    end
                    
                end
                
            end
            
            %% compute smoothness for RCT2
            index = 0;
            if size(RCT2_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                index = 0;
                for iter = timestep/2:4:max_size-timestep/2
                    
                    index = index + 1;
                    
                    RCT2_random{iMouse,iDay}(:,1) = [min(RCT2_tmp{iMouse,iDay}(:,1)) + (max(RCT2_tmp{iMouse,iDay}(:,1))-min(RCT2_tmp{iMouse,iDay}(:,1))).*rand(size(RCT2_tmp{iMouse,iDay}(:,1),1),1)];
                    RCT2_random{iMouse,iDay}(:,2) = [min(RCT2_tmp{iMouse,iDay}(:,2)) + (max(RCT2_tmp{iMouse,iDay}(:,2))-min(RCT2_tmp{iMouse,iDay}(:,2))).*rand(size(RCT2_tmp{iMouse,iDay}(:,1),1),1)];
                    RCT2_random{iMouse,iDay}(:,3) = [min(RCT2_tmp{iMouse,iDay}(:,3)) + (max(RCT2_tmp{iMouse,iDay}(:,3))-min(RCT2_tmp{iMouse,iDay}(:,3))).*rand(size(RCT2_tmp{iMouse,iDay}(:,1),1),1)];
                    
                    
                    [n3_2{iMouse,iDay}{index},~,total_residual3_2{iMouse,iDay}(index)] = fitNormal(RCT2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+timestep/2,:),0,'b');
                    
                    [n3_2_rand{iMouse,iDay}{index},~,total_residual3_2_rand{iMouse,iDay}(index)] = fitNormal(RCT2_random{iMouse,iDay}(iter-(timestep/2-1):iter+timestep/2,:),0,'b');
                    
                    if index > 1
                        
                        angle_rct2{iMouse,iDay}(index-1) = rad2deg(acos(dot(n3_2{iMouse,iDay}{index},n3_2{iMouse,iDay}{index-1})));
                        angle_rct2_rand{iMouse,iDay}(index-1) = rad2deg(acos(dot(n3_2_rand{iMouse,iDay}{index},n3_2_rand{iMouse,iDay}{index-1})));
                        
                        angleabs_rct2{iMouse,iDay}(index-1) = rad2deg(acos(dot(n3_2{iMouse,iDay}{index},n3_2{iMouse,iDay}{1})));
                        angleabs_rct2_rand{iMouse,iDay}(index-1) = rad2deg(acos(dot(n3_2_rand{iMouse,iDay}{index},n3_2_rand{iMouse,iDay}{1})));
                                                
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

%% summmary

localness = [cell2mat(cellflat(cyl_localness_v2)'); cell2mat(cellflat(cyl2_localness_v2)'); cell2mat(cellflat(rct_localness_v2)'); cell2mat(cellflat(rct2_localness_v2)')];
smoothness = [cell2mat(cellflat(angle_cyl1)'); cell2mat(cellflat(angle_cyl2)'); cell2mat(cellflat(angle_rct1)'); cell2mat(cellflat(angle_rct2)')];

localness_rand = [cell2mat(cellflat(cyl_localness_v2_rand)'); cell2mat(cellflat(cyl2_localness_v2_rand)'); cell2mat(cellflat(rct_localness_v2_rand)'); cell2mat(cellflat(rct2_localness_v2_rand)')];
smoothness_rand = [cell2mat(cellflat(angle_cyl1_rand)'); cell2mat(cellflat(angle_cyl2_rand)'); cell2mat(cellflat(angle_rct1_rand)'); cell2mat(cellflat(angle_rct2_rand)')];


smoothness_abs = [cell2mat(cellflat(angleabs_cyl1)'); cell2mat(cellflat(angleabs_cyl2)'); cell2mat(cellflat(angleabs_rct1)'); cell2mat(cellflat(angleabs_rct2)')];
smoothness_rand_abs = [cell2mat(cellflat(angleabs_cyl1_rand)'); cell2mat(cellflat(angleabs_cyl2_rand)'); cell2mat(cellflat(angleabs_rct1_rand)'); cell2mat(cellflat(angleabs_rct2_rand)')];

%% figure localness

figure()
t1 = stdshade(localness,0.2,'red',time)
hold on
t2 = stdshade(localness_rand,0.2,'black',time)
xlabel('Time (s)')
ylabel('Localness (a.u.)')
legend([t1 t2],'Experimental','Random')
ax = gca
ax.LineWidth = 10
set(gca,'fontname','times');
set(gca,'Fontsize',80);
ax.FontSize = 80;
legend off
%%%%%%%%%%%%%%%%%%%

%% smoothness computation

time_2 = linspace(timestep/2,max_size+timestep/2,size(smoothness,2));

figure()
t1 = stdshade(smoothness,0.2,'red',time_2)
hold on
t2 = stdshade(smoothness_rand,0.2,'black',time_2)
xlabel('Time (s)')
ylabel('Smoothness ($^o$)','interpreter','latex')
legend([t1 t2],'Experimental','Random')
ax = gca
ax.LineWidth = 10
set(gca,'fontname','times');
set(gca,'Fontsize',80);
ax.FontSize = 80;
legend off

%% smoothness computation - absolute value

time_2 = linspace(timestep/2,max_size+timestep/2,size(smoothness,2));

figure()
t1 = stdshade(smoothness_abs,0.2,'red',time_2)
hold on
t2 = stdshade(smoothness_rand_abs,0.2,'black',time_2)
xlabel('Time (s)')
ylabel('Smoothness ($^o$)','interpreter','latex')
legend([t1 t2],'Experimental','Random')
ax = gca
ax.LineWidth = 10
set(gca,'fontname','times');
set(gca,'Fontsize',80);
ax.FontSize = 80;
legend off

%% example individual localness

figure()
plot(time,localness(50,:),'color','red','LineWidth',4)
hold on
plot(time,localness_rand(50,:),'color','black','LineWidth',3)

xlabel('Time (s)')
ylabel('Localness (a.u.)')
legend([t1 t2],'Experimental','Random')
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',40);
ax.FontSize = 40;
legend off

%% example individual smoothness

figure()
plot(time_2,angleabs_cyl1{3,5},'color','red','LineWidth',4)
hold on
plot(time_2,angleabs_cyl1_rand{3,5},'--','color','red','LineWidth',2)
plot(time_2+max(time_2)+timestep/2,angleabs_rct1{3,5},'color','[0.5 0.5 0.5]','LineWidth',4)
plot(time_2+max(time_2)+timestep/2,angleabs_rct1_rand{3,5},'--','color','[0.5 0.5 0.5]','LineWidth',2)
xlabel('Time (s)')
ylabel('Smoothness ($^o$)','interpreter','latex')
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',40);
ax.FontSize = 40;

%% stats summary

localness_flatten = reshape(localness,[],1);
mean_local = mean(localness_flatten)
sem_local = std(localness_flatten)/sqrt(length(localness_flatten))

localness_flatten_rand = reshape(localness_rand,[],1);
mean_local_rand = mean(localness_flatten_rand)
sem_local_rand = std(localness_flatten_rand)/sqrt(length(localness_flatten_rand))


smoothness_flatten = reshape(smoothness,[],1);
mean_smooth = mean(smoothness_flatten)
sem_smooth = std(smoothness_flatten)/sqrt(length(smoothness_flatten))


smoothness_flatten_rand = reshape(smoothness_rand,[],1);
mean_smooth_rand = mean(smoothness_flatten_rand)
sem_smooth_rand = std(smoothness_flatten_rand)/sqrt(length(smoothness_flatten_rand))

