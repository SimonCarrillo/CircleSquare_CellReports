clc
clear all

set(0,'DefaultFigureWindowStyle','docked')

%% introduce variables

residual_normalized_cyl1 = cell(6,9);
residual_normalized_cyl2 = cell(6,9);
residual_normalized_rct1 = cell(6,9);
residual_normalized_rct2 = cell(6,9);
residual_normalized_hmc1 = cell(6,9);
residual_normalized_hmc2 = cell(6,9);

derivative_residual_cyl1 = cell(6,9);
derivative_residual_cyl2 = cell(6,9);
derivative_residual_rct1 = cell(6,9);
derivative_residual_rct2 = cell(6,9);
derivative_residual_hmc1 = cell(6,9);
derivative_residual_hmc2 = cell(6,9);

residual_normalized_cyl1(:) = {nan(1,12)};
residual_normalized_cyl2(:) = {nan(1,12)};
residual_normalized_rct1(:) = {nan(1,12)};
residual_normalized_rct2(:) = {nan(1,12)};
residual_normalized_hmc1(:) = {nan(1,12)};
residual_normalized_hmc2(:) = {nan(1,12)};

derivative_residual_cyl1(:) = {nan(1,11)};
derivative_residual_cyl2(:) = {nan(1,11)};
derivative_residual_rct1(:) = {nan(1,11)};
derivative_residual_rct2(:) = {nan(1,11)};
derivative_residual_hmc1(:) = {nan(1,11)};
derivative_residual_hmc2(:) = {nan(1,11)};

%% experiment to analyze

mouselist = ["M29" "M34" "M39" "M35" "M20" "M19"];
ratiolist = ["8" "8" "8" "8" "8" "8"]; % "4" "4" "4" "4" "4" "4"

save_output = 0 ; % 1 for saving

for iMouse = 1:length(mouselist)
    
    animal = mouselist(iMouse); % mouse number
    
    time_bin = '1000'; % in ms
    expt = 'CircleSquare'; % experiment setup
    loadStr = 'TempCorr'; % 'Zlinear' or 'TempCorr'
    
    save_output = 0; % 1 for saving angles and proportion of overlap
    
    discard_tech= ["kc" "neg" "randPos"];
    
    for idiscard = 1:1
        
        discardPop = discard_tech(idiscard);
        
        ratio = ratiolist(iMouse);
        remove_outliers = 0; %1 for removing outliers
        alpharadius = Inf; % for alphashape
        ISO_output = 1; % 1 for plot
        
        path_ISO = strcat('/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/Response_Jan2022/Analysis_codes/',animal,'/','ISO_','CellDiscarded_',discardPop,'Ratio',ratio,loadStr,time_bin,expt,'_',animal,'_','.mat');    load(path_ISO);
        
        %% retrieve daily data
        %%% remove first week of data
        daystoremove = 0;
        
        %days = linspace(1,size(ISO,1),size(ISO,1));
        
        days{iMouse} = linspace(1+daystoremove,size(ISO,1),size(ISO,1)-daystoremove);
        
        for iDay = days{iMouse}
            
            ISOmap{iMouse,iDay} = ISO{iDay,2}; % extract dim reduced ISOMAP
            spikes{iMouse,iDay} = ISO{iDay,1}; % extract dim reduced ISOMAP
            sessList{iMouse,iDay} = ISO{iDay,4}; % extract session information
            sessStrList{iMouse,iDay} = ISO{iDay,5}; % extract session information
            
        end
        
        for iDay = days{iMouse}
            
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
            
            
            %%%%%%%%%%
            %% CYL1
            
            index = 0;
            index_timestep = 0;
            
            if size(CYL1_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                for timestep = time_windows
                    
                    
                    index_timestep = index_timestep+1;
                    
                    
                    index = 0;
                    for iter = timestep/2:4:size(CYL1_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [~,~,total_residual_cyl1{iMouse,iDay}(index,index_timestep)] = fitNormal(CYL1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        
                    end
                    
                    number_cyl1{iMouse,iDay}(index_timestep) = nnz(total_residual_cyl1{iMouse,iDay}(:,index_timestep));
                    residual_normalized_cyl1{iMouse,iDay}(index_timestep) = sum(total_residual_cyl1{iMouse,iDay}(:,index_timestep))/number_cyl1{iMouse,iDay}(index_timestep);
                    timestamp_save_cyl1{iMouse,iDay}(index_timestep) = timestep;
                end
                
                
            end
            
            %% CYL2
            
            index = 0;
            index_timestep = 0;
            
            if size(CYL2_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                for timestep = time_windows
                    
                    
                    index_timestep = index_timestep+1;
                    
                    
                    index = 0;
                    for iter = timestep/2:4:size(CYL2_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [~,~,total_residual_cyl2{iMouse,iDay}(index,index_timestep)] = fitNormal(CYL2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        
                    end
                    
                    number_cyl2{iMouse,iDay}(index_timestep) = nnz(total_residual_cyl2{iMouse,iDay}(:,index_timestep));
                    residual_normalized_cyl2{iMouse,iDay}(index_timestep) = sum(total_residual_cyl2{iMouse,iDay}(:,index_timestep))/number_cyl2{iMouse,iDay}(index_timestep);
                    timestamp_save_cyl2{iMouse,iDay}(index_timestep) = timestep;
                end
                
                
            end
            
            %% RCT1
            
            index = 0;
            index_timestep = 0;
            
            if size(RCT1_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                for timestep = time_windows
                    
                    
                    index_timestep = index_timestep+1;
                    
                    
                    index = 0;
                    for iter = timestep/2:4:size(RCT1_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [~,~,total_residual_rct1{iMouse,iDay}(index,index_timestep)] = fitNormal(RCT1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        
                    end
                    
                    number_rct1{iMouse,iDay}(index_timestep) = nnz(total_residual_rct1{iMouse,iDay}(:,index_timestep));
                    residual_normalized_rct1{iMouse,iDay}(index_timestep) = sum(total_residual_rct1{iMouse,iDay}(:,index_timestep))/number_rct1{iMouse,iDay}(index_timestep);
                    timestamp_save_rct1{iMouse,iDay}(index_timestep) = timestep;
                end
                
                
            end
            
            %% RCT2
            
            index = 0;
            index_timestep = 0;
            
            if size(RCT2_tmp{iMouse,iDay},1) < max(time_windows)
                
                continue
                
            else
                
                for timestep = time_windows
                    
                    
                    index_timestep = index_timestep+1;
                    
                    
                    index = 0;
                    for iter = timestep/2:4:size(RCT2_tmp{iMouse,iDay},1)-timestep/2
                        index = index+1;
                        
                        [~,~,total_residual_rct2{iMouse,iDay}(index,index_timestep)] = fitNormal(RCT2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
                        
                        
                    end
                    
                    number_rct2{iMouse,iDay}(index_timestep) = nnz(total_residual_rct2{iMouse,iDay}(:,index_timestep));
                    residual_normalized_rct2{iMouse,iDay}(index_timestep) = sum(total_residual_rct2{iMouse,iDay}(:,index_timestep))/number_rct2{iMouse,iDay}(index_timestep)';
                    timestamp_save_rct2{iMouse,iDay}(index_timestep) = timestep;
                end
                
                
            end
            %%%%%%
            
%              %% HMC1
%             
%             index = 0;
%             index_timestep = 0;
%             
%             if size(HMC1_tmp{iMouse,iDay},1) < max(time_windows)
%                 
%                 continue
%                 
%             else
%                 
%                 for timestep = time_windows
%                     
%                     
%                     index_timestep = index_timestep+1;
%                     
%                     
%                     index = 0;
%                     for iter = timestep/2:4:size(HMC1_tmp{iMouse,iDay},1)-timestep/2
%                         index = index+1;
%                         
%                         [~,~,total_residual_hmc1{iMouse,iDay}(index,index_timestep)] = fitNormal(HMC1_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
%                         
%                         
%                     end
%                     
%                     number_hmc1{iMouse,iDay}(index_timestep) = nnz(total_residual_hmc1{iMouse,iDay}(:,index_timestep));
%                     residual_normalized_hmc1{iMouse,iDay}(index_timestep) = sum(total_residual_hmc1{iMouse,iDay}(:,index_timestep))/number_hmc1{iMouse,iDay}(index_timestep);
%                     timestamp_save_hmc1{iMouse,iDay}(index_timestep) = timestep;
%                 end
%                 
%                 
%             end
%             
%             %% HMC2
%             
%             index = 0;
%             index_timestep = 0;
%             
%             if size(HMC2_tmp{iMouse,iDay},1) < max(time_windows)
%                 
%                 continue
%                 
%             else
%                 
%                 for timestep = time_windows
%                     
%                     
%                     index_timestep = index_timestep+1;
%                     
%                     
%                     index = 0;
%                     for iter = timestep/2:4:size(HMC2_tmp{iMouse,iDay},1)-timestep/2
%                         index = index+1;
%                         
%                         [~,~,total_residual_hmc2{iMouse,iDay}(index,index_timestep)] = fitNormal(HMC2_tmp{iMouse,iDay}(iter-(timestep/2-1):iter+(timestep/2),:),0,'g');
%                         
%                         
%                     end
%                     
%                     number_hmc2{iMouse,iDay}(index_timestep) = nnz(total_residual_hmc2{iMouse,iDay}(:,index_timestep));
%                     residual_normalized_hmc2{iMouse,iDay}(index_timestep) = sum(total_residual_hmc2{iMouse,iDay}(:,index_timestep))/number_hmc2{iMouse,iDay}(index_timestep)';
%                     timestamp_save_hmc2{iMouse,iDay}(index_timestep) = timestep;
%                 end
%                 
%                 
%             end
%             
            
            
            
            
            
            derivative_residual_cyl1{iMouse,iDay} = diff(residual_normalized_cyl1{iMouse,iDay})./diff(timestamp_save_cyl1{iMouse,iDay});
            derivative_residual_rct1{iMouse,iDay} = diff(residual_normalized_rct1{iMouse,iDay})./diff(timestamp_save_rct1{iMouse,iDay});
            derivative_residual_cyl2{iMouse,iDay} = diff(residual_normalized_cyl2{iMouse,iDay})./diff(timestamp_save_cyl2{iMouse,iDay});
            derivative_residual_rct2{iMouse,iDay} = diff(residual_normalized_rct2{iMouse,iDay})./diff(timestamp_save_rct2{iMouse,iDay});
%             derivative_residual_hmc1{iMouse,iDay} = diff(residual_normalized_hmc1{iMouse,iDay})./diff(timestamp_save_hmc1{iMouse,iDay});
%             derivative_residual_hmc2{iMouse,iDay} = diff(residual_normalized_hmc2{iMouse,iDay})./diff(timestamp_save_hmc2{iMouse,iDay});
%                         
            
            
            clear sess_Str
            
        end
    end 
end

%% unpack value per mouse for every time window

cyl1 = nan(6,length(time_windows));
cyl2 = nan(6,length(time_windows));
rct1 = nan(6,length(time_windows));
rct2 = nan(6,length(time_windows));
% hmc1 = nan(6,length(time_windows));
% hmc2 = nan(6,length(time_windows));

tmp_cyl1 = cell2mat(cellflat(residual_normalized_cyl1)');
tmp_cyl2 = cell2mat(cellflat(residual_normalized_cyl2)');
tmp_rct1 = cell2mat(cellflat(residual_normalized_rct1)');
tmp_rct2 = cell2mat(cellflat(residual_normalized_rct2)');
% tmp_hmc1 = cell2mat(cellflat(residual_normalized_hmc1)');
% tmp_hmc2 = cell2mat(cellflat(residual_normalized_hmc2)');

tmp_dcyl1 = cell2mat(cellflat(derivative_residual_cyl1)');
tmp_dcyl2 = cell2mat(cellflat(derivative_residual_cyl2)');
tmp_drct1 = cell2mat(cellflat(derivative_residual_rct1)');
tmp_drct2 = cell2mat(cellflat(derivative_residual_rct2)');
% tmp_dhmc1 = cell2mat(cellflat(derivative_residual_hmc1)');
% tmp_dhmc2 = cell2mat(cellflat(derivative_residual_hmc2)');


for iMouse = 1:length(mouselist)
    
    iJump = (iMouse-1)*8+iMouse:iMouse*9;
        
    cyl1(iMouse,:) = mean(tmp_cyl1(iJump,:),'omitnan');
    rct1(iMouse,:) = mean(tmp_rct1(iJump,:),'omitnan');
    cyl2(iMouse,:) = mean(tmp_cyl2(iJump,:),'omitnan');
    rct2(iMouse,:) = mean(tmp_rct2(iJump,:),'omitnan');
%     hmc1(iMouse,:) = mean(tmp_hmc1(iJump,:),'omitnan');
%     hmc2(iMouse,:) = mean(tmp_hmc2(iJump,:),'omitnan');
    
    dcyl1(iMouse,:) = mean(tmp_dcyl1(iJump,:),'omitnan');
    drct1(iMouse,:) = mean(tmp_drct1(iJump,:),'omitnan');
    dcyl2(iMouse,:) = mean(tmp_dcyl2(iJump,:),'omitnan');
    drct2(iMouse,:) = mean(tmp_drct2(iJump,:),'omitnan');
%     dhmc1(iMouse,:) = mean(tmp_dhmc1(iJump,:),'omitnan');
%     dhmc2(iMouse,:) = mean(tmp_dhmc2(iJump,:),'omitnan');
    
end
    
total_cyl1 = mean(cyl1,'omitnan');
total_cyl2 = mean(cyl2,'omitnan');
total_rct1 = mean(rct1,'omitnan');
total_rct2 = mean(rct2,'omitnan');
% total_hmc1 = mean(rct1,'omitnan');
% total_hmc2 = mean(rct2,'omitnan');

total_cyl1_sem = std(cyl1,'omitnan')/sqrt(6);
total_cyl2_sem = std(cyl2,'omitnan')/sqrt(6);
total_rct1_sem = std(rct1,'omitnan')/sqrt(6);
total_rct2_sem = std(rct2,'omitnan')/sqrt(6);

total_dcyl1 = mean(dcyl1,'omitnan');
total_dcyl2 = mean(dcyl2,'omitnan');
total_drct1 = mean(drct1,'omitnan');
total_drct2 = mean(drct2,'omitnan');
% total_dhmc1 = mean(rct1,'omitnan');
% total_dhmc2 = mean(rct2,'omitnan');

%% figure

figure(30)
left_color = [0 0 0];
right_color = [0 0 0];
set(figure(30),'defaultAxesColorOrder',[left_color; right_color]);
hold on
yyaxis left
s = plot(time_windows,total_rct1,'-','LineWidth',3,'Color','[0.5 0.5 0.5]');
s.Color(4)=0.3;
s = plot(time_windows,total_rct2,'-','LineWidth',3,'Color','[0.5 0.5 0.5]');
s.Color(4)=0.7;
s = plot(time_windows,total_cyl1,'-','LineWidth',3,'Color','red');
s.Color(4)=0.3;
s = plot(time_windows,total_cyl2,'-','LineWidth',3,'Color','red');
s.Color(4)=0.7;
% s = plot(time_windows,total_hmc1,'-','LineWidth',3,'Color','green');
% s.Color(4)=0.3;
% s = plot(time_windows,total_hmc2,'-','LineWidth',3,'Color','green');
% s.Color(4)=0.7;
ylabel('Normalized residual')
yyaxis right
s2 = plot(time_windows(1:end-1),total_drct1,'--','LineWidth',3,'Color','[0.5 0.5 0.5]');
s2.Color(4)=0.3;
s2 = plot(time_windows(1:end-1),total_drct2,'--','LineWidth',3,'Color','[0.5 0.5 0.5]');
s2.Color(4)=0.7;
s2 = plot(time_windows(1:end-1),total_dcyl1,'--','LineWidth',3,'Color','red');
s2.Color(4)=0.3;
s2 = plot(time_windows(1:end-1),total_dcyl2,'--','LineWidth',3,'Color','red');
s2.Color(4)=0.7;
% s2 = plot(time_windows(1:end-1),total_dhmc1,'--','LineWidth',3,'Color','green');
% s2.Color(4)=0.3;
% s2 = plot(time_windows(1:end-1),total_dhmc2,'--','LineWidth',3,'Color','green');
% s2.Color(4)=0.7;
xlabel('Timestep in seconds')
ylabel('First derivative of the normalized residual')
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',20);
ax.FontSize = 20;
xlim([0 120])
ylim([0 2*10^-2])

%% print mean and sem for timewindow 60 s

mean_summary = [total_cyl1(9) total_cyl2(9) total_rct1(9) total_rct2(9)]

sem_summary = [total_cyl1_sem(9) total_cyl2_sem(9) total_rct1_sem(9) total_rct2_sem(9)]
