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
    %animal = 'M39'; % mouse number
    loadStr = 'TempCorr'; % 'Zlinear' or 'TempCorr'
    
    save_output = 0; % 1 for saving angles and proportion of overlap
    
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
            
            
            CYL1{iMouse,iDay} = CYL{iMouse,iDay}(1:round(size(CYL{iMouse,iDay},1)/2),:);
            RCT1{iMouse,iDay} = RCT{iMouse,iDay}(1:round(size(RCT{iMouse,iDay},1)/2),:);
            CYL2{iMouse,iDay} = CYL{iMouse,iDay}(round(size(CYL{iMouse,iDay},1)/2)+1:end,:);
            RCT2{iMouse,iDay} = RCT{iMouse,iDay}(round(size(RCT{iMouse,iDay},1)/2)+1:end,:);
            HMC1{iMouse,iDay} = HMC{iMouse,iDay}(1:round(size(HMC{iMouse,iDay},1)/2),:);
            HMC2{iMouse,iDay} = HMC{iMouse,iDay}(round(size(HMC{iMouse,iDay},1)/2)+1:end,:);
            
            
            %% Participation ratio CYL
            
            CNeuron = cov(CYL1{iMouse,iDay});
            
            PRNumer = sum((diag(CNeuron))).^2;
            
            
            for i = 1:size(CNeuron,1)
                PRDenom1(i,i) = (CNeuron(i,i)).^2;
            end
            
            PRDenom = sum(PRDenom1,'all');
            PR_cyl1(iMouse,iDay) = PRNumer/PRDenom
            
            clear PRDenom1 PRNumer CNeuron PRDenom
            
            %% Participation ratio CYL
            
            CNeuron = cov(CYL2{iMouse,iDay});
            
            PRNumer = sum((diag(CNeuron))).^2;
            
            
            for i = 1:size(CNeuron,1)
                PRDenom1(i,i) = (CNeuron(i,i)).^2;
            end
            
            PRDenom = sum(PRDenom1,'all');
            PR_cyl2(iMouse,iDay) = PRNumer/PRDenom
            
            clear PRDenom1 PRNumer CNeuron PRDenom
            %% Participation ratio RCT
            
            CNeuron = cov(RCT1{iMouse,iDay});
            
            PRNumer = sum((diag(CNeuron))).^2;
            
            
            for i = 1:size(CNeuron,1)
                PRDenom1(i,i) = (CNeuron(i,i)).^2;
            end
            
            PRDenom = sum(PRDenom1,'all');
            PR_rct1(iMouse,iDay)  = PRNumer/PRDenom
            
            clear PRDenom1 PRNumer CNeuron PRDenom
            
            %% Participation ratio RCT
            
            CNeuron = cov(RCT2{iMouse,iDay});
            
            PRNumer = sum((diag(CNeuron))).^2;
            
            
            for i = 1:size(CNeuron,1)
                PRDenom1(i,i) = (CNeuron(i,i)).^2;
            end
            
            PRDenom = sum(PRDenom1,'all');
            PR_rct2(iMouse,iDay)  = PRNumer/PRDenom
            
            clear PRDenom1 PRNumer CNeuron PRDenom
            
            %% Participation ratio HMC
            
            CNeuron = cov(HMC1{iMouse,iDay});
            PRNumer = sum((diag(CNeuron))).^2;
            
            for i = 1:size(CNeuron,1)
                PRDenom1(i,i) = (CNeuron(i,i)).^2;
            end
            
            PRDenom = sum(PRDenom1,'all');
            PR_hmc1(iMouse,iDay) = PRNumer/PRDenom
            
            clear PRDenom1 PRNumer CNeuron PRDenom
            
            %% Participation ratio HMC
            
            CNeuron = cov(HMC2{iMouse,iDay});
            PRNumer = sum((diag(CNeuron))).^2;
            
            for i = 1:size(CNeuron,1)
                PRDenom1(i,i) = (CNeuron(i,i)).^2;
            end
            
            PRDenom = sum(PRDenom1,'all');
            PR_hmc2(iMouse,iDay) = PRNumer/PRDenom
            
            clear PRDenom1 PRNumer CNeuron PRDenom
            
            %% Participation ratio CYL - spikes
            
            CNeuron = cov(spikes{1,6}(1:300,:));
            PRNumer = sum((diag(CNeuron))).^2;
            
            for i = 1:size(CNeuron,1)
                PRDenom1(i,i) = (CNeuron(i,i)).^2;
            end
            
            PRDenom = sum(PRDenom1,'all');
            PR_spikes = PRNumer/PRDenom
            
            clear PRDenom1 PRNumer CNeuron PRDenom
        end
        
    end
    
end


HMC1_PR = reshape(PR_hmc1,[],1);
RCT1_PR = reshape(PR_rct1,[],1);
CYL1_PR = reshape(PR_cyl1,[],1);
HMC2_PR = reshape(PR_hmc2,[],1);
RCT2_PR = reshape(PR_rct2,[],1);
CYL2_PR = reshape(PR_cyl2,[],1);

ind_zeros = find(HMC1_PR>0);
HMC_PR = [HMC1_PR(ind_zeros); HMC2_PR(ind_zeros)];
RCT_PR = [RCT1_PR(ind_zeros); RCT2_PR(ind_zeros)] ;
CYL_PR = [CYL1_PR(ind_zeros); CYL2_PR(ind_zeros)] ;


% %% account for 50% of the data
% 
% [HMC_PR,TF_hmc] = rmoutliers(HMC_PR,'percentile',[25 75]);
% [CYL_PR,TF_cyl] = rmoutliers(CYL_PR,'percentile',[25 75]);
% [RCT_PR,TF_rct] = rmoutliers(RCT_PR,'percentile',[25 75]);


%%%%%%%%%%%%%%%%%%%%%%%%
%% control group I

for iter = 1:100
    
    t = 0:1:300;
    xt1 = t+rand;
    yt1 = 2*t;
    zt1 = 5*t;
    
    xt1rand = t+0.1*rand(size(t));
    yt1rand = 2*t+0.1*rand(size(t));
    zt1rand = 5*t+0.1*rand(size(t));
    
    %%% Participation ratio = 1
    
    PR1 = [xt1rand; yt1rand; zt1rand];
    
    CNeuron = cov(PR1');
    PRNumer = sum((diag(CNeuron))).^2;
    
    for i = 1:size(CNeuron,1)
        PRDenom1(i,i) = (CNeuron(i,i)).^2;
    end
    
    PRDenom = sum(PRDenom1,'all');
    PR_control1(iter)= PRNumer/PRDenom
    
    clear PRDenom1 PRNumer CNeuron PRDenom
    %% control group II
    
    n3 = [0 0 1];
    n1 = [0 1 0.5];
    n2 = [1 0 0.5];
    
    x3 = linspace(min(RCT{iMouse,iDay}(:,1)),max(RCT{iMouse,iDay}(:,1)),15);
    y3 = linspace(min(RCT{iMouse,iDay}(:,2)),max(RCT{iMouse,iDay}(:,2)),15);
    [xg3, yg3] = meshgrid(x3, y3);
    z3 = -n3(1)/n3(3)*(xg3)-n3(2)/n3(3)*(yg3);
    
    xnew = reshape(xg3,[],1).*rand(225,1);
    ynew = reshape(yg3,[],1).*rand(225,1);
    znew = reshape(z3,[],1).*rand(225,1);
    
    x2 = linspace(min(RCT{iMouse,iDay}(:,1)),max(RCT{iMouse,iDay}(:,1)),15);
    y2 = linspace(min(RCT{iMouse,iDay}(:,2)),max(RCT{iMouse,iDay}(:,2)),15);
    [xg2, yg2] = meshgrid(x2, y2);
    z2 = -n2(1)/n2(3)*(xg2)-n2(2)/n2(3)*(yg2);
    
    xnew2 = reshape(xg2,[],1).*rand(225,1);
    ynew2 = reshape(yg2,[],1).*rand(225,1);
    znew2 = reshape(z2,[],1).*rand(225,1);
    
    x1 = linspace(min(RCT{iMouse,iDay}(:,1)),max(RCT{iMouse,iDay}(:,1)),15);
    y1 = linspace(min(RCT{iMouse,iDay}(:,2)),max(RCT{iMouse,iDay}(:,2)),15);
    [xg1, yg1] = meshgrid(x1, y1);
    z1 = -n1(1)/n1(3)*(xg1)-n1(2)/n1(3)*(yg1);
    
    xnew1 = reshape(xg1,[],1).*rand(225,1);
    ynew1 = reshape(yg1,[],1).*rand(225,1);
    znew1 = reshape(z1,[],1).*rand(225,1);
    
    %     figure(10)
    %     surf(xg3,yg3,z3)
    %     hold on
    %     surf(xg2,yg2,z2)
    %     surf(xg1,yg1,z1)
    %     view(3)
    %     scatter3(xnew(1:3:end),ynew(1:3:end),znew(1:3:end))
    %     scatter3(xnew2(1:3:end),ynew2(1:3:end),znew2(1:3:end))
    %     scatter3(xnew1(1:3:end),ynew1(1:3:end),znew1(1:3:end))
    %
    %%% Participation ratio = 1
    
    PR2 = [xnew ynew znew]';
    PR2_2 = [[xnew(1:2:end);xnew2(1:2:end)] [ynew(1:2:end);ynew2(1:2:end)] [znew(1:2:end);znew2(1:2:end)]]';
    PR2_3 = [[xnew(1:3:end);xnew2(1:3:end);xnew1(1:3:end)] [ynew(1:3:end);ynew2(1:3:end);ynew1(1:3:end)] [znew(1:3:end);znew2(1:3:end);znew1(1:3:end)]]';
    
    CNeuron = cov(PR2');
    PRNumer = sum((diag(CNeuron))).^2;
    
    for i = 1:size(CNeuron,1)
        PRDenom1(i,i) = (CNeuron(i,i)).^2;
    end
    
    PRDenom = sum(PRDenom1,'all');
    PR_control2(iter)= PRNumer/PRDenom
    
    clear PRDenom1 PRNumer CNeuron PRDenom
    
    CNeuron = cov(PR2_2');
    PRNumer = sum((diag(CNeuron))).^2;
    
    for i = 1:size(CNeuron,1)
        PRDenom1(i,i) = (CNeuron(i,i)).^2;
    end
    
    PRDenom = sum(PRDenom1,'all');
    PR_control2_2(iter)= PRNumer/PRDenom
    
    clear PRDenom1 PRNumer CNeuron PRDenom
    
    CNeuron = cov(PR2_3');
    PRNumer = sum((diag(CNeuron))).^2;
    
    for i = 1:size(CNeuron,1)
        PRDenom1(i,i) = (CNeuron(i,i)).^2;
    end
    
    PRDenom = sum(PRDenom1,'all');
    PR_control2_3(iter)= PRNumer/PRDenom
    
    clear PRDenom1 PRNumer CNeuron PRDenom
    %% control group III
    
    RCT1_random(:,1) = [5*rand(300,1)];
    RCT1_random(:,2) = [5*rand(300,1)];
    RCT1_random(:,3) = [5*rand(300,1)];
    
    %%% Participation ratio = 3
    
    PR3 = RCT1_random';
    
    %     figure()
    %     scatter3(RCT1_random(:,1), RCT1_random(:,2), RCT1_random(:,3))
    %
    CNeuron = cov(PR3');
    PRNumer = sum((diag(CNeuron))).^2;
    
    for i = 1:size(CNeuron,1)
        PRDenom1(i,i) = (CNeuron(i,i)).^2;
    end
    
    PRDenom = sum(PRDenom1,'all');
    PR_control3(iter)= PRNumer/PRDenom
    
    clear PRDenom1 PRNumer CNeuron PRDenom
    
    %% control group swissroll
    
    [swissroll, ~, ~] = generate_data('swiss',1000,0.1);
    
    %%% Participation ratio = 3
    
    CNeuron = cov(swissroll);
    PRNumer = sum((diag(CNeuron))).^2;
    
    for i = 1:size(CNeuron,1)
        PRDenom1(i,i) = (CNeuron(i,i)).^2;
    end
    
    PRDenom = sum(PRDenom1,'all');
    PR_control_swiss(iter)= PRNumer/PRDenom
    
    clear PRDenom1 PRNumer CNeuron PRDenom
    
    
end

%% testing color codes scatter

angle_threshold = [2 5 10 15 20 30 40 50 60 90]; % in degrees

for iThreshold = 1:length(angle_threshold)
    
    path_load = strcat('shifts_all',num2str(angle_threshold(iThreshold)),'.mat');
    load(path_load);
    
    cyl1 = reshape(PR_cyl1,[],1);
    cyl1_color = reshape(shifts_cyl1kc,[],1);
    cyl2 = reshape(PR_cyl2,[],1);
    cyl2_color = reshape(shifts_cyl2kc,[],1);
    rct1 = reshape(PR_rct1,[],1);
    rct1_color = reshape(shifts_rct1kc,[],1);
    rct2 = reshape(PR_rct2,[],1);
    rct2_color = reshape(shifts_rct2kc,[],1);
    hmc1 = reshape(PR_hmc1,[],1);
    hmc1_color = reshape(shifts_hmc1kc,[],1);
    hmc2 = reshape(PR_hmc2,[],1);
    hmc2_color = reshape(shifts_hmc2kc,[],1);
    
    cyl_data = [cyl1; cyl2];
    cyl_colordata = [cyl1_color; cyl2_color];
    rct_data = [rct1; rct2];
    rct_colordata = [rct1_color; rct2_color];
    hmc_data = [hmc1; hmc2];
    hmc_colordata = [hmc1_color; hmc2_color];
    
    ind_zeros = find(cyl_data>0);
    
    cyl_data = cyl_data(ind_zeros);
    cyl_colordata = cyl_colordata(ind_zeros);
    rct_data = rct_data(ind_zeros);
    rct_colordata = rct_colordata(ind_zeros);
    hmc_data = hmc_data(ind_zeros);
    hmc_colordata = hmc_colordata(ind_zeros);
   
%         %% remove 50% "outliers"
%         
%     cyl_data = cyl_data(~TF_cyl);
%     cyl_colordata = cyl_colordata(~TF_cyl);
%     rct_data = rct_data(~TF_rct);
%     rct_colordata = rct_colordata(~TF_rct);
%     hmc_data = hmc_data(~TF_hmc);
%     hmc_colordata = hmc_colordata(~TF_hmc);
%     
    
    total = [cyl_data rct_data hmc_data];
    total_color = [cyl_colordata rct_colordata hmc_colordata];
    
    figure()
    h = boxplot(total,'Notch','off','Labels',{'CYL','RCT','HMC'},'Whisker',1,'MedianStyle','line');
    hold on
    s_cyl = scatter(1*ones(size(total(:,1))).*(1+(rand(size(total(:,1)))-0.5)/5),total(:,1),20*round(total_color(:,1))+1,'filled');
    s_rct =scatter(2*ones(size(total(:,2))).*(1+(rand(size(total(:,2)))-0.5)/10),total(:,2),20*round(total_color(:,2))+1,'filled');
    s_hmc =scatter(3*ones(size(total(:,3))).*(1+(rand(size(total(:,3)))-0.5)/15),total(:,3),20*round(total_color(:,3))+1,'filled');
    hold off
    ax = gca
    ax.LineWidth = 4
    set(gca,'fontname','times');
    set(gca,'Fontsize',40);
    ax.FontSize = 40;
    ylabel('Participation Ratio','interpreter','latex','FontSize',40)
    set(h,{'linew'},{4})
    ax.LineWidth = 4
    s_rct.MarkerFaceColor = [0.5 0.5 0.5];
    s_rct.MarkerEdgeColor = [0.5 0.5 0.5];
    s_cyl.MarkerFaceColor = [1 0 0];
    s_cyl.MarkerEdgeColor = [1 0 0];
    s_hmc.MarkerFaceColor = [0 1 0];
    s_hmc.MarkerEdgeColor = [0 1 0];
    ylim([1 3])
    
    %% correlation between PR and angle shifts
    
    [R_cyl{iThreshold},P_cyl{iThreshold}] = corrcoef(cyl_data, cyl_colordata);
    [R_rct{iThreshold},P_rct{iThreshold}] = corrcoef(rct_data, rct_colordata);
    [R_hmc{iThreshold},P_hmc{iThreshold}] = corrcoef(hmc_data, hmc_colordata);
    
    corr_cyl(iThreshold) = R_cyl{iThreshold}(1,2);
    corr_rct(iThreshold) = R_rct{iThreshold}(1,2);
    corr_hmc(iThreshold) = R_hmc{iThreshold}(1,2);
    
    p_cyl(iThreshold) = P_cyl{iThreshold}(1,2);
    p_rct(iThreshold) = P_rct{iThreshold}(1,2);
    p_hmc(iThreshold) = P_hmc{iThreshold}(1,2);
    
    
    shifts_mean1(iThreshold) = mean(cyl_colordata);
    shifts_mean2(iThreshold) = mean(rct_colordata);
    shifts_mean3(iThreshold) = mean(hmc_colordata);
end
%%% correlation analysis

figure()
plot(angle_threshold,corr_cyl,'Color','red','LineWidth',4)
hold on
plot(angle_threshold,corr_rct,'color','[0.5 0.5 0.5]','LineWidth',4)
plot(angle_threshold,corr_hmc,'color', '[0 1 0]','LineWidth',4)
% plot(angle_threshold,p_cyl,'--','Color','red','LineWidth',2)
% plot(angle_threshold,p_rct,'--','color','[0.5 0.5 0.5]','LineWidth',2)
% plot(angle_threshold,p_hmc,'--','color', '[0 1 0]','LineWidth',2)
% plot(angle_threshold,0.05*ones(size(p_hmc)),'--','color', '[0 1 0]','LineWidth',4)
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',40);
ax.FontSize = 40;
xlabel('Angle Threshold [degrees]')
ylabel('Correlation coefficient')

%%% figure

derivative1_shifts= diff(shifts_mean1)./diff(angle_threshold);
derivative2_shifts= diff(shifts_mean2)./diff(angle_threshold);
derivative3_shifts= diff(shifts_mean3)./diff(angle_threshold);

total_derivative = mean([derivative1_shifts; derivative2_shifts; derivative3_shifts]);
total_shifts = mean([shifts_mean1; shifts_mean2; shifts_mean3]);


figure(30)
left_color = [0 0 0];
right_color = [0 0 0];
set(figure(30),'defaultAxesColorOrder',[left_color; right_color]);
hold on
yyaxis left
plot(angle_threshold, round(300./shifts_mean1),'-','Color','red','LineWidth',3)
plot(angle_threshold, round(300./shifts_mean2),'-','color','[0.5 0.5 0.5]','LineWidth',3)
plot(angle_threshold, round(300./shifts_mean3),'-','color','[0 1 0]','LineWidth',3)
xlabel('Angle Threshold [degrees]')
ylabel('Time between registrations [s]')
yyaxis right
plot(angle_threshold(1:end-1), derivative1_shifts,':','color','red','LineWidth',3)
plot(angle_threshold(1:end-1), derivative2_shifts,':','color','[0.5 0.5 0.5]','LineWidth',3)
plot(angle_threshold(1:end-1), derivative3_shifts,':','color','[0 1 0]','LineWidth',3)
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',24);
ax.FontSize = 24;
xlabel('Angle Threshold [degrees]')
ylabel('Derivative of the time delay')
xlim([angle_threshold(1) angle_threshold(end-1)])

figure(40)
left_color = [0 0 0];
right_color = [1 0 0];
set(figure(40),'defaultAxesColorOrder',[left_color; right_color]);
hold on
yyaxis left
plot(angle_threshold, round(300./total_shifts),'Color','black','LineWidth',3)
xlabel('Angle Threshold [degrees]')
ylabel('Time between registrations [s]')
yyaxis right
plot(angle_threshold(1:end-1), total_derivative,':','color','red','LineWidth',3)
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',24);
ax.FontSize = 24;
xlabel('Angle Threshold [degrees]')
ylabel('Derivative of the time delay')
xlim([angle_threshold(1) angle_threshold(end-1)])