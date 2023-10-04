clc
clear all

mouselist = ["M29" "M34" "M39" "M35" "M20" "M19"];
ratiolist = ["4" "4" "4" "4" "4" "4"]; % "4" "4" "4" "4" "4" "4"

for iMouse = 1:length(mouselist)

    title_expt = strcat(mouselist(iMouse),'ratio',ratiolist(iMouse));
    
    animal = mouselist(iMouse); % mouse number
    ratio = ratiolist(iMouse); %discarding ratio
    
    path_angle = strcat('/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/Response_Jan2022/Analysis_codes/',animal,'/','angle','_',animal,'_ratio',ratio,'.mat');
    path_overlap = strcat('/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/Response_Jan2022/Analysis_codes/',animal,'/','overlap','_',animal,'_ratio',ratio,'.mat');
    
    angle_tmp = load(path_angle);
    proportion_tmp = load(path_overlap);
    
    
    if animal == "M20"
        
        angle_rct_cyl(1,:) = [nan nan nan];
        angle_rct_cyl(2:3,:) = angle_tmp.angle_rct_cyl(1:2,:);
        angle_rct_cyl(4,:) = [nan nan nan];
        angle_rct_cyl(5:9,:) = angle_tmp.angle_rct_cyl(3:7,:);
        
        proportion_rctcyl(1,:) = [nan nan nan];
        proportion_rctcyl(2:3,:) = proportion_tmp.proportion_rctcyl(1:2,:);
        proportion_rctcyl(4,:) = [nan nan nan];
        proportion_rctcyl(5:9,:) = proportion_tmp.proportion_rctcyl(3:7,:);
        
        proportion_cylrct(1,:) = [nan nan nan];
        proportion_cylrct(2:3,:) = proportion_tmp.proportion_cylrct(1:2,:);
        proportion_cylrct(4,:) = [nan nan nan];
        proportion_cylrct(5:9,:) = proportion_tmp.proportion_cylrct(3:7,:);        
 
    elseif animal == "M19"
        
        angle_rct_cyl(1:2,:) = angle_tmp.angle_rct_cyl(1:2,:);
        angle_rct_cyl(3,:) = [nan nan nan];
        angle_rct_cyl(4,:) = angle_tmp.angle_rct_cyl(3,:);        
        angle_rct_cyl(5,:) = [nan nan nan];
        angle_rct_cyl(6:9,:) = angle_tmp.angle_rct_cyl(4:7,:); 
        
        proportion_rctcyl(1:2,:) = proportion_tmp.proportion_rctcyl(1:2,:);
        proportion_rctcyl(3,:) = [nan nan nan];
        proportion_rctcyl(4,:) = proportion_tmp.proportion_rctcyl(3,:);        
        proportion_rctcyl(5,:) = [nan nan nan];
        proportion_rctcyl(6:9,:) = proportion_tmp.proportion_rctcyl(4:7,:); 
        
        proportion_cylrct(1:2,:) = proportion_tmp.proportion_cylrct(1:2,:);
        proportion_cylrct(3,:) = [nan nan nan];
        proportion_cylrct(4,:) = proportion_tmp.proportion_cylrct(3,:);        
        proportion_cylrct(5,:) = [nan nan nan];
        proportion_cylrct(6:9,:) = proportion_tmp.proportion_cylrct(4:7,:);  
        
    else
        
        proportion_rctcyl = proportion_tmp.proportion_rctcyl;
        proportion_cylrct = proportion_tmp.proportion_cylrct;  
        angle_rct_cyl = angle_tmp.angle_rct_cyl;

    end
    %% plot of all results
    
    figure()
    t1 =plot(proportion_cylrct(:,1),'-o','LineWidth',2.5,'Color','g')
    hold on
    t2 = plot(proportion_cylrct(:,2),'-o','LineWidth',2.5,'Color','r')
    t3 =plot(proportion_cylrct(:,3),'-o','LineWidth',2.5,'Color','[0 0.5 0]')
    grid on
    legend([t1 t2 t3],'kc','neg','randPos')
    ylabel('proportion','interpreter','latex','FontSize',40)
    xlabel('day','interpreter','latex','FontSize',40)
    t1 =plot(proportion_rctcyl(:,1),'--','LineWidth',4,'Color','g')
    hold on
    t2 =plot(proportion_rctcyl(:,2),'--','LineWidth',4,'Color','r')
    t3 =plot(proportion_rctcyl(:,3),'--','LineWidth',4,'Color','[0 0.5 0]')
    title(title_expt)
    grid on
    hold off
    ax = gca
    ax.LineWidth = 4
    set(gca,'fontname','times');
    set(gca,'Fontsize',40);
    ax.FontSize = 40;
    legend([t1 t2 t3],'kc','neg','randPos')
    ylabel('proportion','interpreter','latex','FontSize',40)
    xlabel('day','interpreter','latex','FontSize',40)
    
    mean_proportion(:,1) = (proportion_cylrct(:,1)+proportion_rctcyl(:,1))/2;
    mean_proportion(:,2) = (proportion_cylrct(:,2)+proportion_rctcyl(:,2))/2;
    mean_proportion(:,3) = (proportion_cylrct(:,3)+proportion_rctcyl(:,3))/2;
    
    angle_kc(iMouse,:) = angle_rct_cyl(:,1);
    angle_neg(iMouse,:) = angle_rct_cyl(:,2);
    angle_randPos(iMouse,:) = angle_rct_cyl(:,3);
    overlap_kc(iMouse,:) = mean_proportion(:,1);
    overlap_neg(iMouse,:) = mean_proportion(:,2);
    overlap_randPos(iMouse,:) = mean_proportion(:,3);
    
    
    figure()
    t1 =plot(mean_proportion(:,1),'-o','LineWidth',2.5,'Color','g')
    hold on
    t2 = plot(mean_proportion(:,2),'-o','LineWidth',2.5,'Color','r')
    t3 =plot(mean_proportion(:,3),'-o','LineWidth',2.5,'Color','[0 0.5 0]')
    grid on
    hold off
    ax = gca
    ax.LineWidth = 4
    set(gca,'fontname','times');
    set(gca,'Fontsize',40);
    ax.FontSize = 40;
    legend([t1 t2 t3],'kc','neg','randPos')
    ylabel('proportion','interpreter','latex','FontSize',40)
    xlabel('day','interpreter','latex','FontSize',40)
    title(title_expt)
    
    
    
    figure()
    t1 =plot(angle_rct_cyl(:,1),'-o','LineWidth',4,'Color','g')
    hold on
    t2 =plot(angle_rct_cyl(:,2),'-o','LineWidth',4,'Color','r')
    t3 =plot(angle_rct_cyl(:,3),'-o','LineWidth',4,'Color','[0 0.5 0]')
    title(title_expt)
    grid on
    hold off
    ax = gca
    ax.LineWidth = 4
    set(gca,'fontname','times');
    set(gca,'Fontsize',40);
    ax.FontSize = 40;
    legend([t1 t2 t3],'kc','neg','randPos')
    ylabel('angle','interpreter','latex','FontSize',40)
    xlabel('day','interpreter','latex','FontSize',40)

    clear mean_proportion
end

%% average and standard error representation with shade

figure()
t1 = stdshade(overlap_kc,0.2,'g',[1 2 3 8 9 10 15 16 17],10)
hold on
t2 = stdshade(overlap_neg,0.2,'r',[1 2 3 8 9 10 15 16 17],10)
t3 = stdshade(overlap_randPos,0.2,'black',[1 2 3 8 9 10 15 16 17],10)
hold off
legend([t1 t2 t3],'All cells','Anti-correlated cells removed','Random cells removed', 'location','best')

%% normal data check and equal variance

mean_kc = mean(overlap_kc,1,'omitnan');
mean_neg = mean(overlap_neg,1,'omitnan');
mean_randPos = mean(overlap_randPos,1,'omitnan');

raw_data = [mean_kc' mean_neg' mean_randPos']';
variancetest = vartestn(raw_data);
[normaltest1,p1] = kstest(raw_data(1,:));
[normaltest2,p2] = kstest(raw_data(2,:));
[normaltest3,p3] = kstest(raw_data(3,:));
[normaltest,p] = kstest(overlap_neg)

%% Kruskal Wallis

[pcells,~,stats_kruskalcells] = kruskalwallis(raw_data')

c_cells_nonparametric = multcompare(stats_kruskalcells,'ctype','dunn-sidak')

%% mean of every week per mouse per type - week1 [mean(mouse1) mean(mouse2) mean(mouse3)]
%% kruskal wallis for weeks

for iMouse = 1:6
    Inew = 1
    for Ii= 1:3     
        
        display(Inew:Inew+2)
        
        week_data_kc(iMouse,Ii) = mean(overlap_kc(iMouse,Inew:Inew+2),'omitnan')
        week_data_neg(iMouse,Ii) = mean(overlap_neg(iMouse,Inew:Inew+2),'omitnan')
        week_data_randPos(iMouse,Ii) = mean(overlap_randPos(iMouse,Inew:Inew+2),'omitnan')

        Inew = Inew+3;
    end
end

[pweeks,~,stats_kruskalweeks] = kruskalwallis(week_data_kc);
c_weeks_nonparametric = multcompare(stats_kruskalweeks,'ctype','dunn-sidak');

[pweeks2,~,stats_kruskalweeks2] = kruskalwallis(week_data_neg);
c_weeks_nonparametric2 = multcompare(stats_kruskalweeks2,'ctype','dunn-sidak');

[pweeks3,~,stats_kruskalweeks3] = kruskalwallis(week_data_randPos);
c_weeks_nonparametric3 = multcompare(stats_kruskalweeks3,'ctype','dunn-sidak');

%% two way anova for weeks

[pweeksanova,~,stats_weeks] = anova1(week_data_kc);
c_weeks_anova = multcompare(stats_weeks,'ctype','tukey-kramer');

[pweeks2anova,~,stats_weeks2] = anova1(week_data_neg);
c_weeks_anova2 = multcompare(stats_weeks2,'ctype','tukey-kramer');

[pweeks3anova,~,stats_weeks3] = anova1(week_data_randPos);
c_weeks_anova3 = multcompare(stats_weeks3,'ctype','tukey-kramer');

%% Dunnet's test

p = dunnett(stats_weeks)

%% boxplot cell type

figure()
h =boxplot(raw_data','Notch','off','Labels',{'all cells','anti-cofiring cells removed','random cells removed'},'Whisker',1,'MedianStyle','line')
ax = gca
box on
set(h,{'linew'},{4})
ax.LineWidth = 4
set(gca,'fontname','times');
ylabel('overlap','interpreter','latex','FontSize',40)
set(gca,'Fontsize',40);

%% boxplot week effect for kc

figure()
h =boxplot(week_data_kc','Notch','off','Labels',{'Week 1','Week 2','Week 3'},'Whisker',1,'MedianStyle','line')
hold on 
h2 =boxplot(week_data_neg','Notch','off','Labels',{'Week 1','Week 2','Week 3'},'Whisker',1,'MedianStyle','line')
ax = gca
box on
set(h,{'linew'},{4})
ax.LineWidth = 4
set(gca,'fontname','times');
ylabel('overlap','interpreter','latex','FontSize',40)
set(gca,'Fontsize',40);

%% boxplot week effect for kc

figure()
data = {week_data_kc', week_data_neg', week_data_randPos'}; 
h = boxplotGroup(data, 'PrimaryLabels', {'All' 'Anti' 'Rand'}, ...
  'SecondaryLabels',{'Week 1', 'Week 2','Week 3'}, 'InterGroupSpace', 2)
ax = gca
box on
ax.LineWidth = 4
set(gca,'fontname','times');
ylabel('overlap','interpreter','latex','FontSize',20)
set(gca,'Fontsize',20);