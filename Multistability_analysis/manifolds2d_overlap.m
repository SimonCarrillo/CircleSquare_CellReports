clc
clear all

set(0,'DefaultFigureWindowStyle','docked')

%% experiment to analyze

time_bin = '1000'; % in ms
expt = 'CircleSquare'; % experiment setup
animal = 'M39'; % mouse number
loadStr = 'TempCorr'; % 'Zlinear' or 'TempCorr'

save_output = 0; % 1 for saving angles and proportion of overlap

discard_tech= ["kc" "neg" "randPos"];

for idiscard = 1:3
    
    discardPop = discard_tech(idiscard);
    
    ratio = '8';
    remove_outliers = 0; %1 for removing outliers
    alpharadius = Inf; % for alphashape
    ISO_output = 1; % 1 for plot
    
    path_ISO = strcat('ISO_','CellDiscarded_',discardPop,'Ratio',ratio,loadStr,time_bin,expt,'_',animal,'_','.mat');
    load(path_ISO);
    
    %% retrieve daily data
    
    days = [5];
    
    for iDay = days
        
        ISOmap{1,iDay} = ISO{iDay,2}; % extract dim reduced ISOMAP
        sessList{1,iDay} = ISO{iDay,4}; % extract session information
        sessStrList{1,iDay} = ISO{iDay,5}; % extract session information
        
    end
    
    for iDay = days
        
        ISO_data = ISOmap{1,iDay};
        sess = sessList{1,iDay};
        sess_Str = sessStrList{1,iDay};
        
        %% figures ISO
        
        if ISO_output == 1
            
            index1 = 0;
            index2 = 0;
            index3 = 0;
            
            figure()
            for Ii = 1:length(ISO_data(:,1))
                
                if sess_Str(Ii,1) == 'H' && (sess_Str(Ii,4) == '1' || sess_Str(Ii,4) == '2')
                    
                    index1 = index1 + 1;
                    
                    HMC_store{iDay}(index1,:,:,:) = ISO_data(Ii,:);

                    
                elseif sess_Str(Ii,1) == 'R' && (sess_Str(Ii,4) == '1' || sess_Str(Ii,4) == '2')
                    
                    index2 = index2 + 1;
                    
                    RCT_store{iDay}(index2,:,:,:) = ISO_data(Ii,:);

                    
                elseif sess_Str(Ii,1) == 'C' && (sess_Str(Ii,4) == '1' || sess_Str(Ii,4) == '2')
                    
                    index3 = index3 + 1;
                    
                    CYL_store{iDay}(index3,:,:,:) = ISO_data(Ii,:);

                    
                end
                
            end

        end
        
        %% least square fit and plot normal vector to best fit plane
        
        %%% for plotting
        
        title_expt = strcat('DiscardedCells',discardPop,'Ratio',ratio,animal,'Day',int2str(iDay));
        
        %% best fit plane
        
        % remove outliers or not
        if remove_outliers == 1
            
            HMC = rmoutliers(HMC_store{1,iDay},'quartiles');
            CYL = rmoutliers(CYL_store{1,iDay},'quartiles');
            RCT = rmoutliers(RCT_store{1,iDay},'quartiles');
            
        else
            
            HMC = HMC_store{1,iDay};
            CYL = CYL_store{1,iDay};
            RCT = RCT_store{1,iDay};
            
        end

 %% planar fit
 
        figure()
        [n,center] = fitNormal(HMC,1,'[0 1 0]')
        [n2,center2] = fitNormal(CYL,1,'[1 0 0]')
        [n3,center3] = fitNormal(RCT,1,'[0.5 0.5 0.5]')
        hold on
        
        x = linspace(min(HMC(:,1)),max(HMC(:,1)),100);
        y = linspace(min(HMC(:,2)),max(HMC(:,2)),100);
        [xg, yg] = meshgrid(x, y);
        
        x2 = linspace(min(CYL(:,1)),max(CYL(:,1)),100);
        y2 = linspace(min(CYL(:,2)),max(CYL(:,2)),100);
        [xg2, yg2] = meshgrid(x2, y2);
        
        x3 = linspace(min(RCT(:,1)),max(RCT(:,1)),100);
        y3 = linspace(min(RCT(:,2)),max(RCT(:,2)),100);
        [xg3, yg3] = meshgrid(x3, y3);
        
        z = -n(1)/n(3)*(xg)-n(2)/n(3)*(yg);
        z2 = -n2(1)/n2(3)*(xg2)-n2(2)/n2(3)*(yg2);
        z3 = -n3(1)/n3(3)*(xg3)-n3(2)/n3(3)*(yg3);
        
        
        h1 =surf(xg,yg,z,'FaceColor','[0 1 0]','FaceAlpha',0.5)
        h2 =surf(xg2,yg2,z2,'FaceColor','[1 0 0]','FaceAlpha',0.5)
        h3 = surf(xg3,yg3,z3,'FaceColor','[0.5 0.5 0.5]','FaceAlpha',0.5)
        legend([h1 h2 h3],'HomeCage','Cylinder','Rectangle')
        grid off
        box off
        view(3)
        
        title(title_expt)
        hold off
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'ztick',[])
        set(gcf,'color','w')
        axis off
        
        %% angle between normal vectors
        
        cells_kept = size(ISO{iDay,1},2)
        angle_rct_cyl(iDay,idiscard) = rad2deg(acos(dot(n2,n3)))
              
        %% 2D manifold fitting
        
        manifold{iDay} = alphaShape(HMC(:,1),HMC(:,2),HMC(:,3),alpharadius);
        manifold2{iDay} = alphaShape(RCT(:,1),RCT(:,2),RCT(:,3),alpharadius);
        manifold3{iDay} = alphaShape(CYL(:,1),CYL(:,2),CYL(:,3),alpharadius);
        
        
        %% intersection of 2D manifolds rct and cyl
        
        id1=inShape(manifold2{1,iDay},CYL_store{1,iDay}(:,1),CYL_store{1,iDay}(:,2),CYL_store{1,iDay}(:,3));
        id2=inShape(manifold3{1,iDay},RCT_store{1,iDay}(:,1),RCT_store{1,iDay}(:,2),RCT_store{1,iDay}(:,3));
        
        x1 = CYL_store{1,iDay}(:,1);
        y1 = CYL_store{1,iDay}(:,2);
        z1 = CYL_store{1,iDay}(:,3);
        
        x2 = RCT_store{1,iDay}(:,1);
        y2 = RCT_store{1,iDay}(:,2);
        z2 = RCT_store{1,iDay}(:,3);
        
        shp3_intersection{iDay}=alphaShape([x1(id1); x2(id2)], [y1(id1); y2(id2)], [z1(id1); z2(id2)],alpharadius);
      
        proportion_cylrct(iDay,idiscard) = sum(id1(:) == 1)/length(id1)
        proportion_rctcyl(iDay,idiscard) = sum(id2(:) == 1)/length(id2)
        
        clear id1 id2
        
        %% three environment and intersection plotted per day
        
        figure(10+iDay+idiscard)
        plot(manifold{1,iDay}, 'FaceAlpha',0.5, 'FaceColor','[0 1 0]')
        hold on
        plot(manifold2{1,iDay}, 'FaceAlpha',0.5, 'FaceColor','[1 0 0]')
        plot(manifold3{1,iDay}, 'FaceAlpha',0.5, 'FaceColor','[0.5 0.5 0.5]')
        view(3)
        box off
        grid off
        hold off
        title(title_expt)
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'ztick',[])
        set(gcf,'color','w')
        axis off
        
        
    end
    
    clear manifold manifold2 manifold3 shp3_intersection HMC_store CYL_store RCT_store
    
end

%% plot of all results

figure()
t1 =plot(proportion_cylrct(:,1),'-o','LineWidth',4,'Color','g')
hold on
t2 = plot(proportion_cylrct(:,2),'-o','LineWidth',4,'Color','r')
t3 =plot(proportion_cylrct(:,3),'-o','LineWidth',4,'Color','[0 0.5 0]')
%title(title_expt)
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

figure()
t1 =plot(proportion_rctcyl(:,1),'-o','LineWidth',4,'Color','g')
hold on
t2 =plot(proportion_rctcyl(:,2),'-o','LineWidth',4,'Color','r')
t3 =plot(proportion_rctcyl(:,3),'-o','LineWidth',4,'Color','[0 0.5 0]')
%title(title_expt)
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

figure()
t1 =plot(angle_rct_cyl(:,1),'-o','LineWidth',4,'Color','g')
hold on
t2 =plot(angle_rct_cyl(:,2),'-o','LineWidth',4,'Color','r')
t3 =plot(angle_rct_cyl(:,3),'-o','LineWidth',4,'Color','[0 0.5 0]')
%title(title_expt)
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

%% save files

if save_output == 1
    path_angle = strcat('angle','_',animal,'_ratio',ratio,'.mat');
    path_overlap = strcat('overlap','_',animal,'_ratio',ratio,'.mat');
    
    
    save(path_angle,'angle_rct_cyl')
    
    save(path_overlap,'proportion_cylrct','proportion_rctcyl')
    
end
