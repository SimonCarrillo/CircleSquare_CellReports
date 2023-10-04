clc
clear all

set(0,'DefaultFigureWindowStyle','docked')

%% experiment to analyze

n_dimensions = ["3" "5" "8" "10" "25" "50" "100" "1"];
save_output = 1; % 1 for save

for idim = 1:length(n_dimensions)
    
    display(idim)
    
    mouselist = ["M29" "M34" "M39" "M35" "M20" "M19"];
    ratiolist = ["2" "8" "2" "8" "20" "20"];
    
    dimension = n_dimensions(idim);
    
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
            
            path_ISO = strcat('/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/Response_Jan2022/Analysis_codes/Final_figures/Raw_data/Participation_ratio_ndimensions_parametric/','ISO_',dimension,'_','CellDiscarded_',discardPop,'Ratio',ratio,loadStr,time_bin,expt,'_',animal,'_','.mat');
            load(path_ISO);
            
            %% retrieve daily data
            
            days{iMouse} = linspace(1,size(ISO,1),size(ISO,1));
            
            for iDay = days{iMouse}
                
                ISOmap{iMouse,iDay}= ISO{iDay,2}; % extract dim reduced ISOMAP
                spikes{iMouse,iDay} = ISO{iDay,1}; % extract dim reduced ISOMAP
                sessList{iMouse,iDay} = ISO{iDay,4}; % extract session information
                sessStrList{iMouse,iDay} = ISO{iDay,5}; % extract session information
                
            end
            
            for iDay = days{iMouse}
                
                if dimension ~= "1"
                    ISO_data = ISOmap{iMouse,iDay};
                else
                    ISO_data = spikes{iMouse,iDay};
                    total_N(iMouse,iDay) = size(ISO_data,2);
                end
                
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
                PR_cyl1{iMouse,iDay} = PRNumer/PRDenom;
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                
                %% Participation ratio CYL
                
                CNeuron = cov(CYL2{iMouse,iDay});
                
                PRNumer = sum((diag(CNeuron))).^2;
                
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_cyl2{iMouse,iDay} = PRNumer/PRDenom;
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                %% Participation ratio RCT
                
                CNeuron = cov(RCT1{iMouse,iDay});
                
                PRNumer = sum((diag(CNeuron))).^2;
                
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_rct1{iMouse,iDay}  = PRNumer/PRDenom;
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                
                %% Participation ratio RCT
                
                CNeuron = cov(RCT2{iMouse,iDay});
                
                PRNumer = sum((diag(CNeuron))).^2;
                
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_rct2{iMouse,iDay}  = PRNumer/PRDenom;
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                
                %% Participation ratio HMC
                
                CNeuron = cov(HMC1{iMouse,iDay});
                PRNumer = sum((diag(CNeuron))).^2;
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_hmc1{iMouse,iDay} = PRNumer/PRDenom;
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                
                %% Participation ratio HMC
                
                CNeuron = cov(HMC2{iMouse,iDay});
                PRNumer = sum((diag(CNeuron))).^2;
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_hmc2{iMouse,iDay}= PRNumer/PRDenom;
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                
                
                
            end
            
        end
    end
    
    
    for iMouse = 1:length(mouselist)
        
        for iDay = days{iMouse}
            
            if dimension ~= "1"
                
                HMC1_PR_tmp{iMouse,iDay}= PR_hmc1{iMouse,iDay}/str2num(dimension);
                RCT1_PR_tmp{iMouse,iDay} = PR_rct1{iMouse,iDay}/str2num(dimension);
                CYL1_PR_tmp{iMouse,iDay} = PR_cyl1{iMouse,iDay}/str2num(dimension);
                HMC2_PR_tmp{iMouse,iDay} = PR_hmc2{iMouse,iDay}/str2num(dimension);
                RCT2_PR_tmp{iMouse,iDay}= PR_rct2{iMouse,iDay}/str2num(dimension);
                CYL2_PR_tmp{iMouse,iDay} = PR_cyl2{iMouse,iDay}/str2num(dimension);
                
                
            else
                
                HMC1_PR_tmp{iMouse,iDay} = PR_hmc1{iMouse,iDay}/total_N(iMouse,iDay);
                RCT1_PR_tmp{iMouse,iDay}= PR_rct1{iMouse,iDay}/total_N(iMouse,iDay);
                CYL1_PR_tmp{iMouse,iDay} = PR_cyl1{iMouse,iDay}/total_N(iMouse,iDay);
                HMC2_PR_tmp{iMouse,iDay}= PR_hmc2{iMouse,iDay}/total_N(iMouse,iDay);
                RCT2_PR_tmp{iMouse,iDay} = PR_rct2{iMouse,iDay}/total_N(iMouse,iDay);
                CYL2_PR_tmp{iMouse,iDay}= PR_cyl2{iMouse,iDay}/total_N(iMouse,iDay);
                
            end
            
           clear sess_Str 
        end
    end
    
    HMC1_PR = cell2mat(cellflat(HMC1_PR_tmp))';
    HMC2_PR = cell2mat(cellflat(HMC2_PR_tmp))';
    CYL1_PR = cell2mat(cellflat(CYL1_PR_tmp))';
    CYL2_PR = cell2mat(cellflat(CYL2_PR_tmp))';
    RCT1_PR = cell2mat(cellflat(RCT1_PR_tmp))';
    RCT2_PR = cell2mat(cellflat(RCT2_PR_tmp))';
    
    ind_zeros = find(HMC1_PR>0);
    HMC_PR(:,idim) = [HMC1_PR(ind_zeros); HMC2_PR(ind_zeros)];
    RCT_PR(:,idim) = [RCT1_PR(ind_zeros); RCT2_PR(ind_zeros)];
    CYL_PR(:,idim) = [CYL1_PR(ind_zeros); CYL2_PR(ind_zeros)];
    
    clear HMC1_PR HMC2_PR CYL1_PR CYL2_PR RCT1_PR RCT2_PR HMC_store CYL_store RCT_store PR_hmc1 PR_hmc2 PR_cyl1 PR_cyl2 PR_rct1 PR_rct2
    
end

%% summary data

ndim_vector = [3 5 8 10 25 50 100 370];

figure()
t1 = stdshade(RCT_PR,0.2,[0.5 0.5 0.5],ndim_vector)
hold on
stdshade(CYL_PR,0.2,'red',ndim_vector)
stdshade(HMC_PR,0.2,'green',ndim_vector)
t2 = plot(ndim_vector,ones(8,1),'--','color','black')
xlabel('Number of dimensions - IsoMap')
ylabel('Participation Ratio')
legend([t1 t2],'Experimental','Chance','location','southwest')
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',40);
ax.FontSize = 40;
xlim([0 max(ndim_vector)])
ylim([0 1.1])

axes('Position',[.6 .6 .2 .2])
box on
stdshade(RCT_PR(:,1:4),0.2,[0.5 0.5 0.5],ndim_vector(1:4))
hold on
stdshade(CYL_PR(:,1:4),0.2,'red',ndim_vector(1:4))
stdshade(HMC_PR(:,1:4),0.2,'green',ndim_vector(1:4))
xlabel('')
ylabel('')
xticks([3 5 8 10])
% yticks([3 5 8 10])
% xlim([2 11])
ylim([0.6 0.9])

%% save data

if save_output == 1
    
    path_PR = strcat('PR_all','300','.mat');
    save(path_PR,'RCT_PR','HMC_PR','CYL_PR','n_dimensions');
    
    
end
