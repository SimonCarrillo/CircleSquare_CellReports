clc
clear all

set(0,'DefaultFigureWindowStyle','docked')

%% experiment to analyze

n_dimensions = ["3" "5" "8" "10" "25" "50" "100" "1"];

for idim = 1:length(n_dimensions)
    
    mouselist = ["M39"];
    ratiolist = ["2"];
    
    dimension = n_dimensions(idim);
    
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
            
            path_ISO = strcat('/Volumes/GoogleDrive/My Drive/Fenton_Lab/CellPress_CircleSquare/Response_Jan2022/Analysis_codes/Final_figures/Raw_data/Participation_ratio_ndimensions_parametric/','ISO_',dimension,'_','CellDiscarded_',discardPop,'Ratio',ratio,loadStr,time_bin,expt,'_',animal,'_','.mat');
            load(path_ISO);
            
            %% retrieve daily data
            
            days = linspace(1,size(ISO,1),size(ISO,1));
            
            for iDay = days
                
                ISOmap{1,iDay} = ISO{iDay,2}; % extract dim reduced ISOMAP
                spikes{1,iDay} = ISO{iDay,1}; % extract dim reduced ISOMAP
                sessList{1,iDay} = ISO{iDay,4}; % extract session information
                sessStrList{1,iDay} = ISO{iDay,5}; % extract session information
                
            end
            
            for iDay = days
                
                if dimension ~= "1"
                    ISO_data = ISOmap{1,iDay};
                else
                    ISO_data = spikes{1,iDay};
                    total_N(iDay) = size(ISO_data,2);
                end
                
                sess = sessList{1,iDay};
                sess_Str = sessStrList{1,iDay};
                
                %% figures ISO
                
                if ISO_output == 1
                    
                    index1 = 0;
                    index2 = 0;
                    index3 = 0;
                    
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
                
                HMC1 = HMC(1:300,:);
                HMC2 = HMC(301:end,:);
                CYL1 = CYL(1:300,:);
                CYL2 = CYL(301:end,:);
                RCT1 = RCT(1:300,:);
                RCT2 = RCT(301:end,:);
                
                %% Participation ratio CYL
                
                CNeuron = cov(CYL1);
                
                PRNumer = sum((diag(CNeuron))).^2;
                
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_cyl1(iMouse,iDay) = PRNumer/PRDenom
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                
                %% Participation ratio CYL
                
                CNeuron = cov(CYL2);
                
                PRNumer = sum((diag(CNeuron))).^2;
                
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_cyl2(iMouse,iDay) = PRNumer/PRDenom
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                %% Participation ratio RCT
                
                CNeuron = cov(RCT1);
                
                PRNumer = sum((diag(CNeuron))).^2;
                
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_rct1(iMouse,iDay)  = PRNumer/PRDenom
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                
                %% Participation ratio RCT
                
                CNeuron = cov(RCT2);
                
                PRNumer = sum((diag(CNeuron))).^2;
                
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_rct2(iMouse,iDay)  = PRNumer/PRDenom
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                
                %% Participation ratio HMC
                
                CNeuron = cov(HMC1);
                PRNumer = sum((diag(CNeuron))).^2;
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_hmc1(iMouse,iDay) = PRNumer/PRDenom
                
                clear PRDenom1 PRNumer CNeuron PRDenom
                
                %% Participation ratio HMC
                
                CNeuron = cov(HMC2);
                PRNumer = sum((diag(CNeuron))).^2;
                
                for i = 1:size(CNeuron,1)
                    PRDenom1(i,i) = (CNeuron(i,i)).^2;
                end
                
                PRDenom = sum(PRDenom1,'all');
                PR_hmc2(iMouse,iDay) = PRNumer/PRDenom
                
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
    HMC_PR(:,idim) = [HMC1_PR(ind_zeros); HMC2_PR(ind_zeros)];
    RCT_PR(:,idim) = [RCT1_PR(ind_zeros); RCT2_PR(ind_zeros)];
    CYL_PR(:,idim) = [CYL1_PR(ind_zeros); CYL2_PR(ind_zeros)];
    
    clear HMC1_PR HMC2_PR CYL1_PR CYL2_PR RCT1_PR RCT2_PR HMC_store CYL_store RCT_store
    
end

%% summary data

ndim_vector = [3 5 8 10 25 50 100 mean(nonzeros(total_N))];

figure()
t1 = stdshade(RCT_PR,0.2,[0.5 0.5 0.5],ndim_vector)
hold on
stdshade(CYL_PR,0.2,'red',ndim_vector)
stdshade(HMC_PR,0.2,'green',ndim_vector)
t2 = scatter(ndim_vector,ndim_vector,100,'MarkerFaceColor','black')
xlabel('Number of dimensions - IsoMap')
ylabel('Participation Ratio')
legend([t1 t2],'Experimental','Chance','location','northwest')
ax = gca
ax.LineWidth = 4
set(gca,'fontname','times');
set(gca,'Fontsize',40);
ax.FontSize = 40;
xlim([0 mean(nonzeros(total_N))+5])
ylim([0 mean(nonzeros(total_N))+5])

axes('Position',[.4 .7 .2 .2])
box on
stdshade(RCT_PR(:,1:4),0.2,[0.5 0.5 0.5],ndim_vector(1:4))
hold on
stdshade(CYL_PR(:,1:4),0.2,'red',ndim_vector(1:4))
stdshade(HMC_PR(:,1:4),0.2,'green',ndim_vector(1:4))
scatter(ndim_vector(1:4),ndim_vector(1:4),100,'MarkerFaceColor','black')
xlabel('')
ylabel('')
xticks([3 5 8 10])
yticks([3 5 8 10])
xlim([2 11])
ylim([2 11])
