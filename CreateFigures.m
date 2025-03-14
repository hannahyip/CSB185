% add path for data
dir_data = './Data';
addpath(dir_data)

%Load experimental ISGF3 activity and mRNA data
load([dir_data,'/BetaISGF3.mat']);
load([dir_data,'/LambdaISGF3.mat']);
load([dir_data,'/BetaTime.mat']);
load([dir_data,'/BetaCommonGenes.mat']);
load([dir_data,'/LambdaCommonGenes.mat']);
load([dir_data,'/GeneTime.mat']);
load([dir_data,'/BetaBetaGenes.mat']);
load([dir_data,'/LambdaBetaGenes.mat']);
load([dir_data,'/Beta_Pri_C1.mat']);
load([dir_data,'/Beta_Pri_C2.mat']);
load([dir_data,'/Beta_Pri_Lambda_C1.mat']);
load([dir_data,'/Beta_Pri_Lambda_C2.mat']);

%% Parameters and Data to Change
BetaCommonGenes=Beta_Pri_C1;
LambdaCommonGenes=Beta_Pri_Lambda_C1;
optParams = [0.0015,0.1045112568,0.001557772306,5.893558568];
optRMSD = 0.5692592241;

BetaCommonGenes1=Beta_Pri_C2;
LambdaCommonGenes1=Beta_Pri_Lambda_C2;
optParams1 = [0.0015,9.930434494,0.6815296349,0.6199464567];
optRMSD1 = 0.6917609861;

% Redefines GeneTime
GeneTime = [0; 30; 60; 120; 240; 480; 720; 1440; 2160];

%% Normalize ISGF3 and RNA experimental data PART 1

totalISGF3=1; %max nuclear ISGF3
maxPercentage=1;
minPercentage=0.0025;

%scale EMSA proportional to amount of ISGF3 in nucleus (0.25% basal, %max)
allEMSA=[BetaISGF3,LambdaISGF3];
minEMSA=min(allEMSA,[],'all'); 
allEMSA=allEMSA-minEMSA;
maxEMSA=max(allEMSA,[],'all');
EMSAScaled=(allEMSA./maxEMSA)*(maxPercentage*totalISGF3);
EMSAScaled=EMSAScaled+(minPercentage*totalISGF3);

EMSABetaScaled=EMSAScaled(:,1);
EMSALambdaScaled=EMSAScaled(:,2);

%Normalize RNA    
allCommonGenes_data=[BetaCommonGenes,LambdaCommonGenes];
minGene_data=min(allCommonGenes_data,[],'all');
allCommonGenes_data=allCommonGenes_data-minGene_data;
maxGene_data=max(allCommonGenes_data,[],'all'); 

BetaRNA=allCommonGenes_data(:,1)./maxGene_data;
BetaRNA(BetaRNA<0)=0;

LambdaRNA=allCommonGenes_data(:,2)./maxGene_data;
LambdaRNA(LambdaRNA<0)=0;

%Calculating basal RNA concentration by taking average of Beta and IFN
% Lambda basal conditions
avgBasalRNA=mean([BetaRNA(1),LambdaRNA(1)]);

%Calculating basal ISGF3 concentration by taking average of IFN Beta and
% Lambda basal conditions
avgBasalISGF3=mean([EMSABetaScaled(1),EMSALambdaScaled(1)]);

%Interpolation of IFNBeta-induced ISGF3 data
BetamakimaFit=interp1(BetaTime,EMSABetaScaled,[0:800],'makima');

%Interpolation of IFNLambda-induced ISGF3 data
LambdamakimaFit=interp1(BetaTime,EMSALambdaScaled,[0:800],'makima');

x=avgBasalRNA;
time=[0:10:84000];
ISGF3=avgBasalISGF3;

%Run steady state model with optimized parameter set
[t_ss_new,y_ss_new]=ode15s(@(t,x) GeneSteadyState(t,x,ISGF3,optParams), ...
                                time,x);
t=[1:1:2500]; %time
SSInitial_new=y_ss_new(end); 

%Run ODE with best fit parameter set
[t_beta_new,y_beta_new]=ode15s(@(t,x) ISGF3GeneReg(t,x, ...
                                    BetamakimaFit,optParams),t, ...
                                    SSInitial_new);

[t_lambda_new,y_lambda_new]=ode15s(@(t,x) ISGF3GeneReg(t,x, ...
                                        LambdamakimaFit,optParams),t, ...
                                        SSInitial_new);

%Scale model proportional to amount of RNA
allRNA_sim=[y_beta_new,y_lambda_new];
minRNA=min(allRNA_sim,[],'all');
allRNA_sim=allRNA_sim-minRNA;
maxRNA_sim=max(allRNA_sim,[],'all');

%Calculating percent max of model RNA output
normmRNA_beta_new=allRNA_sim(:,1)./maxRNA_sim;

normmRNA_lambda_new=allRNA_sim(:,2)./maxRNA_sim;


%% Normalize ISGF3 and RNA experimental data PART 2

%Normalize RNA    
allCommonGenes_data=[BetaCommonGenes1,LambdaCommonGenes1];
minGene_data=min(allCommonGenes_data,[],'all');
allCommonGenes_data=allCommonGenes_data-minGene_data;
maxGene_data=max(allCommonGenes_data,[],'all'); 

BetaRNA1=allCommonGenes_data(:,1)./maxGene_data;
BetaRNA1(BetaRNA1<0)=0;

LambdaRNA1=allCommonGenes_data(:,2)./maxGene_data;
LambdaRNA1(LambdaRNA1<0)=0;

%Calculating basal RNA concentration by taking average of Beta and IFN
% Lambda basal conditions
avgBasalRNA=mean([BetaRNA1(1),LambdaRNA1(1)]);

%Run steady state model with optimized parameter set
[t_ss_new,y_ss_new]=ode15s(@(t,x) GeneSteadyState(t,x,ISGF3,optParams1), ...
                                time,x);
t=[1:1:2500]; %time
SSInitial_new=y_ss_new(end); 

%Run ODE with best fit parameter set
[t_beta_new1,y_beta_new1]=ode15s(@(t,x) ISGF3GeneReg(t,x, ...
                                    BetamakimaFit,optParams1),t, ...
                                    SSInitial_new);

[t_lambda_new1,y_lambda_new1]=ode15s(@(t,x) ISGF3GeneReg(t,x, ...
                                        LambdamakimaFit,optParams1),t, ...
                                        SSInitial_new);

%Scale model proportional to amount of RNA
allRNA_sim=[y_beta_new1,y_lambda_new1];
minRNA=min(allRNA_sim,[],'all');
allRNA_sim=allRNA_sim-minRNA;
maxRNA_sim=max(allRNA_sim,[],'all');

%Calculating percent max of model RNA output
normmRNA_beta_new1=allRNA_sim(:,1)./maxRNA_sim;

normmRNA_lambda_new1=allRNA_sim(:,2)./maxRNA_sim;

%% Plot Graphs
  figure
     subplot(2,2,1)
        plot(t_beta_new, y_beta_new, '-m', 'LineWidth', 2.5) 
        hold on
        plot(t_beta_new1, y_beta_new1, '-g', 'LineWidth', 2.5) 
        title('IFNBeta Model (Not Normalized)','FontSize',18, 'FontWeight','bold')
        xlabel('Time (minutes)', 'FontSize', 18, 'FontWeight', 'bold')
        ylabel('RNA (A.U.)', 'FontSize', 18, 'FontWeight', 'bold')
        ax = gca;
        autoY = get(gca, 'YLim');
        ax.YLim = autoY;
        legend({'Beta mRNA model, Cluster 1', 'Beta mRNA model, Cluster 2'}, 'FontSize', 12, 'FontWeight', 'bold', 'Location', 'best')  % Add legend
        hold off

    subplot(2,2,2)
        plot(t_lambda_new, y_lambda_new, '-r', 'LineWidth', 2.5) 
        hold on
        plot(t_lambda_new1, y_lambda_new1, '-b', 'LineWidth', 2.5) 
        title('IFNBeta Model (Not Normalized)','FontSize',18, 'FontWeight','bold')
        xlabel('Time (minutes)', 'FontSize', 18, 'FontWeight', 'bold')
        ylabel('RNA (A.U.)', 'FontSize', 18, 'FontWeight', 'bold')
        ax = gca;
        autoY = get(gca, 'YLim');
        ax.YLim = autoY;
        legend({'Lambda mRNA model, Cluster 1', 'Lambda mRNA model, Cluster 2'}, 'FontSize', 12, 'FontWeight', 'bold', 'Location', 'best')  % Add legend
        hold off

   subplot(2,2,3)
        plot(GeneTime,BetaRNA,'mx',t_beta_new,normmRNA_beta_new, ...
                ':m','LineWidth',2.5,'MarkerSize',12)
        hold on
        plot(GeneTime,BetaRNA1,'gx',t_beta_new,normmRNA_beta_new1, ...
                ':g','LineWidth',2.5,'MarkerSize',12)
        title('Beta RNA with Optimized Parameters (Normalized)','FontSize',18, ...
                'FontWeight','bold')
        xlabel('Time (minutes)','FontSize',18,'FontWeight','bold')
        ylabel('Percent Max mRNA Expression','FontSize',18,'FontWeight','bold')
        legend({'Beta mRNA data, Cluster 1','Beta mRNA model, Cluster 1', ...
                'Beta mRNA data, Cluster 2','Beta mRNA model, Cluster 2'}, ...
                'FontSize',12,'FontWeight','bold')
        ax=gca;
        ax.YLim=[0,1.2]; 
        hold off

    subplot(2,2,4)
       plot(GeneTime,LambdaRNA,'ro',t_lambda_new, normmRNA_lambda_new, ...
               ':r','LineWidth',2.5,'MarkerSize',12)
       hold on
       plot(GeneTime,LambdaRNA1,'bo',t_lambda_new, normmRNA_lambda_new1, ...
               ':b','LineWidth',2.5,'MarkerSize',12)
       title('Lambda RNA with Optimized Parameters (Normalized)','FontSize',18, ...
                'FontWeight','bold')
       xlabel('Time (minutes)','FontSize',18,'FontWeight','bold')
       ylabel('Percent Max mRNA Expression','FontSize',18,'FontWeight','bold')
       legend({'Lambda mRNA data, Cluster 1','Lambda mRNA model, Cluster 1', ...
                'Lambda mRNA data, Cluster 2','Lambda mRNA model, Cluster 2'}, ...
                'FontSize',12,'FontWeight','bold')
       ax=gca;
       ax.YLim=[0,1.2]; 
       hold off




