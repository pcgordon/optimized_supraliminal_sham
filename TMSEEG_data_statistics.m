%% TEP STATISTICS
%Load data

n_sample=; %insert sample size

for i=1:n_sample

load() %Load data: EEG_interp (see TMSEEG_data_processing)   

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 45;
cfg.keeptrials = 'no';
cfg.trials = find(EEG_interp.trialinfo(:,1)==1) ;
TEPs_ACTIVE{i} = ft_timelockanalysis(cfg, EEG_interp);
cfg.trials = find(EEG_interp.trialinfo(:,1)==0); 
TEPs_SHAM{i} = ft_timelockanalysis(cfg, EEG_interp);
end

 % GRANDAVERAGE
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
GA_TEPs_SHAM    =  ft_timelockgrandaverage(cfg,TEPs_SHAM{:});  
GA_TEPs_ACTIVE    =  ft_timelockgrandaverage(cfg,TEPs_ACTIVE{:});
cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'individual';
GA_subt=ft_math(cfg, GA_TEPs_ACTIVE, GA_TEPs_SHAM);


% Prepare neighbours (for plotting and statistics)
EEGchannels = EEG_interp.label;
elec       = ft_convert_units(ft_read_sens('C:\Users\BNP Lab User\Desktop\fieldtrip-20160508\template\electrode\standard_1005.elc'),'mm');
index = find(any(cell2mat(cellfun(@strcmpi,EEGchannels',num2cell(repmat(elec.label,1,numel(EEGchannels)),1),'uniformoutput',false)),2));
elec.chanpos    = elec.chanpos(index,:);
elec.elecpos    = elec.elecpos(index,:);
elec.label      = elec.label(index,:);
neighbour_layout.elec       = elec;
cfg  = [];
neighbour_layout.layout  = ft_prepare_layout(cfg,neighbour_layout);
cfg             = [];
cfg.method      = 'triangulation';
cfg.elec        = elec;
cfg.channels    = 'EEG';
neighbour_layout.neighbor   = ft_prepare_neighbours(cfg);

%% Butterfly Plot

x_lim=[-0.04 0.3];
y_lim=[-35 7.5 ];

figure
hold on
plot(GA_subt.time, squeeze(mean(GA_TEPs_ACTIVE.individual(:, :,:),1)),'Linewidth',0.5,'Color',[144,144,144]./255)
plot(GA_subt.time, squeeze(mean(GA_TEPs_ACTIVE.individual(:,30,:),1)),'Linewidth',1,'Color',[53,94,182]./255)

plot(GA_subt.time, -15+squeeze(mean(GA_TEPs_SHAM.individual(:, :,:),1)),'Linewidth',0.5,'Color',[144,144,144]./255)
plot(GA_subt.time, -15+squeeze(mean(GA_TEPs_SHAM.individual(:,30,:),1)),'Linewidth',1,'Color',[53,94,182]./255)

plot(GA_subt.time, -30+squeeze(mean(GA_subt.individual(:, :,:),1)),'Linewidth',0.5,'Color',[144,144,144]./255)
plot(GA_subt.time, -30+squeeze(mean(GA_subt.individual(:,30,:),1)),'Linewidth',1,'Color',[53,94,182]./255)

plot([0.020 0.020],[y_lim(1) y_lim(2)], 'k:')
plot([0.040 0.04],[y_lim(1) y_lim(2)], 'k:')
plot([0.06 0.06],[y_lim(1) y_lim(2)], 'k:')
plot([0.09 0.09],[y_lim(1) y_lim(2)], 'k:')
plot([0.13 0.13],[y_lim(1) y_lim(2)], 'k:')
plot([0.25 0.25],[y_lim(1) y_lim(2)], 'k:')

xlim(x_lim)
plot([0 0],ylim,'-k','Linewidth',2)
text(0.25, 4,'REAL', 'FontWeight', 'Bold')
text(0.25, 4-15,'SHAM', 'FontWeight', 'Bold')
text(0.18, 4-32,['REAL ' char(8211) ' SHAM'], 'FontWeight', 'Bold')

ylim(y_lim)
yticks([-35:2.5:5])
yticklabels({  '-5' '-2.5' '0' '2.5' '5' '' '-5' '-2.5' '0' '2.5' '5' '' '-5' '-2.5' '0' '2.5' '5' ''})
set(gca,'fontsize',14, 'FontWeight', 'Bold')
set(gcf,'units','points','position',[100,100,400,600])
%title('Audit Only')

%% T test - Identify the time course of the clusters

cfg = [];
cfg.latency = [0 1.0];
cfg.avgovertime = 'no';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.025;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbour_layout.neighbor; 
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization =1000;
design = zeros(2,2*n_sample);
for i = 1:n_sample
  design(1,i) = i;
end
for i = 1:n_sample
  design(1,n_sample+i) = i;
end
design(2,1:n_sample)        = 1;
design(2,n_sample+1:2*n_sample) = 2;
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;
stat = ft_timelockstatistics(cfg, GA_TEPs_ACTIVE, GA_TEPs_SHAM);


for i=1:length(stat.posclusters)
    dummy=find(sum(stat.posclusterslabelmat==i)~=0)/1000;
    posclust_latency(i,:)=[dummy(1) dummy(end)];
end
for i=1:length(stat.negclusters)
    dummy=find(sum(stat.negclusterslabelmat==i)~=0)/1000;
    negclust_latency(i,:)=[dummy(1) dummy(end)];
end

 
 %% T test - fix time course of clusters
cfg = [];
cfg.latency = []; %Insert latency of significant clusters ( posclust_latency / negclust_latency)
cfg.avgovertime = 'yes';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.025;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbour_layout.neighbor; 
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization =1000;
design = zeros(2,2*n_sample);
for i = 1:n_sample
  design(1,i) = i;
end
for i = 1:n_sample
  design(1,n_sample+i) = i;
end
design(2,1:n_sample)        = 1;
design(2,n_sample+1:2*n_sample) = 2;
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;
stat = ft_timelockstatistics(cfg, GA_TEPs_ACTIVE, GA_TEPs_SHAM);

% Plot STAT
pval_lim=0.05;
cfg = [];
cfg.style = 'straight';
cfg.alpha  = pval_lim;
cfg.parameter = 'stat';
cfg.zlim   = [-4 4];
cfg.highlightsymbolseries  =  ['.','.','.','.','.'] ;
cfg.highlightsizeseries =[2 2 2 ].*10;
cfg.highlightcolorneg=[	221,125,30]./255;
cfg.highlightcolorpos=[	221,125,30]./255;
cfg.subplotsize = [1 2];
cfg.layout = neighbour_layout.layout;
ft_clusterplot(cfg, stat);
 colormap(mymap)
 
 
 
%% Plot time course of ROIs

chan=ismember(data_eeg.label, {'C1' 'C3' 'FC1' 'CF3'});
chan=ismember(data_eeg.label, {'F2' 'F4' 'FC2' 'FC4'});
chan=ismember(data_eeg.label, {'P5' 'P3' 'CP3' 'CP5'});

x_lim=[-0.04 0.30];
y_lim=[-4.5 4.5];

hold on
stdshade(squeeze(mean(GA_TEPs_ACTIVE.individual(:,chan,:),2)) ,0.3,[183,0,255]./255,GA_TEPs_ACTIVE.time,[],'sem',3)
stdshade(squeeze(mean(GA_TEPs_SHAM.individual(:,chan,:),2)) ,0.3,[7,109,3]./255,GA_TEPs_ACTIVE.time,[],'sem',3)
lineA=plot([0 0], [0 0], 'Color', [183,0,255]./255);
lineB=plot([0 0], [0 0], 'Color', [7,109,3]./255);
xlim(x_lim); ylim(y_lim); 
plot([0 0],ylim,'-k','Linewidth',4); 
plot(xlim,[0 0],'--k','Linewidth',2)
%xlabel('s'); ylabel(['\mu' 'V'])

patch([0.025 0.035 0.035 0.025], [y_lim(1) y_lim(1) y_lim(2) y_lim(2)],[0.2 0.2 0.2],'LineStyle', 'none','FaceAlpha',0.2)
patch([0.100 0.15 0.15 0.100], [y_lim(1) y_lim(1) y_lim(2) y_lim(2)],[0.2 0.2 0.2],'LineStyle', 'none','FaceAlpha',0.2)
patch([0.480 0.520 0.520 0.480], [y_lim(1) y_lim(1) y_lim(2) y_lim(2)],[0.2 0.2 0.2],'LineStyle', 'none','FaceAlpha',0.2)

patch([0.068 0.090 0.090 0.068], [y_lim(1) y_lim(1) y_lim(2) y_lim(2)],[0.2 0.2 0.2],'LineStyle', 'none','FaceAlpha',0.2)
patch([0.040 0.052 0.052 0.040], [y_lim(1) y_lim(1) y_lim(2) y_lim(2)],[0.2 0.2 0.2],'LineStyle', 'none','FaceAlpha',0.2)

legend( [lineA  lineB], {'ACTIVE' 'SHAM'}, 'Location','southeast')
set(gca,'fontsize',14, 'FontWeight', 'Bold')
set(gcf,'position',[200,200,600,200])



%% TRF STATISTICS

n_sample=; %insert sample size

for i=1:n_sample
load() %Load data: INDUCED_active (see TMSEEG_data_processing)   
load() %Load data: INDUCED_sham (see TMSEEG_data_processing)  
TOTAL_INDUCED_ACTIVE{i}=INDUCED_active;
TOTAL_INDUCED_SHAM{i}=INDUCED_sham;
end

%GRANDAVERAGE
cfg = [];
cfg.channel   = 'all';
cfg.keepindividual   = 'yes';
cfg.parameter = 'powspctrm';
GA_TFR_ACTIVE = ft_freqgrandaverage(cfg,TOTAL_INDUCED_ACTIVE{:});  
GA_TFR_SHAM = ft_freqgrandaverage(cfg,TOTAL_INDUCED_SHAM{:});  

cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'powspctrm';
GA_TFR_subt=ft_math(cfg, GA_TFR_ACTIVE, GA_TFR_SHAM);

%% PLOTS

z_lim=[-0.25 0.25];
xLIM=[-0.2 0.5];
yLIM=[3 45];
yTICKS=[3 7 13 21 30 45];
interp_scaling=2;

 subplot(3,2,1)
[T, F, DATA, SCALE]=newTFRplot(GA_TFR_SHAM.time, GA_TFR_SHAM.freq, squeeze(mean(mean(GA_TFR_SHAM.powspctrm(:,:,:,:),2),1)), interp_scaling);
plot_title=[ ];
newTFRplot_plot(T, F, DATA, SCALE, z_lim,{[-0.01 0.01]},1.5,16,yLIM, yTICKS,xLIM,...
    'jet', plot_title);
 subplot(3,2,2)
[T, F, DATA, SCALE]=newTFRplot(GA_TFR_SHAM.time, GA_TFR_SHAM.freq, squeeze(mean(mean(GA_TFR_SHAM.powspctrm(:,30,:,:),2),1)), interp_scaling);
plot_title=[];
newTFRplot_plot(T, F, DATA, SCALE, z_lim,{[-0.01 0.01]},1.5,16,yLIM, yTICKS,xLIM,...
    'jet', plot_title);

 subplot(3,2,3)
[T, F, DATA, SCALE]=newTFRplot(GA_TFR_ACTIVE.time, GA_TFR_ACTIVE.freq, squeeze(mean(mean(GA_TFR_ACTIVE.powspctrm(:,:,:,:),2),1)), interp_scaling);
plot_title=[];
newTFRplot_plot(T, F, DATA, SCALE, z_lim,{[-0.01 0.01]},1.5,16,yLIM, yTICKS,xLIM,...
    'jet', plot_title);
 subplot(3,2,4)
[T, F, DATA, SCALE]=newTFRplot(GA_TFR_ACTIVE.time, GA_TFR_ACTIVE.freq, squeeze(mean(mean(GA_TFR_ACTIVE.powspctrm(:,30,:,:),2),1)), interp_scaling);
plot_title=[ ];
newTFRplot_plot(T, F, DATA, SCALE, z_lim,{[-0.01 0.01]},1.5,16,yLIM, yTICKS,xLIM,...
    'jet', plot_title);

 subplot(3,2,5)
[T, F, DATA, SCALE]=newTFRplot(GA_TFR_subt.time, GA_TFR_subt.freq, squeeze(mean(mean(GA_TFR_subt.powspctrm(:,:,:,:),2),1)), interp_scaling);
plot_title=[ ];
newTFRplot_plot(T, F, DATA, SCALE, z_lim,{[-0.01 0.01]},1.5,16,yLIM, yTICKS,xLIM,...
    'jet', plot_title);
 subplot(3,2,6)
[T, F, DATA, SCALE]=newTFRplot(GA_TFR_subt.time, GA_TFR_subt.freq, squeeze(mean(mean(GA_TFR_subt.powspctrm(:,30,:,:),2),1)), interp_scaling);
plot_title=[ ];
newTFRplot_plot(T, F, DATA, SCALE, z_lim,{[-0.01 0.01]},1.5,16,yLIM, yTICKS,xLIM,...
    'jet', plot_title);

set(gcf,'position',[200,200,1300,1000])
 
 %% T test - fix time course of clusters
 
 %select a frequency of interest
FR=[4 7];
FR=[8 12];
FR=[13 20];
FR=[21 29];
FR=[30 45];

cfg = [];
cfg.latency = [0 1];
cfg.frequency=FR;
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'yes';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.01;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbour_layout.neighbor; 
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.01;
cfg.numrandomization = 1000;
cfg.computestat    = 'yes';
cfg.computecritval = 'yes';
cfg.computeprob    = 'yes';

design = zeros(2,2*n_sample);
for i = 1:n_sample
  design(1,i) = i;
end
for i = 1:n_sample
  design(1,n_sample+i) = i;
end
design(2,1:n_sample)        = 1;
design(2,n_sample+1:2*n_sample) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

stat = ft_freqstatistics(cfg, GA_TFR_ACTIVE,GA_TFR_SHAM);
 
 for i=1:length(stat.posclusters)
    dummy=find(sum(stat.posclusterslabelmat==i)~=0)/100;
    posclust_latency(i,:)=[dummy(1) dummy(end)];
end
for i=1:length(stat.negclusters)
    dummy=find(sum(stat.negclusterslabelmat==i)~=0)/100;
    negclust_latency(i,:)=[dummy(1) dummy(end)];
end
 

%% T test - fix time course of clusters

cfg = [];
cfg.latency = []; %Insert latency of significant clusters ( posclust_latency / negclust_latency)
cfg.frequency=[]; %select a frequency of interest
cfg.avgovertime = 'yes';
cfg.avgoverfreq = 'yes';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.01;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbour_layout.neighbor; 
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.01;
cfg.numrandomization = 1000;
cfg.computestat    = 'yes';
cfg.computecritval = 'yes';
cfg.computeprob    = 'yes';

design = zeros(2,2*n_sample);
for i = 1:n_sample
  design(1,i) = i;
end
for i = 1:n_sample
  design(1,n_sample+i) = i;
end
design(2,1:n_sample)        = 1;
design(2,n_sample+1:2*n_sample) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

stat = ft_freqstatistics(cfg, GA_TFR_ACTIVE,GA_TFR_SHAM);

% Plot STAT
pval_lim=0.010;
cfg = [];
cfg.style = 'straight';
cfg.alpha  = pval_lim;
cfg.parameter = 'stat';
cfg.zlim   = [-4 4];
cfg.highlightsymbolseries  =  ['.','.','.','.','.'] ;
cfg.highlightsizeseries =[2 2 2 ].*10;
cfg.highlightcolorneg=[	161,75,10]./255;
cfg.highlightcolorpos=[	161,75,10]./255;
cfg.subplotsize = [1 2];
cfg.layout = neighbour_layout.layout;
ft_clusterplot(cfg, stat);
 colormap(mymap)
