%% plot the simulation results of non-stationary noise removal 
% using Mosaic-T and Mosaic-TF
% complementary code for the publication 
% "Non-stationary Noise Removal from Repeated Sweep Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to JASA Express Letters
% on 30.04.2024
%% plots figure 1 from the manuscript 
close all; clear all; clc
%% read the audiofile with our synthesized sweeps 
[ref_sweep, fs] = audioread("sweep27s.wav");
ref_sweep = ref_sweep(fs+1:6*fs);
ref_sweep_ = ref_sweep(1:220416);

%% create signals to be contaminated
sweep_con = repmat(ref_sweep_, [3,1]);
sweep_con_2 = repmat(ref_sweep_, [3,1]);
time = (0:length(sweep_con)-1)./fs;
%% crate non-stationary noise signals

%transients
click_1 = zeros(length(sweep_con), 1);
click_2 = zeros(length(sweep_con), 1);


click_1(1.8*fs : 2*fs) = 0.4;
click_2(13.5*fs : 13.7*fs) = 0.3;

% tonal disturbances
tone =  zeros(length(sweep_con), 1);
tone(7.3*fs: 8*fs) = 0.1*sin(2*pi*1000*time(7.3*fs: 8*fs));
tone_2 = tone;
tone_2(1*fs: 4*fs) = 0.2*sin(2*pi*2000*time(1*fs: 4*fs));
tone_2(12*fs: 14*fs) = 0.4*sin(2*pi*5600*time(12*fs: 14*fs));

% contaminate the sweeps
sweep_con = sweep_con  + click_1 + click_2 + tone;
sweep_con_2 = sweep_con_2  +  tone_2;
%%
seg_sweep = reshape(sweep_con, [length(ref_sweep_), 3]);
seg_sweep_2 = reshape(sweep_con_2, [length(ref_sweep_), 3]);
%% Mosaic-T
med_sweep_time = median(seg_sweep, 2);
%% parameters for STFT
overlap = 128;
window_length = 2*overlap;
n_time = floor((size(seg_sweep, 1)-overlap)/(window_length - overlap));
nfft = 2^10;
FF = logspace(log10(20),log10(fs/2),2*nfft);
%% calculate STFT of the contaminated series
[sweep_con_spec, f_, t_] = stft(sweep_con_2,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');
%% Mosaic-TF - calculate the STFT of each sweep 
for n = 1: size(seg_sweep_2,2)
    [sw_spec(:, :,n), f,t] = stft(seg_sweep_2(:,n),fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');
end

%% Mosiac-TF - calculate the median
med_sweep_tf = median(real(sw_spec), 3) + 1i*median(imag(sw_spec),3);
sw_med = real(istft( med_sweep_tf,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided'));
%% colors for plotting
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840; 0 0.298 0.6; 0 0.4 0.2;0.8 0 0.4 ; 0.5, 0.5, 0.5; 0.75, 0 ,0.75 ; 0,0,0];


%% plot the results
set(groot,'defaultAxesTickLabelInterpreter','latex');
fig1 = figure(1); clf
sb1 = subplot(2,4, 1:2);
plot(time, sweep_con, 'k'); hold on
plot(time, click_1 + click_2 + tone, 'r', 'Linewidth', 2)
ylim([-0.65 0.95])
set(gca,'Fontsize',12,'YTick', -0.6:0.3:0.9, 'Xtick', 0:3:15)
xlabel({'Time (s)'; '(a)'}, 'interpreter', 'latex');
ylabel('Signal value', 'interpreter', 'latex')
sb1.Position(1) = 0.05;
sb1.Position(3) = 0.44;

sb2 = subplot(2,4, 3);
plot(time(1: length(ref_sweep_)), med_sweep_time, 'k')
ylim([-0.65 0.95])
xlim([0 5])
set(gca,'Fontsize',12,'YTick', -0.6:0.3:0.9)
xlabel({'Time (s)'; '(b)'}, 'interpreter', 'latex');
sb2.Position(2) = sb1.Position(2);
sb2.Position(4) = sb1.Position(4);

sb2_ = subplot(2,4, 4); 
histogram(seg_sweep(2.3*fs: 3*fs,[1,3]), 'DisplayStyle','stairs', 'Normalization','count', 'LineWidth',1.5,'edgecolor', 'k', HandleVisibility='off')
hold on
histogram(seg_sweep(2.3*fs: 3*fs,:), 'DisplayStyle','stairs', 'Normalization','count', 'LineWidth',1.5,'edgecolor', 'r', HandleVisibility='off')
set(gca, 'YDir', 'normal','Yscale', 'log')
xlabel({'Sample value'; '(c)'}, 'interpreter', 'latex');
ylabel('Sample count', 'interpreter', 'latex')
ylim([1 10^6])

plot([inf inf],[inf inf],'k', 'LineWidth',1.5)
plot([inf inf],[inf inf], 'r','LineWidth',1.5)

xline(median(seg_sweep(2.3*fs: 3*fs, [1,3]), 'all'),'color', colors(1,:) , 'LineWidth',2)
xline(median(seg_sweep(2.3*fs: 3*fs, :), 'all'), 'r--',  'LineWidth',2)
lgd = legend('Sweeps 1, 3', 'All sweeps', 'Median, sweeps 1, 3','Median, all sweeps', 'Location','north', 'numcolumns', 1, 'interpreter', 'latex', 'Fontsize', 10);
set(gca,'Fontsize',12)
sb2_.Position(1) = 0.82;

sb3 = subplot(2,4, 5:6);
s = pcolor(t_,f_, db(sweep_con_spec));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-40 45])
s.HandleVisibility = 'off';

xlabel({'Time (s)'; '(d)'}, 'interpreter', 'latex'); ylabel('Frequency (kHz)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

hold on
% xline(1.8, 'r')
% xline(13.5, 'r')
plot(time(7.3*fs: 8*fs), time(7.3*fs: 8*fs)*0+1000, 'r', 'LineWidth',2)
plot(time(12*fs: 14*fs), time(12*fs: 14*fs)*0+5600, 'r', 'LineWidth',2)
plot(time(1*fs: 4*fs), time(1*fs: 4*fs)*0+2000, 'r', 'LineWidth',2)
sb3.Position(1) = 0.05;
sb3.Position(3) = 0.44;

%
sb4 = subplot(2,4, 7);
s = pcolor(t,f, db(med_sweep_tf));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-40 45])

s.HandleVisibility = 'off';
clb = colorbar;
set(get(clb,'Label'),'String','Energy (dB)', 'Fontsize',12,'interpreter', 'latex')
clb.TickLabelInterpreter = 'latex';
xlabel({'Time (s)', '(e)'}, 'interpreter', 'latex');% ylabel('Frequency (kHz)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')
sb4.Position(2) = sb3.Position(2);
sb4.Position(4) = sb3.Position(4);
sb4.Position(3) = sb2.Position(3);


sb4_ = subplot(2,4, 8); 
histogram(db(sw_spec(47, :, [1, 3])), 'DisplayStyle','stairs', 'Normalization','count', 'LineWidth',1.5,'edgecolor', 'k', HandleVisibility='off')
hold on
histogram(db(sw_spec(47, :, :)), 'DisplayStyle','stairs', 'Normalization','count', 'LineWidth',1.5,'LineStyle','-','edgecolor', 'r', HandleVisibility='off')
set(gca, 'YDir', 'normal','Yscale', 'log')
xlabel({'Bin magnitude (dB)'; '(f)'}, 'interpreter', 'latex');
ylabel('Bin count', 'interpreter', 'latex')

plot([inf inf],[inf inf],'k', 'LineWidth',1.5)
plot([inf inf],[inf inf], 'r','LineWidth',1.5, 'LineStyle','-')


xline(median(db(sw_spec(47, :, [1, 3])), 'all'),'color', colors(1,:) , 'LineWidth',2)
xline(median(db(sw_spec(47, :, :)), 'all'), 'r--',  'LineWidth',2)
lgd = legend('Sweeps 1, 3', 'All sweeps', 'Median, sweeps 1, 3','Median, all sweeps', 'Location','north', 'numcolumns', 1, 'interpreter', 'latex', 'Fontsize', 10);
ylim([1 10^5])
set(gca,'Fontsize',12, 'YTick', [1 100 100000])
sb4_.Position(1) = 0.82;
%%
fig1.Position(3) = 1200;
fig1.Position(4) = 570;

%% print figure
set(fig1,'Units','Inches');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[fig1.Position(3), fig1.Position(4)])
print(fig1,'Median_filter','-dpdf','-r600')
