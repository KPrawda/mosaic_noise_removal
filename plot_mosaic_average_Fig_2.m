%% compare Mosaic-T and Mosiac-TF to averaging 
% in terms of non-stationary noise removal 
% complementary code for the publication 
% "Non-stationary Noise Removal from Repeated Sweep Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to JASA Express Letters
% on 30.04.2024
%% plots figure 2 from the manuscript 
clear all; close all; clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
%% read measured sweeps
[measured_sweep, fs] = audioread('.\Case - measured sweeps - Fig.2\3_sweeps_speech_transient.wav');
%% some sweep pamaeters
numSweep = 3;
durSweep = 2*fs;
nMedian = 3;

%% divide the series into singular sweeps
sw1 = zeros(round(durSweep*1.1), numSweep);

for n = 1:numSweep
    sw1(:, n) = measured_sweep((n-1)*durSweep+1:round(n*durSweep+0.1*durSweep));
end

%% upsample for time-aligning of sweeps
upFactor = 10; % upsampling factor
sw_up = zeros(size(sw1, 1)*upFactor, numSweep);

% upsampling
for n = 1: numSweep
    sw_up(:, n)  = interp(sw1(:, n) , upFactor);
end
%% time-align sweeps
sw_up_shifted = sw_up;

% cross-correlation to find the shift between sweeps
[c, lags] = xcorr(sw_up, 'normalized');
[val, loc] = max(c(:, 1:nMedian));
loc = loc - loc(1); % sweep 1 is the reference (arbitrary choice)

% shift the sweeps by an appropriate number of samples
for it = 2:3
        sw_up_shifted(:, it) =circshift(sw_up(:, it), loc(it));
end
%% decimate
sw_sh = zeros(size(sw1, 1), numSweep);

for pos = 1:size(sw1, 3)
    for n = 1: numSweep
        sw_sh(:,n, pos) = decimate(sw_up_shifted(:,n, pos), upFactor);
    end
end
%% STFT parameters
overlap = 256;
window_length = 2*overlap;
n_time = floor((size(sw1, 1)-overlap)/(window_length - overlap));
nfft = 2^10;
FF = logspace(log10(1),log10(fs/2),nfft);

%% calculate the STFT of the sweeps
sw_spec = zeros(nfft, n_time, numSweep);

for n = 1: numSweep
    [sw_spec(:, :,n), f,t] = stft(sw_sh(:,n),fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');
end

%% perform different non-stationary noise removal techniques
% Mosaic-TF
med_spec = median(real(sw_spec), 3) + 1i* median(imag(sw_spec), 3);
sw_med = real(istft( med_spec,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided'));

% Mosaic-T
sw_med_time= median(sw_sh(1:size(sw_med, 1),  :), 2);
med_spec_time = stft(sw_med_time,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');

% Mean
sw_mean = mean(sw_sh(1:size(sw_med, 1),  :), 2);
mean_spec =stft(sw_mean,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');
%% plot the figure
col1 = [0, 0.4470, 0.7410];

fig = figure(1); clf

%%%%%%%%%%%%%%%%%%%%%%%
sb1 = subplot(2,3,1);

s = pcolor(t,f,db(sw_spec(:, :, 1)));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-60 0])



ylabel('Frequency (kHz)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')


hold on
rectangle('Position', [0.3 200 0.6 10500], 'EdgeColor','r', 'LineStyle','--', 'LineWidth',2)
xlabel({'(a)'}, 'interpreter', 'latex');

sb1.Position(1) = 0.05;
sb1.Position(3) = 0.24;
%%%%%%%%%%%%%
sb2 = subplot(2,3,2);

s = pcolor(t,f,db(sw_spec(:, :, 2)));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-60 0])

hold on
rectangle('Position', [1 800 0.2 15000], 'EdgeColor',col1, 'LineStyle','-.', 'LineWidth',2)
rectangle('Position', [2 2000 0.15 15000], 'EdgeColor',col1, 'LineStyle','-.', 'LineWidth',2)

xlabel({'(b)'}, 'interpreter', 'latex');

set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

sb2.Position(1) =sb2.Position(1) -0.07;
sb2.Position(3) = sb1.Position(3);

%%%%%%%%%%%%%
sb3 = subplot(2,3,3);

s = pcolor(t,f,db(sw_spec(:, :, 3)));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-60 0])

hold on
rectangle('Position', [2 25 0.19 15000], 'EdgeColor',col1, 'LineStyle','-.', 'LineWidth',2)
rectangle('Position', [0 2000 0.2 9000], 'EdgeColor',col1, 'LineStyle','-.', 'LineWidth',2)
rectangle('Position', [0.4 500 0.4 12000], 'EdgeColor',col1, 'LineStyle','-.', 'LineWidth',2)
rectangle('Position', [1.1 1000 0.2 9000], 'EdgeColor',col1, 'LineStyle','-.', 'LineWidth',2)

xlabel({ '(c)'}, 'interpreter', 'latex');

clb = colorbar;
set(get(clb,'Label'),'String','Magnitude (dB)', 'Fontsize',12,'interpreter', 'latex')
clb.TickLabelInterpreter = 'latex';

set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

sb3.Position(1) =sb3.Position(1) -0.07;
sb3.Position(3) = sb1.Position(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
sb4 =  subplot(2,3,4);
s = pcolor(t,f,db(mean_spec));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-60 0])


xlabel({'Time (s)', '(d)'}, 'interpreter', 'latex'); ylabel('Frequency (kHz)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

sb4.Position(1) =sb1.Position(1);
sb4.Position(3) = sb1.Position(3);
sb4.Position(2)=0.13;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sb5 = subplot(2,3,5);
s = pcolor(t,f,db(med_spec_time));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-60 0])


xlabel({'Time (s)', '(e)'}, 'interpreter', 'latex'); 
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

sb5.Position(1) =sb5.Position(1) -0.07;
sb5.Position(3) = sb1.Position(3);
sb5.Position(2)=sb4.Position(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sb6 = subplot(2,3,6);
s = pcolor(t,f,db(med_spec));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-60 0])

s.HandleVisibility = 'off';
clb = colorbar;
set(get(clb,'Label'),'String','Magnitude (dB)', 'Fontsize',12,'interpreter', 'latex')
clb.TickLabelInterpreter = 'latex';

xlabel({'Time (s)', '(f)'}, 'interpreter', 'latex');
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

sb6.Position(1) =sb6.Position(1) -0.07;
sb6.Position(3) = sb1.Position(3);
sb6.Position(2)=sb4.Position(2);

%%
fig.Position(4) = 520;
fig.Position(3) = 1200;
%% print figure
set(fig,'Units','Inches');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[fig.Position(3), fig.Position(4)])
% print(fig,'Spectrograms_sweeps_new','-dpdf','-r600')

