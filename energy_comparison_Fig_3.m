%% compare Mosiac-TF to averaging 
% in terms of energy loss due to time variance
% complementary code for the publication 
% "Non-stationary Noise Removal from Repeated Sweep Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to JASA Express Letters
% on 30.04.2024
%% plots figure 3 from the manuscript 
clear all; close all; clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
%% read measured sweeps
[measured_sweep, fs] = audioread('3_sweeps_energy_comparison.wav');
%% some sweep pamaeters
numSweep = 3;

%% STFT parameters
overlap = 256;
window_length = 2*overlap;
n_time = floor((size(measured_sweep, 1)-overlap)/(window_length - overlap));
nfft = 2^10;
FF = logspace(log10(1),log10(fs/2),nfft);
%% calculate the STFT of the measured sweeps
sw_spec = zeros(nfft, n_time, numSweep);

for n = 1: numSweep
    [sw_spec(:, :,n), f,t] = stft(measured_sweep(:,n),fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');
end
%% conduct different non-statioanry noise removal methods
% Mosaic-TF    
med_spec = median(real(sw_spec), 3) + 1i* median(imag(sw_spec), 3);
sw_med = real(istft( med_spec,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided'));

% Mean in time domain
sw_mean = mean(measured_sweep(1:size(sw_med, 1),  :), 2);
mean_spec =stft(sw_mean,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');
%% signal energy
sig_energy = abs(measured_sweep.^2);                % single sweep energy
sig_energy = sig_energy(1: length(sw_med), :, :);
sig_energy_med = abs(sw_med.^2);                    % energy after Mosaic-TF
sig_energy_mean = abs(sw_mean.^2);                  % energy after averaging


med_energy = median(sig_energy, 2);                 % median energy, not energy of median
mean_energy = mean(sig_energy, 2);                  % mean of energy, not energy of mean


diff_energy_med =  10*log10(med_energy) - 10*log10(sig_energy_med);     % calculate energy difference
diff_energy_mean =  10*log10(med_energy) - 10*log10(sig_energy_mean);   

mfilt_diff_med = medfilt1(diff_energy_med, nfft, [], 1);                % median filter the energy difference to smoothen out the curves
mfilt_diff_mean = medfilt1(diff_energy_mean, nfft, [], 1);


%% plot the figure
col1 = [0, 0.4470, 0.7410];

tsw = (0:length(sig_energy)-1)./fs;
tsw_m = (0:length(sig_energy_med)-1)./fs;
fig2 = figure(2); clf
s1 = subplot(1,2,1);
plot(tsw, medfilt1(db(sig_energy(:, end)), nfft), 'k', 'LineWidth',3); hold on
plot(tsw_m, medfilt1(db(sig_energy_med), nfft),'-.', 'color', col1, 'LineWidth',1.2)
plot(tsw_m, medfilt1(db(sig_energy_mean), nfft), 'r--', 'LineWidth',1.2)
xlim([500/fs, 2.2])
xline(2, 'k--')
Xtix = get(gca, 'Xtick');
set(gca, 'Xtick', [0.011, Xtix], 'XTickLabel', [0 Xtix], 'FontSize', 12)
ylabel('Energy (dB)','Interpreter','latex')
legend('Single sweep', 'Mosaic-TF', 'Time-domain Averaging', 'location', 'south', 'numcolumns', 1, 'interpreter', 'latex',  'FontSize', 12)
xlabel('Time (s)','Interpreter','latex')

%%%%
s1.Position(1) = 0.05;
s1.Position(2) = 0.18;
s1.Position(3) = 0.45;
%%%%%
s2 = subplot(1,2,2);
hold on
plot(tsw_m, medfilt1(diff_energy_med, nfft),'-.', 'color', col1, 'LineWidth',2)
plot(tsw_m, medfilt1(diff_energy_mean, nfft), 'r--', 'LineWidth',2)
xlim([500/fs, 2.2])
xline(2, 'k--')
yline(0, 'k:')
ylim([-0.5 6.5])

box on
set(gca, 'Xtick', [0.011, Xtix], 'XTickLabel', [0 Xtix], 'FontSize', 12)
xlabel('Time (s)','Interpreter','latex')
ylabel('Energy difference (dB)','Interpreter','latex')
legend( 'Mosaic-TF', 'Time-domain Averaging', 'location', 'north', 'numcolumns', 1, 'interpreter', 'latex',  'FontSize', 12)

s2.Position(2) = s1.Position(2) ;s2.Position(3) = 0.4;
%%
fig2.Position(4) = 250;
fig2.Position(3) = 1200;
%% print figure
set(fig2,'Units','Inches');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[fig2.Position(3), fig2.Position(4)])
% print(fig2,'Energy_loss','-dpdf','-r600')

