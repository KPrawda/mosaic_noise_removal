%% compare the non-stationary noise removal methods
% complementary code for the publication 
% "Non-stationary Noise Removal from Repeated Sweep Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to JASA Express Letters
% on 30.04.2024
%% plots figure 4 from the manuscript 
%% comparison to the following methods:
% [1]
% D. Ćirić, A. Pantić, and D. Radulović, “Transient noise effects in measurement of room impulse response 
% by swept sine technique,” in Proc. 10th Int. Conf. Telecommunication in Modern Satellite Cable and
% Broadcasting Services (TELSIKS) (2011), pp. 269–272

% [2]
% A. Telle and P. Vary, “A novel approach for impulse response measurements in environments 
% with time-varying noise,” in Proc. 20th Int. Cong. Acoust. (ICA) (2010), pp. 1–5.

% BONUS (not in the paper)
% [3]
% P. Massé, T. Carpentier, O. Warusfel, and M. Noisternig, “A robust denoising process for spatial room
% impulse responses with diffuse reverberation tails,” J. Acoust. Soc. Am. 147(4), 2250–2260 (2020) 
% doi:10.1121/10.0001070.312
% reasons for not including this in the paper:
% - the treatment of phase information is not mentioned in the paper
% - the \alpha parameter is not explained sufficiently 
% - averaging does NOT remove the contamination completely, thus this method will always leave some non-stationary noise in the sweep 

%% housekeeping
clear all; close all; clc
set(groot,'defaultAxesTickLabelInterpreter','latex');
%% read audiofiles
[sw, fs] = audioread("impulse_2nd_sweep.wav");  % read series of 5 sweeps from Arni with an impulse in 2nd sweep
ref_sw = audioread("sweep27s.wav");             % read the excitation sweep
%% split the 5 recorded sweeps into single signals
sw= sw(fs+1:end-fs);
sw_ = reshape(sw, [5*fs, 5]);
ref_sw_ = ref_sw(fs+1:6*fs);                    % cut the excitation signal to include 1 sweep only

%% calculate PCC between the clean (1st sweep) and the contaminated sweep (2nd sweep)
PCC_meas = corrcoef(sw_(:, 1), sw_(:, 2))
%% calculate the detection threshold (Rule of Two)
bx = 1.4826^2;
noise = sw_(1:fs,:);
noiseEnergy =  squeeze(bx*median(noise.^2));
sig_energy = mean(sw_.^2);
snr_thresh_sw = mean(sig_energy - 2*noiseEnergy)./mean(sig_energy)
%% procedure A from [1]
sw_A = sw_(:, 2);
rms_A = rms(sw_A(2*fs:2.6*fs));                 % contamination is here, RMS for normalization
rms_ref = rms(ref_sw_(2*fs:2.6*fs));
sw_A(2*fs:2.6*fs) = (rms_A/rms_ref)*ref_sw_(2*fs:2.6*fs); % swapping the contaminated part with the scaled part of the refeence sweep
%% 
PCC_A = corrcoef(sw_(:, 1), sw_A)
noiseA = [sw_(1:fs,1),sw_A(1:fs)] ;
noiseEnergyA =  squeeze(bx*median(noiseA.^2));
sig_energyA = mean([sw_(:,1),sw_A].^2);
snr_thresh_A = mean(sig_energyA - 2*noiseEnergyA)./mean(sig_energyA)
%% procedure B from [1]
sw_B =  sw_(:, 2);

% isolate parts of contaminated and reference sweeps
sw_contamination = sw_B(2*fs:2.6*fs);           % contamination is here
sw_ref_part = ref_sw_(2*fs:2.6*fs);

% calculate FFts
sw_ref_band_fft = fft(sw_ref_part);
sw_band_fft = fft(sw_contamination);

% multiply the magnitudes of contaminated and references 
new_mag = abs(sw_band_fft.*sw_ref_band_fft);

% reconstruct the new spectrum with new magnitude and keeping the phase
sw_new_fft = new_mag.*exp(1i*angle(sw_band_fft));

% swap the part of the contaminated sweep
sw_new_band = real(ifft(sw_new_fft));
sw_new_band = (sw_new_band./max(abs(sw_new_band)))*max(abs(sw_contamination));

sw_B(2*fs:2.6*fs) = sw_new_band;
%%
PCC_B = corrcoef(sw_(:, 1), sw_B)

noiseB = [sw_(1:fs,1),sw_B(1:fs)] ;
noiseEnergyB =  squeeze(bx*median(noiseB.^2));
sig_energyB = mean([sw_(:,1),sw_B].^2);
snr_thresh_B = mean(sig_energyB - 2*noiseEnergyB)./mean(sig_energyB)
%% procedure C from [1]

sw_C = sw_(:, 2);
% deconvolution - the RIR will serve as transient noise 
S = fft(ref_sw_);
X = fft(sw_(:,3));
H = X./S;
h = ifft(H); % deconvolved RIR

na = sw_(1:fs, :);
na = reshape(na, [5*fs,1]);
Na = fft(na);

X_C = fft(sw_C);
H_C = X_C./(S + H + Na);

X_C_ = H_C.*S;

sw_C_ = ifft(X_C_);

%%
PCC_C = corrcoef(sw_(:, 1), sw_C_)
noiseC = [sw_(1:fs,1),sw_C_(1:fs)] ;
noiseEnergyC =  squeeze(bx*median(noiseC.^2));
sig_energyC = mean([sw_(:,1),sw_C_].^2);
snr_thresh_C = mean(sig_energyC - 2*noiseEnergyC)./mean(sig_energyC)
%% weighted average from [2]
pow_i = 0.2*(sum(abs(sw_).^2, 1)); %signal power
c_i = 1./(pow_i/sum(pow_i)); % weigths
c_i = c_i./sum(c_i); % weight normalizatiopn to 1

sw_mean_w = 0;
for i = 1:5
    sw_mean_w = sw_mean_w + c_i(i)*sw_(:, i);
end
%%
PCC_mean_w = corrcoef(sw_(:, 1), sw_mean_w)

noise_w = [sw_(1:fs,1),sw_mean_w(1:fs)] ;
noiseEnergy_w =  squeeze(bx*median(noise_w.^2));
sig_energy_w = mean([sw_(:,1),sw_mean_w].^2);
snr_thresh_w = mean(sig_energy_w - 2*noiseEnergy_w)./mean(sig_energy_w)
%% spectrogram parameters
overlap = 256;
window_length = 2*overlap;
% n_time = floor((size(seg_sweep, 1)-overlap)/(window_length - overlap));
nfft = 2^10;
fs = 44100;
FF = logspace(log10(20),log10(fs/2),nfft);

%% spectrograms

[sw_spec, f,t] = stft(sw_,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');

    
sw_spec_med = median(real(sw_spec), 3) + 1i*median(imag(sw_spec), 3); % MOSAIC-TF: median in time and freq


sw_med_tf = real(istft( sw_spec_med,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided'));

sw_A_spec = stft(sw_A,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');
sw_B_spec = stft(sw_B,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');
sw_C_spec = stft(sw_C_,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');

sw_mean_w_spec = stft(sw_mean_w,fs,'Window',hamming(window_length),'OverlapLength',overlap, 'FFTlength', 2*nfft-1, 'FrequencyRange', 'onesided');

%%
PCC_TF = corrcoef(sw_(1:length(sw_med_tf), 1), sw_med_tf)

noise_tf = [sw_(1:fs,[1, 3:5]),sw_med_tf(1:fs)] ;
noiseEnergy_tf =  squeeze(bx*median(noise_tf.^2));
sig_energy_tf = mean([sw_(1:length(sw_med_tf),[1, 3:5]),sw_med_tf].^2);
snr_thresh_tf = mean(sig_energy_tf - 2*noiseEnergy_tf)./mean(sig_energy_tf)
%% plot
fig = figure(1); clf;

sb1 = subplot(2,3,1);
s = pcolor(t,f, db(sw_spec(:,:, 2)));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-50 30])
s.HandleVisibility = 'off';


xlabel({'(a)'}, 'interpreter', 'latex'); ylabel('Frequency (kHz)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')
annotation('textbox','Position',  [0.052, 0.85, 1, 0.1], 'String', "PCC = 0.99830"  ,'FitBoxToText','on', 'BackgroundColor', 'w', 'Interpreter','latex', 'EdgeColor','none')

hold on
rectangle('position', [2,  25, 0.7 20000],'linestyle', '--', 'edgecolor', [0 0.4470 0.7410], 'linewidth', 2)

sb1.Position(1) = 0.05;
sb1.Position(2) = 0.65;
sb1.Position(3) = 0.24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb2 = subplot(2,3,2);
s = pcolor(t,f, db(sw_spec_med));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-50 30])

set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

xlabel({'(b)'}, 'interpreter', 'latex');

sb2.Position(1) = sb2.Position(1) -0.07;
sb2.Position(3) = sb1.Position(3);
sb2.Position(2) =sb1.Position(2);

hold on
annotation('textbox','Position',  [0.342, 0.85, 1, 0.1], 'String', "PCC = 0.99981"  ,'FitBoxToText','on', 'BackgroundColor', 'w', 'Interpreter','latex', 'EdgeColor','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb3 = subplot(2,3,3);
s = pcolor(t,f, db(sw_mean_w_spec));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-50 30])

clb = colorbar;
set(get(clb,'Label'),'String','Magnitude (dB)', 'Fontsize',12,'interpreter', 'latex')
clb.TickLabelInterpreter = 'latex';
clb.Label.FontSize = 12;

xlabel({'(c)'}, 'interpreter', 'latex');

set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

sb3.Position(1) = sb3.Position(1) -0.07;
sb3.Position(3) = sb1.Position(3);
sb3.Position(2) =sb1.Position(2);

hold on
annotation('textbox','Position',  [0.625, 0.85, 1, 0.1], 'String', "PCC = 0.99979"  ,'FitBoxToText','on', 'BackgroundColor', 'w', 'Interpreter','latex', 'EdgeColor','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb4 = subplot(2,3,4);
s = pcolor(t,f, db(sw_A_spec));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-50 30])

sb4.Position(1) = 0.05;
sb4.Position(3) = sb1.Position(3);
sb4.Position(2)= 0.20;

annotation('textbox','Position',  [0.052, 0.40, 1, 0.1], 'String', "PCC = 0.46423"  ,'FitBoxToText','on', 'BackgroundColor', 'w', 'Interpreter','latex', 'EdgeColor','none')


xlabel({'Time (s)'; '(d)'}, 'interpreter', 'latex');ylabel('Frequency (kHz)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb5 = subplot(2,3,5);
s = pcolor(t,f, db(sw_B_spec));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-50 30])

xlabel({'Time (s)'; '(e)'}, 'interpreter', 'latex');
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

sb5.Position(1) = sb2.Position(1);
sb5.Position(3) = sb1.Position(3);
sb5.Position(2)=   sb4.Position(2) ;

annotation('textbox','Position',  [0.342, 0.40, 1, 0.1], 'String', "PCC = 0.99069"  ,'FitBoxToText','on', 'BackgroundColor', 'w', 'Interpreter','latex', 'EdgeColor','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb6 = subplot(2,3,6);
s = pcolor(t,f, db(sw_C_spec));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-50 30])

clb = colorbar;
set(get(clb,'Label'),'String','Magnitude (dB)', 'Fontsize',12,'interpreter', 'latex')
clb.TickLabelInterpreter = 'latex';
clb.Label.FontSize = 12;

xlabel({'Time (s)'; '(f)'}, 'interpreter', 'latex');
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

sb6.Position(1) = sb3.Position(1);
sb6.Position(3) = sb3.Position(3);
sb6.Position(2)=   sb4.Position(2) ;

annotation('textbox','Position',  [0.625, 0.40, 1, 0.1], 'String', "PCC = 0.99830"  ,'FitBoxToText','on', 'BackgroundColor', 'w', 'Interpreter','latex', 'EdgeColor','none')

%%
fig.Position(3) = 1200;
fig.Position(4) = 500;


%% BONUS - average from [3]

sw_spec_mean = mean(abs(sw_spec), 3); % mean magnitude spectrogram
sw_spec_std = std(abs(sw_spec),0, 3); % standard deviation magnitude spectrogram

alpha = 1;                            % deviation factor, values not mentioned in the paper so very hand-tuned parameter

max_deviation = sw_spec_mean + alpha * sw_spec_std; % calculating max deviation

mag_contaminated = abs(sw_spec(:,:,2));

sw_cont_deviation =  mag_contaminated-sw_spec_mean; % the deviation between the contaminated sweep and average sweep

sum(sw_cont_deviation > max_deviation, 'all') % does it work?

for i = 1: size(sw_spec, 1)
    for j = 1: size(sw_spec, 2)
        if sw_cont_deviation(i,j) > max_deviation(i,j)
            mag_contaminated(i,j) = sw_spec_mean(i,j);
        end
    end
end
            

%% plot
fig2 = figure(2); clf;

sb21 = subplot(1,3,1);
s = pcolor(t,f, db(sw_spec(:,:, 2)));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-50 30])
s.HandleVisibility = 'off';


xlabel({'(a)'}, 'interpreter', 'latex'); ylabel('Frequency (kHz)', 'interpreter', 'latex')
set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')
% annotation('textbox','Position',  [0.052, 0.85, 1, 0.1], 'String', "PCC = 0.99830"  ,'FitBoxToText','on', 'BackgroundColor', 'w', 'Interpreter','latex', 'EdgeColor','none')

hold on
% rectangle('position', [2,  25, 0.7 20000],'linestyle', '--', 'edgecolor', [0 0.4470 0.7410], 'linewidth', 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb22 = subplot(1,3,2);
s = pcolor(t,f, db(sw_spec_med));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-50 30])

set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')

xlabel({'(b)'}, 'interpreter', 'latex');



hold on
% annotation('textbox','Position',  [0.342, 0.85, 1, 0.1], 'String', "PCC = 0.99981"  ,'FitBoxToText','on', 'BackgroundColor', 'w', 'Interpreter','latex', 'EdgeColor','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb23 = subplot(1,3,3);
s = pcolor(t,f, db( mag_contaminated));
s.EdgeColor = 'none';
set(gca, 'YDir', 'normal','Yscale', 'log')
colormap(flipud(gray))
caxis([-50 30])

clb = colorbar;
set(get(clb,'Label'),'String','Magnitude (dB)', 'Fontsize',12,'interpreter', 'latex')
clb.TickLabelInterpreter = 'latex';
clb.Label.FontSize = 12;

xlabel({'(c)'}, 'interpreter', 'latex');

set(gca,'Fontsize',12,'YTick', [ 10 100  1000  10000], 'YTicklabel',{ '0.01', '0.1'  '1'  '10'})
box on
set(gca, 'layer','top')


hold on
% annotation('textbox','Position',  [0.625, 0.85, 1, 0.1], 'String', "PCC = 0.99979"  ,'FitBoxToText','on', 'BackgroundColor', 'w', 'Interpreter','latex', 'EdgeColor','none')


%%
fig2.Position(3) = 1200;
fig2.Position(4) = 250;