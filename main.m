%-------------------------------%%--------------------------------%
%  Authors: Alessandra Galli, Sabrina Brigadoi, Giada Giorgi,     %
%              Giovanni Sparacino and Claudio Narduzzi.           %
%        Department of Information Engineering                    %
%                    University of Padova                         %
%                    galliale@dei.unipd.it                        %
%-------------------------------%%--------------------------------%
%% ------------------------------------------------------------------------
% Code related and described in the paper: 
% A. Galli et al. "Accurate Hemodynamic Response Estimation by Removal of 
% Stimulus-Evoked Superficial Response in fNIRS Signals." Journal of Neural
% Engineering.
% PLESE CITHE THIS PAPER UF YOU USE THE CODE. 
%% ------------------------------------------------------------------------

clc
clear all
close all

%% ------------------------------------------------------------------------
%%                        HOW  TO PREPARE THE DATA
%% ------------------------------------------------------------------------
% Prepare your fNIRS data in this way:
% Struct file named "raw" containing:
% - raw.triggers size Nx1 : with 1 for trigger t1 and 2 for trigger t2;
% - raw.data size Nx40 : measurement unit nM, odd columns (1:2:C-1) concentration of HbO (HbO-channels)
%   and even columns (2:2:C) concentration of HbR (HbR-channels);
% N = numbero of samples; C = number of channels.
%% ------------------------------------------------------------------------
%%                        START TUNING FROM THE USER
%% ------------------------------------------------------------------------
load subject12 %% INSERT HERE the name of the file to load

Fs = 7.8125;  % Sampling frequency.

ch_short_HbO = [13 15 33 35];                                    % short HbO channels
ch_long_HbO = [1 3 5 7 9 11 17 19 21 23 25 27 29 31 37 39];      % long HbO channels

ch_short_HbR = [14 16 34 36];                                    % short HbR channels
ch_long_HbR = [2 4 6 8 10 12 18 20 22 24 26 28 30 32 38 40];     % long HbR channels

nHRF = 39;                                                       % number of hr

task_selected = 2;                                               % select the task to analyze (1 or 2)

%% ------------------------------------------------------------------------
%%                         END TUNING FROM THE USER
%% ------------------------------------------------------------------------

t1 = find(raw.triggers == 1);     % time-instant of t1 task

t2 = find(raw.triggers == 2);     % time-instant of t2 task

if task_selected == 1
    t_sel = t1;
else
    t_sel = t2;
end

%%
Ts = 1/Fs;
L_time = 64;
Npt = round(L_time*Fs)+1;
data = raw.data;

load Dictionaries_Concentration
duration_hrf = 12;                                    % duration of hr
distance_out = t_sel;
n_punti_in = 400;
Nord = 100;
F_low = 0;
F_high = 0.3;
temp = fdesign.lowpass('N,Fc', Nord, F_high, Fs);
Filtro = design(temp);
HR_estimations = zeros(size(raw.data,2),Npt);
%% HbO
for ind = 1:length(ch_long_HbO)
    ind_l = ch_long_HbO(ind);
    fine_hrf = distance_out+ round(duration_hrf*Fs);
    centro_hrf = (distance_out+fine_hrf)/2;
    
    sig_long_HbO = data(:,ind_l)-mean(data(:,ind_l));
    sig_short_all_HbO = data(:,ch_short_HbO)-mean(data(:,ch_short_HbO));
    
    long_filt_HbO = filtfilt(Filtro.Numerator,1,sig_long_HbO);
    
    short_temp_HbO = zeros(size(sig_short_all_HbO));
    short_filt_all_HbO = zeros(size(sig_short_all_HbO));
    
    for i = 1:size(ch_short_HbO,2)
        short_temp_HbO(:,i) = filtfilt(Filtro.Numerator,1,sig_short_all_HbO(:,i));
        SS1 = short_temp_HbO(1:n_punti_in,i)';
        LS1 = long_filt_HbO(1:n_punti_in)';
        alpha1 = (LS1*SS1')/(SS1*SS1');
        short_filt_all_HbO(:,i)  = alpha1*short_temp_HbO(:,i)';
        C1(i) = corr(long_filt_HbO, short_filt_all_HbO(:,i));
    end
    [~, ind_sh_HbO] = max(C1);
    short_filt_HbO = short_filt_all_HbO(:,ind_sh_HbO);
    
    D1_t0 = D1_t0_type1(1:Npt,:);
    D1 = D1_type1(1:Npt,:);
    D2 = D2_type1(1:Npt,1:Npt);
    
    [noise_HbO, sig_noise_HbO, r_HbO(1,:)] = componentsEstimation(short_filt_HbO', distance_out, duration_hrf, nHRF, D1_t0, D1, D2, Fs, Npt);
    [tot_HbO, sig_tot_HbO, r_HbO(2,:)] = componentsEstimation(long_filt_HbO', distance_out, duration_hrf, nHRF, D1_t0, D1, D2, Fs, Npt);
    
    sig_HbO = sig_tot_HbO - sig_noise_HbO;
    
    all = tot_HbO - noise_HbO;
    
    n_HbO = 0;
    lim_min_HbO = 0;
    lim_max_HbO = 2.5*median(max(all'));
    ind_HbO = [];
    all_HbO = [];
    
    for j = 1: nHRF
        temp = all(j,:);
        [M,posM] = max(abs(temp));
        if M > lim_min_HbO && M < lim_max_HbO
            n_HbO = n_HbO + 1;
            all_HbO(n_HbO,:) = temp;
            ind_HbO =  [ind_HbO, j];
        end
    end
    
    if n_HbO > 1
        HR_estimations(ind_l,:) = mean(all_HbO);
    elseif n_HbO == 1
        HR_estimations(ind_l,:) = (all_HbO);
    else
        HR_estimations(ind_l,:) = zeros(1,Npt);
    end
end


%% HbR
for ind = 1:length(ch_long_HbR)
    ind_l = ch_long_HbR(ind);
    fine_hrf = distance_out+ round(duration_hrf*Fs);
    centro_hrf = (distance_out+fine_hrf)/2;
    sig_long_HbR = data(:,ind_l)-mean(data(:,ind_l));
    sig_short_all_HbR = data(:,ch_short_HbR)-mean(data(:,ch_short_HbR));
    long_filt_HbR = filtfilt(Filtro.Numerator,1,sig_long_HbR);
    
    short_temp_HbR = zeros(size(sig_short_all_HbR));
    short_filt_all_HbR = zeros(size(sig_short_all_HbR));
    
    for i = 1:size(ch_short_HbR,2)
        short_temp_HbR(:,i) = filtfilt(Filtro.Numerator,1,sig_short_all_HbR(:,i));
        SS1 = short_temp_HbR(1:n_punti_in,i)';
        LS1 = long_filt_HbR(1:n_punti_in)';
        alpha1 = (LS1*SS1')/(SS1*SS1');
        short_filt_all_HbR(:,i)  = alpha1*short_temp_HbR(:,i)';
        C1(i) = corr(long_filt_HbR, short_filt_all_HbR(:,i));
    end
    [~, ind_sh_HbR] = max(C1);
    short_filt_HbR = short_filt_all_HbR(:,ind_sh_HbR);
    
    D1_t0 = D1_t0_type2(1:Npt,:);
    D1 = D1_type2(1:Npt,:);
    D2 = D2_type2(1:Npt,1:Npt);
    
    [noise_HbR, sig_noise_HbR, r_HbR(1,:)] = componentsEstimation(short_filt_HbR', distance_out, duration_hrf, nHRF, D1_t0, D1, D2, Fs, Npt);
    [tot_HbR, sig_tot_HbR, r_HbR(2,:)] = componentsEstimation(long_filt_HbR', distance_out, duration_hrf, nHRF, D1_t0, D1, D2, Fs, Npt);
    
    sig_HbR = sig_tot_HbR - sig_noise_HbR;
    
    all = tot_HbR - noise_HbR;
    
    n_HbR = 0;
    lim_min_HbR = 0;
    lim_max_HbR = 2.5*median(max(abs(all')));
    ind_HbR = [];
    all_HbR = [];
    
    for j = 1: nHRF
        temp = all(j,:);
        [M,posM] = max(abs(temp));
        if M > lim_min_HbR && M < lim_max_HbR
            n_HbR = n_HbR + 1;
            all_HbR(n_HbR,:) = temp;
            ind_HbR =  [ind_HbR, j];
        end
    end
    if n_HbR > 1
        HR_estimations(ind_l,:) = mean(all_HbR);
    elseif n_HbR == 1
        HR_estimations(ind_l,:) = (all_HbR);
    else HR_estimations(ind_l,:) = zeros(1,Npt);
    end
end

% ESTIMATED HEMODYNAMIC RESPONSE SAVED IN: HR_estimations
