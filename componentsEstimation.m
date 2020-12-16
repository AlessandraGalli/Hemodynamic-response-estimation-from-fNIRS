function [all, tot, Psig] = componentsEstimation(sig, distance_out, duration_hrf, nHRF, D1_t0, D1, D2, Fs, Npt)

Ts = 1/Fs;
max_shift = 2;
lm_sx = floor(Npt/2)-round(max_shift/Ts)-floor(duration_hrf/Ts/2);
lm_dx = floor(Npt/2)+round(max_shift/Ts)-floor(duration_hrf/Ts/2);
inizio_hrf = distance_out(1:nHRF);
fine_hrf = distance_out(1:nHRF) + round(duration_hrf*Fs);
centro_hrf = (inizio_hrf+fine_hrf)/2;

tot = zeros(size(sig));
for i_hrf = 1:nHRF  
    inizio(i_hrf) = centro_hrf(i_hrf)-floor(Npt/2);
    fine(i_hrf) = centro_hrf(i_hrf)+floor(Npt/2);
    sig_bl = sig(inizio(i_hrf):fine(i_hrf));
    residual = sig_bl';
    Ensig = rms(sig_bl);
    Psig(i_hrf) = Ensig;
    cont_iter = 1;
    D = [];
    hrf = zeros(size(sig_bl));
    while Psig(i_hrf) > 0.02 && cont_iter<50
        cont_iter = cont_iter+1;
        %% iter 1 - physiological noise
        ysin = abs(D1_t0'*residual);
        [~, loc_sin] = max(ysin);
        in_start = (loc_sin-1)*6+1;
        in_end = in_start + 6 - 1;
        Dsin = [D D1(:,(in_start:in_end))];
        s_sin = pinv(Dsin)*sig_bl';
        sig_estimated_sin = real(Dsin*s_sin);                       %stima andamento ECG
        residual_sin = sig_bl' - sig_estimated_sin;
        P_res_sin = rms(residual_sin);
        %% iter 2 - hemodynamic response
        yhrf = (((D2'*residual)));
        [~, loc_temp] = max((yhrf(lm_sx:lm_dx)));
        loc_hrf = loc_temp + lm_sx - 1;
        Dhrf = [D D2(:,loc_hrf)];
        s_hrf = pinv(Dhrf)*sig_bl';
        sig_estimated_hrf = real(Dhrf*s_hrf);                       %stima andamento ECG
        hrf_temp = real(D2(:,loc_hrf)*s_hrf(end));
        residual_hrf = sig_bl' - sig_estimated_hrf;
        P_res_hrf = rms(residual_hrf);
        
        if P_res_hrf < 100 || P_res_hrf < P_res_sin 
            Psig(i_hrf) = P_res_hrf/Ensig;
            residual = residual_hrf;
            hrf = hrf_temp';
            break; 
        else 
             D = Dsin;
             Psig(i_hrf) = P_res_sin/Ensig;
             residual = residual_sin;
        end 
    end 
    
    all(i_hrf,:) = hrf;
    tot(inizio(i_hrf):fine(i_hrf)) = tot(inizio(i_hrf):fine(i_hrf))+hrf;
end 
end 