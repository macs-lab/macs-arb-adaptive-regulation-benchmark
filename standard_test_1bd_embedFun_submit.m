% Adaptive multiple narrow-band disturbance rejection applied to an active
% suspension.
% 1 narrow band disturbance
% Benchmark project by prof. Ioan Landau.
% ============================================================
%   Copyright (c) 2008-, Xu Chen, University of Washington
%   Author(s): Xu Chen
%              University of Washington
% ============================================================
% init version: 2012-08-03
clear all
close all
%% Define Constants
FLAG_CONST_DIST_FREQ        = 1;
SW_EXPERIMENT               = 1;
FLAG_STEP_CHANGE_DIST_FREQ  = 0;
FLAG_CHIRP_DIST             = 2;
SW_ONE_SIMU_TEST            = 1; % test a single frequency by default
SW_DIST_ON                  = 1; % default turn on const or step change disturbance frequency
SW_CHIRP_DIST               = 0; % default turn off chirp disturbance
SW_TUNE                     = 0;
SW_BASELINE_CONTROL_SYS     = 0; % check the baseline system
SW_ADDITIONAL_PLOT          = 0;
SW_STEADY_STATE_CONTROL_SYS = 0;
Fs=800;     Ts=1/Fs;        Te=Ts;

bode_opt            = bodeoptions;
bode_opt.FreqUnits  = 'Hz';
bode_opt.FreqScale  = 'Linear';
bode_opt.xlim       = {[0 400]};
bode_opt.PhaseWrapping = 'On';
%%
% SELECT TEST OPTIONS HERE each time the test is run:
%       FLAG_DIST_FREQ = 0 ------ step changing disturbance frequency
%       FLAG_DIST_FREQ = 1 ------ constant disturbance frequency
%       FLAG_DIST_FREQ = 2 ------ chirp disturbance
disp('=============Multiple Narrow Band Disturbance Rejection============')
disp('=============       Adaptive Disturbance Observer      ============')
disp('===================================================================')
disp('SELECT TEST OPTIONS:')
disp('0 (default) ---- step changing disturbance frequency')
disp('1           ---- constant disturbance frequency')
disp('2           ---- chirp disturbance')
disp('   ')
disp('Press ENTER for default selection.')
while 1
    FLAG_DIST_FREQ = input(':');
    if isempty(FLAG_DIST_FREQ)
        FLAG_DIST_FREQ      = 0;
    end
    if FLAG_DIST_FREQ ~= 0 && FLAG_DIST_FREQ ~= 1 && FLAG_DIST_FREQ ~= 2
        disp('Wrong selection. Please re-select.')
    else
        break;
    end
end
% in case nothing selected
if ~exist('FLAG_DIST_FREQ','var')
    FLAG_DIST_FREQ = FLAG_STEP_CHANGE_DIST_FREQ;
end

disp('===================================================================')
disp('CHOOSE THE TEST LENGTH:')
disp('1 (default) ---- a quick sample test')
disp('0           ---- the entire frequency profile specified by the benchmark')
disp('   ')
disp('Press ENTER for default selection.')
while 1
    SW_ONE_SIMU_TEST = input(':');
    if isempty(SW_ONE_SIMU_TEST)
        SW_ONE_SIMU_TEST      = 1;
    end
    if SW_ONE_SIMU_TEST ~= 0 && SW_ONE_SIMU_TEST ~= 1
        disp('Wrong selection. Please re-select.')
    else
        break;
    end
end

disp('===================================================================')
disp('CHOOSE WHETHER OR NOT TO SAVE THE TEST DATA.')
while 1
    SW_SAVE_DATA = input('Save the test result?\n   1(default, press ENTER): yes\n   0: no\n:');
    if isempty(SW_SAVE_DATA)
        SW_SAVE_DATA        = 1;
    end
    if SW_SAVE_DATA ~= 0 && SW_SAVE_DATA ~= 1
        disp('Wrong selection. Please re-select.')
    else
        break;
    end
end

disp('===================================================================')
disp('ADAPTATION SCHEME.')
while 1
    SW_UNIFORM_ADAP = input('Uniform adaptation gain? (more conservative performance)\n   1(default, press ENTER): yes\n   0: no\n:');
    if isempty(SW_UNIFORM_ADAP)
        SW_UNIFORM_ADAP        = 1;
    end
    if SW_UNIFORM_ADAP ~= 0 && SW_UNIFORM_ADAP ~= 1
        disp('Wrong selection. Please re-select.')
    else
        break;
    end
end

FLAG_PERFORMANCE_EVAL       = 0;
%%
NBn                         = 1; % number of narrow bands
% chirp distrubance parameters
chirp_dist.level            = 1;
chirp_dist.freq1_seq        = [50];
chirp_dist.freq2_seq        = [95];
chirp_dist.para1            = [chirp_dist.freq1_seq(1);0;0];
chirp_dist.para2            = [chirp_dist.freq2_seq(1);0;0];

chirp_dist.freq1_init_time  = 5;
chirp_dist.chirp_init_time  = 10;
chirp_dist.chirp_dur_time   = 4.5;
chirp_dist.freq2_dur_time   = 5;

% load band_pass_filter_50To95 % 2010-09-26
% denoise_filter = tf(BP_ss_simulink)*tf(BP_ss_simulink);%2012-08-03
load band_pass_coef_50To95
denoise_filter = tf(conv(numBP,numBP),conv(denBP,denBP),Ts);
%% Adaptation parameters
% forgetting factor
adap_method = 2;
% initialize estimation at 0 Hz
theta1_init = -2*cos(0*2*pi*Ts);
theta_init  = theta1_init;
F1          = 2000;
F_init      = F1;

% for exponentially increasing forgetting factor
lambda_init = 0.93;
lambda_end  = 0.99;
%% Band-pass Q filter parameter
alpha_init  = 0.865; 
alpha_end   = 0.8650;
adap_init.alpha_pre = 0.865;
alpha       = alpha_end;

adap_init.alpha_init    = alpha_init;
adap_init.alpha_end     = alpha_end;
adap_init.F             = F_init;
adap_init.theta         = theta_init;
adap_init.lambda_init   = lambda_init;
adap_init.lambda_end    = lambda_end;
adap_init.SW_lambda     = adap_method;
adap_init.SW_2Stage     = 0;
%% define_plant_controllers
%/////////// primary path sys tf
load model_prim.mat Bp Ap %numerator and denominator of the primary path
%/////////// closed loop R/S controller
load RS_contr_sec R S
%/////////// secondary path sys tf
load model_sec.mat B A %numerator and denominator of the secondary path

load hinf_inv_landau_coef
P_inv = tf(numINVP,denINVP,Ts);
%/////////// loading noise values
load bruitbench

if SW_BASELINE_CONTROL_SYS
    figure;
    bodeplot(tf(R,S,Ts),bode_opt)
    grid on,zoom on
    title('Frequency response of the feedback controller')
    S_func = feedback(1,tf(R,S,Ts,'variable','z^-1')*tf(B,A,Ts,'variable','z^-1'));
    figure, bodeplot(S_func,bode_opt)
    figure, bodeplot(tf(B,A,Ts,'variable','z^-1'),bode_opt)

    L       = length(bruitbench);
    NFFT    = 2^nextpow2(L);
    [spec_bruitbench.f,spec_bruitbench.amp] =...
        spectre_psd_rms(bruitbench,Fs,NFFT);
    figure;
    plot(spec_bruitbench.f,spec_bruitbench.amp)
    xlabel('Frequency [Hz]')
    ylabel('dB [Vrms]')
    title('Spectral density of the measurement noise')
end

simuName    = 'simulator_1bd_submit';
%% Narrow band disturbances define
% Frequencies to be tested
if SW_ONE_SIMU_TEST % testing a single frequency
    freq_test   = 55; % change to desired region if required
else                % testing the entire frequency
    freq_test   = 50:5:95;
end % SW_ONE_SIMU_TEST
% Define the figure numbers
FIG_NUMBER_STEP_CHANGE_DIST     = [120;121;122];
FIG_NUMBER2_STEP_CHANGE_DIST    = [125;126;127];
FIG_NUMBER_CONST_DIST_FREQ      = [100:100-1+length(freq_test)];
FIG_NUMBER2_CONST_DIST_FREQ     = [150:150-1+length(freq_test)];
FIG_NUMBER3_CONST_DIST_FREQ     = [250:250-1+length(freq_test)];
FIG_NUMBER_CHIRP_DIST           = [200:202];
FIG_NUMBER2_CHIRP_DIST          = [300:302];
%% Run the test
if SW_UNIFORM_ADAP
    adap_init.SW_2Stage     = 1;
    adap_init.SW_lambda     = 0;
    lambda_end              = 0.99;
    adap_init.alpha_pre     = 0.92;
    adap_init.lambda_end    = lambda_end;
        
    load NF98_2;    
    numNF98 = numNF98_2;
    denNF98 = denNF98_2;
    load NF46;
    load NF98;
end
%% CONSTANT UNKONWN DISTURBANCE FREQUENCY
if FLAG_DIST_FREQ == FLAG_CONST_DIST_FREQ
    adap_init.SW_2Stage     = 1;
    adap_init.SW_lambda     = 0;
    if ~SW_UNIFORM_ADAP
        lambda_end              = 1;
    end
    adap_init.alpha_pre     = 0.92;
    adap_init.lambda_end    = lambda_end;
    data_cont_freq.readme = 'stores data in the case of constant disturbance frequencies';
    data_cont_freq.y{1,1} = 'openLoop';
    data_cont_freq.y{1,2} = 'closedLoop';
    for ii = 1:length(freq_test)
        NBw         = freq_test(ii)*2*pi;
        % true parameters (for result comparision later)
        lb1true     = 2*cos(NBw(1)*Ts);
        theta1_true = -lb1true;
        % simulink parameter definition
        % Narrow band disturbance injection time
        t_NBon      = 5;
        % compensation turn on time
        t_Qon       = t_NBon;
        t_AdapOn    = t_NBon;
        % narrow band disturbance duration
        t_NBdur     = 15/5;
        t_AdapOff   = t_AdapOn+t_NBdur*5;
        % 2010-09-26; For disturbance generator v2
        t_dur_lastDist = t_NBdur;
        
        % freq in Hz
        NBf         = NBw/2/pi;
        % simulink time
        t_sim       = 30;
        f           = NBf; f1 = f;
        dist_seq1   = [NBf; 0; 0];
        dist_seq2   = [NBf; 0; 0];
        dist_seq3   = [NBf; 0; 0];
        
        for jj = 1:2
            if jj == 1
                % run the test without compensation
                SW_COMP_ON = 0; % compensation off
                SW_CLOSE_LOOP = 0;
                % define figure names for later use
                fig_name_spec_residule = ...
                    ['level1_spec_residule_',...
                    num2str(freq_test(ii)),...
                    'Hz_openLoop'];
                fig_name_residule_time_trace = ...
                    ['level1_time_trace_residule_',...
                    num2str(freq_test(ii)),...
                    'Hz_openLoop'];
            else
                % run the test with compensation
                SW_COMP_ON = 1; % compensation on
                SW_CLOSE_LOOP = 1;
                
                fig_name_spec_residule = ...
                    ['level1_spec_residule_',...
                    num2str(freq_test(ii)),...
                    'Hz_closedLoop'];
                fig_name_residule_time_trace = ...
                    ['level1_time_trace_residule_',...
                    num2str(freq_test(ii)),...
                    'Hz_closedLoop'];
            end
            %% simulation start
            sim(simuName)
            data_cont_freq.narrow_band_freq(ii,:) = NBf';
            data_cont_freq.y{ii+1,jj} = y;
            if ii == 1 && jj == 1
                if SW_SAVE_DATA % creat the result folder
                    if ~exist(['level1_test_result_',date],'dir')
                        mkdir(['level1_test_result_',date]);
                    end
                end
            end
            level1_test_data_analysis_submit;
            if SW_SAVE_DATA
                try
                    movefile('*.fig',['level1_test_result_',date])
                catch
                end
            end
        end
        pause(10); % let the CPU take a 10-sec rest
    end
    disp ('===================================================================')
    disp ('test results saved to: level1_freq_dom_result_const_freq')
    disp ('                       level1_time_dom_result_const_freq')
    disp ('raw data saved to:     data_cont_freq')
    if SW_SAVE_DATA
        save(['level1_test_result_',date,'\level1_time_dom_result_const_freq'],...
            'level1_time_dom_result_const_freq');
        save(['level1_test_result_',date,'\level1_freq_dom_result_const_freq'],...
            'level1_freq_dom_result_const_freq');
        save (['level1_test_result_',date,'\data_cont_freq'],...
            'data_cont_freq');
    end
    
    %% STEP CHANGE DISTURBANCE FREQUENCY
elseif FLAG_DIST_FREQ == FLAG_STEP_CHANGE_DIST_FREQ
    t_NBon      = 5;               % NB Dist injection time
    t_Qon       = t_NBon;          % bandpass Q filter on time
    t_AdapOn    = t_NBon;
    t_NBdur     = 3;
    t_AdapOff   = t_AdapOn+5*t_NBdur;
    t_sim       = 40;
    t_dur_lastDist = 3;
    if ~SW_UNIFORM_ADAP
        adap_init.SW_2Stage     = 0;
        lambda_end              = 0.999;
        adap_init.lambda_end    = lambda_end;
    end
    % the center frequencies of the three step changing disturbance
    % frequencies
    center_freq = [60,75,85];
    % the three step changing disturbance frquencies (in Hz)
    freq_seq1   = [60,70,60,50,60];
    freq_seq2   = [75,85,75,65,75];
    freq_seq3   = [85,95,85,75,85];
    freq_table  = [freq_seq1; freq_seq2; freq_seq3];
    
    data_step_freq.readme = 'stores the data for the case of step changing disturbance frequencies';
    data_step_freq.y{1,1} = 'openLoop';
    data_step_freq.y{1,2} = 'closedLoop';
    if SW_ONE_SIMU_TEST % perform just one test
        ITER_STEP = 1;
    else
        ITER_STEP = 3;
    end
    for ii = 1:ITER_STEP
        % define the three sets of test sequences
        if ii == 1
            dist_seq1   = [60; 0; 0];
            dist_seq2   = [70; 0; 0];
            dist_seq3   = [50; 0; 0];
        elseif ii == 2
            dist_seq1   = [75; 0; 0];
            dist_seq2   = [85; 0; 0];
            dist_seq3   = [65; 0; 0];
        else
            dist_seq1   = [85; 0; 0];
            dist_seq2   = [95; 0; 0];
            dist_seq3   = [75; 0; 0];
        end
        for jj = 1:2
            if jj == 1
                SW_COMP_ON = 0; % compensation off
                SW_CLOSE_LOOP = 0;
                fig_name_residule_time_trace = ...
                    ['level1_time_trace_residule_center_freq_',...
                    num2str(center_freq(ii)),...
                    'Hz_openLoop'];
            else
                SW_COMP_ON = 1; % compensation on
                SW_CLOSE_LOOP = 1;
                fig_name_residule_time_trace = ...
                    ['level1_time_trace_residule_center_freq_',...
                    num2str(center_freq(ii)),...
                    'Hz_closedLoop'];
            end
            %% simulation start
            sim(simuName)
            data_step_freq.initial_freq(ii,:) = dist_seq1';
            data_step_freq.y{ii+1,jj} = y;
            if ii == 1 && jj == 1
                if SW_SAVE_DATA
                    if ~exist(['level1_test_result_',date],'dir')
                        mkdir(['level1_test_result_',date]);
                    end
                end
            end
            %% time domain result
            if jj == 2
                level1_time_dom_result_step_change_freq.readme =...
                    {'row: each row represents one step-changing disturbance sequence',...
                    'col: from column 1 to column 5: transient for the 1st to 5th dist '};
                for kk = 1:5
                    level1_time_dom_result_step_change_freq.transi_norm_square_3sec(ii,kk) = ...
                        sum(...
                        y.signals.values(...
                        (t_NBon+(kk-1)*t_NBdur)/Ts+1 :...
                        (t_NBon+(kk-1)*t_NBdur+3)/Ts ).^2);
                    
                    level1_time_dom_result_step_change_freq.max_residule(ii,kk) =...
                        max(...
                        y.signals.values(...
                        (t_NBon+(kk-1)*t_NBdur)/Ts :...
                        (t_NBon+(kk-1)*t_NBdur+1)/Ts )...
                        );
                    [temp_TD, temp_maximum, temp_transient_norm_square] =...
                        transient_duration_step_changes_additional_output(...
                        y.signals.values(1:t_AdapOff*Fs),...
                        Fs,t_AdapOff,t_NBon,t_NBdur,t_NBdur,'PlotOff',kk);
                    level1_time_dom_result_step_change_freq.t_transient(ii,kk) =...
                        temp_TD;
                    level1_time_dom_result_step_change_freq.transient_norm_square(ii,kk) =...
                        temp_transient_norm_square;
                    level1_time_dom_result_step_change_freq.maximum_transient(ii,kk) =...
                        temp_maximum;
                end
                h = figure;
                plot(y.time,y.signals.values);grid;
                xlabel('Time [sec]');ylabel('Residual force [V]');
                figure_specific
                if SW_SAVE_DATA
                    hgsave(h,fig_name_residule_time_trace,'-v6')
                end
            end
            
            h = figure(FIG_NUMBER_STEP_CHANGE_DIST(ii));
            grid on;hold on;
            if jj == 1
                plot(y.time,y.signals.values,'r');
            else
                plot(y.time,y.signals.values,'k:');
                legend('open loop','closed loop')
                xlabel('Time [sec]');ylabel('Residual force [V]');
                figure_specific
                if SW_SAVE_DATA
                    hgsave(h,['level1_time_trace_residule_center_freq_',...
                        num2str(center_freq(ii)),...
                        'Hz_compare'],'-v6')
                end
            end
            
            h = figure(FIG_NUMBER2_STEP_CHANGE_DIST(ii));
            hold on;
            if jj == 1
                subplot(211)
                plot(y.time,y.signals.values,'r');
                legend 'Open loop';
                ylabel('Residual force [V]');
                figure_specific
            else
                subplot(212)
                plot(y.time,y.signals.values,'k');
                legend 'Closed loop';
                xlabel('Time [sec]');ylabel('Residual force [V]');
                figure_specific
                if SW_SAVE_DATA
                    hgsave(h,['level1_time_trace_residule_center_freq_',...
                        num2str(center_freq(ii)),...
                        'Hz_subplot_compare'],'-v6')
                end
            end
            %% parameter convergence
            if jj == 2
                % estimated parameters -2*cos(w*Ts)
                lb1vector   = theta_hat.signals.values;
                % change to frequencies in rad/s
                % (for result evaluation, not needed in controller design)
                w1hat       = abs(acos(-lb1vector/2)/Ts);
                
                figure;
                plot(theta_hat.time,theta_hat.signals.values);
                xlabel('time [sec]');
                ylabel('Estimated parameters')
                grid on;
                if SW_SAVE_DATA
                    hgsave(['level1_para_converge_step_change_center_',...
                        num2str(center_freq(ii)),...
                        'Hz'],'-v6')
                end
                
                % parameter convergence (frequency perspective)
                figure;
                plot(theta_hat.time,w1hat/2/pi,'r');
                xlabel('time (sec)');
                ylabel('Estimated frequency (Hz)');
                grid on;
                if SW_SAVE_DATA
                    hgsave(['level1_freq_converge_step_change_center_',...
                        num2str(center_freq(ii)),...
                        'Hz'],'-v6')
                end
            end
        end
        pause(10); % let the CPU take a 10-sec rest
    end
    level1_time_dom_result_step_change_freq.freq_table =...
        freq_table;
    disp ('===================================================================')
    disp ('test results saved to: level1_time_dom_result_step_change_freq')
    disp ('raw data saved to:     data_step_freq')
    if SW_SAVE_DATA
        try
            movefile('*.fig',['level1_test_result_',date])
        catch
        end
        save(['level1_test_result_',date,'\level1_time_dom_result_step_change_freq'], 'level1_time_dom_result_step_change_freq');
        save (['level1_test_result_',date,'\data_step_freq'],...
            'data_step_freq');
    end
    
    %% TEST FOR CHIRP DISTURBANCE
elseif FLAG_DIST_FREQ == FLAG_CHIRP_DIST
    SW_CHIRP_DIST = 1;
    SW_DIST_ON    = 0;
    t_sim         = chirp_dist.chirp_init_time+chirp_dist.chirp_dur_time*2+chirp_dist.freq2_dur_time*2;
    t_NBon      = 5;
    t_Qon       = t_NBon;
    t_AdapOn    = t_NBon;
    t_AdapOff   = t_sim;
    t_dur_lastDist = 5;
    if ~SW_UNIFORM_ADAP
        adap_init.SW_lambda     = 3;
        adap_init.SW_2Stage     = 0;
        
        alpha_init  = 0.66;
        alpha_end   = 0.865;
        alpha       = alpha_end;
        
        adap_init.alpha_init  = alpha_init;
        adap_init.alpha_end   = alpha_end;
    end
    % for consistency
    t_NBdur     = 0;
    
    dist_seq1   = [50; 0; 0];
    dist_seq2   = [95; 0; 0];
    dist_seq3   = [75; 0; 0];
    
    data_chirp_freq.readme = 'stores the data for the case of chirp changing disturbance frequencies';
    data_chirp_freq.y{1,1} = 'openLoop';
    data_chirp_freq.y{1,2} = 'closedLoop';
    
    ITER_STEP = 1;
    for ii = 1:ITER_STEP
        chirp_dist.para1 = [chirp_dist.freq1_seq(ii); 0;0];
        chirp_dist.para2 = [chirp_dist.freq2_seq(ii); 0;0];
        for jj = 1:2
            if jj == 1
                SW_COMP_ON = 0; % compensation off
                SW_CLOSE_LOOP = 0;
                fig_name_residule_time_trace = ...
                    ['level1_time_trace_residule_chirp_dist_',...
                    num2str(chirp_dist.freq1_seq(ii)),'To',...
                    num2str(chirp_dist.freq2_seq(ii)),...
                    'Hz_openLoop'];
            else
                SW_COMP_ON = 1; % compensation on
                SW_CLOSE_LOOP = 1;
                fig_name_residule_time_trace = ...
                    ['level1_time_trace_residule_chirp_dist_',...
                    num2str(chirp_dist.freq1_seq(ii)),'To',...
                    num2str(chirp_dist.freq2_seq(ii)),...
                    'Hz_closedLoop'];
            end
            %% open simulinnk
            sim(simuName)
            
            data_chirp_freq.initial_freq(ii,:) = chirp_dist.para1';
            data_chirp_freq.y{ii+1,jj} = y;
            if ii == 1 && jj == 1
                if SW_SAVE_DATA
                    if ~exist(['level1_test_result_',date],'dir')
                        mkdir(['level1_test_result_',date]);
                    end
                end
            end
            %% time domain result
            if jj == 2
                level1_time_dom_result_chirp_freq.readme =...
                    {'col 1: init chirp freq;',...
                    'col 2: end chirp freq;',...
                    'col 3: chirp increase freq',...
                    'col 4: chirp decrease freq'};
                
                level1_time_dom_result_chirp_freq.transient_norm(ii,1) = ...
                    chirp_dist.freq1_seq(ii);
                level1_time_dom_result_chirp_freq.transient_norm(ii,2) = ...
                    chirp_dist.freq2_seq(ii);
                level1_time_dom_result_chirp_freq.transient_norm(ii,3) = ...
                    sqrt(...
                    sum(...
                    y.signals.values(...
                    chirp_dist.chirp_init_time/Ts+1 :...
                    (chirp_dist.chirp_init_time + chirp_dist.chirp_dur_time)/Ts ).^2)...
                    );
                level1_time_dom_result_chirp_freq.transient_norm(ii,4) = ...
                    sqrt(...
                    sum(...
                    y.signals.values(...
                    (chirp_dist.chirp_init_time + chirp_dist.chirp_dur_time + 5)/Ts+1 :...
                    (chirp_dist.chirp_init_time + chirp_dist.chirp_dur_time + 10)/Ts ).^2)...
                    );
                
                level1_time_dom_result_chirp_freq.transient_norm_square(ii,1) = ...
                    chirp_dist.freq1_seq(ii);
                level1_time_dom_result_chirp_freq.transient_norm_square(ii,2) = ...
                    chirp_dist.freq2_seq(ii);
                level1_time_dom_result_chirp_freq.transient_norm_square(ii,3) = ...
                    level1_time_dom_result_chirp_freq.transient_norm(ii,3)^2;
                level1_time_dom_result_chirp_freq.transient_norm_square(ii,4) = ...
                    level1_time_dom_result_chirp_freq.transient_norm(ii,4)^2;
                
                level1_time_dom_result_chirp_freq.max_residule(ii,1) = ...
                    chirp_dist.freq1_seq(ii);
                level1_time_dom_result_chirp_freq.max_residule(ii,2) = ...
                    chirp_dist.freq2_seq(ii);
                level1_time_dom_result_chirp_freq.max_residule(ii,3) =...
                    max(y.signals.values(...
                    chirp_dist.chirp_init_time/Ts :...
                    (chirp_dist.chirp_init_time + chirp_dist.chirp_dur_time)/Ts ));
                level1_time_dom_result_chirp_freq.max_residule(ii,4) =...
                    max(y.signals.values(...
                    (chirp_dist.chirp_init_time + chirp_dist.chirp_dur_time + 5)/Ts :...
                    (chirp_dist.chirp_init_time + chirp_dist.chirp_dur_time + 10)/Ts ));
                h = figure;
                plot(y.time,y.signals.values);grid;
                xlabel('Time [sec]');ylabel('Residual force [V]');
                figure_specific
                if SW_SAVE_DATA
                    hgsave(h,fig_name_residule_time_trace,'-v6')
                end
            end
            
            h = figure(FIG_NUMBER_CHIRP_DIST(ii));
            grid on;hold on;
            if jj == 1
                plot(y.time,y.signals.values,'r');
            else
                plot(y.time,y.signals.values,'k:');
                
                legend('open loop','closed loop')
                
                xlabel('Time [sec]');ylabel('Residual force [V]');
                figure_specific
                if SW_SAVE_DATA
                    hgsave(h,['level1_time_trace_residule_chirp_dist_',...
                        num2str(chirp_dist.freq1_seq(ii)),'To',...
                        num2str(chirp_dist.freq2_seq(ii)),...
                        'Hz_compare'],'-v6')
                end
            end
            
            h = figure(FIG_NUMBER2_CHIRP_DIST(ii));
            hold on;
            if jj == 1
                subplot(211)
                plot(y.time,y.signals.values,'r');
                
                legend 'Open loop';
                
                ylabel('Residual force [V]');
                figure_specific
            else
                subplot(212)
                plot(y.time,y.signals.values,'k');
                
                legend 'Closed loop';
                
                xlabel('Time [sec]');ylabel('Residual force [V]');
                figure_specific
                if SW_SAVE_DATA
                    hgsave(h,['level1_time_trace_residule_chirp_dist_',...
                        num2str(chirp_dist.freq1_seq(ii)),'To',...
                        num2str(chirp_dist.freq2_seq(ii)),...
                        'Hz_subplot_compare'],'-v6')
                end
            end
            %% parameter convergence
            if jj == 2
                lb1vector = theta_hat.signals.values;
                % Frequency in rad/s
                w1hat = abs(acos(-lb1vector/2)/Ts);
                % parameter convergence
                figure;
                plot(theta_hat.time,theta_hat.signals.values);
                xlabel('time [sec]');
                ylabel('Estimated parameters')
                grid on;
                xlim([0,t_sim])
                if SW_SAVE_DATA
                    hgsave(['level1_para_converge_chirp_dist_',...
                        num2str(chirp_dist.freq1_seq(ii)),'To',...
                        num2str(chirp_dist.freq2_seq(ii)),...
                        'Hz'],'-v6')
                end
                
                % para convergence (frequency perspective)
                figure;
                plot(theta_hat.time,w1hat/2/pi,'r');
                xlabel('time (sec)');
                ylabel('Estimated frequency (Hz)');
                grid on;
                xlim([0,t_sim])
                if SW_SAVE_DATA
                    hgsave(['level1_freq_converge_chirp_dist_',...
                        num2str(chirp_dist.freq1_seq(ii)),'To',...
                        num2str(chirp_dist.freq2_seq(ii)),...
                        'Hz'],'-v6')
                end
            end
        end
        pause(10); % let the CPU take a 10-sec rest
    end
    disp ('===================================================================')
    disp ('test results saved to: level1_time_dom_result_chirp_freq')
    disp ('raw data saved to:     data_chirp_freq')
    if SW_SAVE_DATA
        try
            movefile('*.fig',['level1_test_result_',date])
        catch
        end
        save(['level1_test_result_',date,'\level1_time_dom_result_chirp_freq'],...
            'level1_time_dom_result_chirp_freq');
        save (['level1_test_result_',date,'\data_chirp_freq'],...
            'data_chirp_freq');
    end
end
