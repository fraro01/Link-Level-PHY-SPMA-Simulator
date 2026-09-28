% This is a secondary auxiliary file for the simulator, hence there must
% not be any code before the 'function' declaration.

%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
* 
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      AGC RF - ADC - AGC DIG                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AGC Gain adds a variable analog gain added to the fixed 90 dB front-end

% takes in input the current time slot analyzed "rxIn" for the current 'k' 
% burst under alaysis inside the main for loop of the RX
function [rxOut, vettore] = agc_model(rxIn, Fs)
    %% PARAMETERS
    target_dBFS = -12; % decibels relative to full scale (0 would be: 
    % 10*log10(input_value / max_scale_value))
    % RF AGC dynamics
    attackTime = 6e-6; % based on the AGC clock discretization
    releaseTime = 36e-6; % based on the AGC clock discretization
    % Gain limits
    rfGainMin_dB = 0;
    rfGainMax_dB = 60;
    digitalGainMin_dB = 0;
    digitalGainMax_dB = 24;
    % ADC
    adcNumBits = 10; % maximum number of available bits to represent values
    adcFullScale = 1; % full scale, hence range: [-1;+1], for clipping, it 
    % is defined in AMPLITUDE.
    % how often the AGC controller observes the time slot, estimates the 
    % level, calculates the error, updates the gain. it's the 
    % discretization of the AGC clock
    updatePeriod = max(1e-7, 1/Fs); % seconds. AGC update period. NOTE! you should  
    % not go below: 1/Fs s, in the next code line there is the reason why! 
    % Number of samples contained in one AGC update window
    samplesPerUpdate = max(1, round(updatePeriod * Fs));

    % LOCK / GAIN FREEZE
    lockThreshold_dB   = 1;
    unlockThreshold_dB = 2;
    % Number of consecutive AGC updates required
    % before declaring Lock
    numUpdatesForLock = 3;

    %% INITIAL STATE
    % initialize the state
    state.GRF_dB = 0;
    state.GDIG_dB = 0;
    % at first the initial state is not locked
    state.locked = false;
    state.lockCounter = 0;

    %% ADC PARAMETERS
    % ADC maximum positive integer code.
    adcMaxCode = 2^(adcNumBits-1) - 1; % maximum representable number of 
    % ADC, in order to be symmetric to 0 --> -511,..,0,...,+511 (total 1023
    % not 1024).
    % Quantization step. distance between two consecutive quantized levels
    adcStep = adcFullScale / adcMaxCode;
    
    % initialize output vector
    rxOut = complex(zeros(size(rxIn))); % same length as the input (time 
    % slot)   
    
    vettore = [];

    % AGC LOOP
    % NOTE! the update pace is on the basis of the discretization clock set
    for n = 1:samplesPerUpdate:length(rxIn) % start:pace:end
        % do not exceed the length of the signal
        idxEnd = min(n + samplesPerUpdate - 1, length(rxIn)); % indices 
        idx = n:idxEnd; % array of sample indices
        x = rxIn(idx); % related samples

        % TODO! BASED ON THE POWER LEVEL, SET A POSSIBLE "IDLE" STATE?
        % 'continue' MAY EXIT THE FOR CYCLE
        % based on the 'vettore' that you get from the 'main.m', plot it,
        % and then set the proper threshold of the follwing IDLE state
        % accordingly.
        vettore = [vettore; sqrt(mean(abs(x).^2))];
    %{%}
        % it also resembles the logic of the Squelch of the radio
        if sqrt(mean(abs(x).^2)) <= 0.05 % NOTE! this strongly improves EVM 
            rxOut(idx) = 0; % if no transmission is found set the level to zero!
            continue
        end
    
        %% 1. APPLY CURRENT RF GAIN
        % IMPORTANT:
        % state.GRF_dB has been calculated during the
        % PREVIOUS AGC window.
        % Hence the model is causal window by window.
        rfGainLin = 10^(state.GRF_dB/20); % referred to amplitude
        xRF = x * rfGainLin; % gain applied
        
        %% 2.1. ADC CLIPPING
        % Separate I and Q components.
        xr = real(xRF);
        xi = imag(xRF);
        % Input clipping. [-1;+1], very clever with only one row of code
        xr = min(max(xr,-adcFullScale),adcFullScale);
        xi = min(max(xi,-adcFullScale),adcFullScale);
        %% 2.2. ADC QUANTIZATION
        % Uniform quantization. Uniform since the quantization step is 
        % equal. how many quantization steps corresponds the input value, 
        % round them to the nearest integer, then go back to the normalized 
        % value
        xr = adcStep * round(xr/adcStep);
        xi = adcStep * round(xi/adcStep);
        % Complex ADC output. after clipping and quantization
        xADC = xr + 1j*xi;
        
        %% 3. APPLY CURRENT DIGITAL GAIN
        % This gain was also calculated during the
        % previous AGC window.
        digitalGainLin = 10^(state.GDIG_dB/20); % reffered to amplitude
        y = xADC * digitalGainLin;
        rxOut(idx) = y; % OUTPUT

        %% 4. MEASURE INPUT LEVEL
        % This measurement will determine the RF gain
        % to be applied starting from the NEXT window.
        rmsIn = sqrt(mean(abs(x).^2)); % Root Mean Square of the AMPLITUDE
        % conversion in decibel relative to Full Scale (referred to amplit)
        inputLevel_dBFS = 20*log10(max(rmsIn/adcFullScale, 1e-12)); % max,
        % protection is set
        % DESIRED RF GAIN, distance
        desiredRFgain_dB = target_dBFS - inputLevel_dBFS;

        %% 5. CURRENT OUTPUT ERROR
        % Used for LOCK / GAIN FREEZE.
        rmsOut = sqrt(mean(abs(y).^2)); % Root Mean Square of the AMPLITUDE
        % conversion in decibel relative to Full Scale (referred to amplit)
        outputLevel_dBFS = 20*log10(max(rmsOut/adcFullScale,1e-12)); % max,
        % protection is set
        % error at the current iteration, calculated as the distance
        % between what we have set now from our target
        outputError_dB = outputLevel_dBFS - target_dBFS;

% Hysteresis is a phenomenon whereby the state of a system depends not only 
% on the current values ​​of the applied forces, but also on its past history 
% The term derives from the Greek hystéresis (delay), indicating the 
% tendency  of a system to react with a delay or to retain a memory of the  
% stresses it has experienced.

        %% 6. LOCK / GAIN FREEZE
        if state.locked
            % -----------------------------------------
            % HYSTERESIS
            % Once locked, we remain locked until the
            % error becomes larger than the wider
            % unlock threshold.
            % -----------------------------------------
            if abs(outputError_dB) > unlockThreshold_dB
                state.locked = false;
                state.lockCounter = 0; % initialize the lock counter
            end
        % if we are not currenlty locked
        else
            % -----------------------------------------
            % LOCK PERSISTENCE
            % One single good measurement is not enough.
            % The error must remain inside the Lock
            % region for several consecutive updates.
            % -----------------------------------------
            if abs(outputError_dB) <= lockThreshold_dB
                % increment the lock counter
                state.lockCounter = state.lockCounter + 1;
            % if the error is larger than the lockThreshold_dB REinitialize
            % the lockCounter from 0
            else
                state.lockCounter = 0;
            end
            % if the resulting lockCounter is greater than the "lock"
            % threshold, we enter the lock state
            if state.lockCounter >= numUpdatesForLock
                state.locked = true;
            end
        end

        %% 7. CALCULATE GAINS FOR NEXT WINDOW
        % if we are not currently locked
        if ~state.locked
            %% 7.1. RF AGC ATTACK / RELEASE
            if desiredRFgain_dB < state.GRF_dB
                % ATTACK
                % Input is too strong.
                % Reduce RF gain quickly.
                % The currently required RF gain is lower than the gain
                % presently applied -> decrease RF gain quickly.
                tau = attackTime;
            else
                % RELEASE
                % Input is too weak.
                % Increase RF gain more slowly.
                tau = releaseTime;
            end
            % First-order exponential dynamics
            % discretization of an exponential response of first order. 
            % Alpha answers: what fraction of the distance we travel per 
            % each update. 
            % e.g.: if updatePeriod = 1 μs; tau = 6 μs ⇒ α = 0.15, hence 
            % per each 'updatePeriod' μs we travel 15% of the residual 
            % distance from the desired gain.
            % 'attackTime' is the time constant τ of the exponential, after 
            % ONE time constant τ, we travel: (1 − e^−1) ≈ 63.2%; hence  
            % after: 'updatePeriod' = τ µs, we are at the 63.2% of the 
            % correction
            alpha = 1 - exp(-updatePeriod/tau);
            % RF gain for NEXT AGC window
            % gain converges gradually with alpha, acting as a transitory
            newGRF_dB = state.GRF_dB + alpha * (desiredRFgain_dB - ...
                state.GRF_dB);
            % RF gain limits
            % update the state at current pace, making it within boundaries
            newGRF_dB = min(max(newGRF_dB,rfGainMin_dB), rfGainMax_dB);

            %% 8 DIGITAL AGC
            % Measure ADC output level
            rmsADC = sqrt(mean(abs(xADC).^2)); % referred to amplitude
            levelADC_dBFS = 20*log10( max(rmsADC/adcFullScale,1e-12));
            % Residual correction required
            % after RF AGC + ADC
            newGDIG_dB = target_dBFS - levelADC_dBFS;
            % Digital gain limits
            newGDIG_dB = min(max(newGDIG_dB, digitalGainMin_dB), ...
                digitalGainMax_dB);

            %% SAVE NEW STATE
            % These gains will be applied starting from
            % the NEXT AGC update window.
            state.GRF_dB = newGRF_dB;
            state.GDIG_dB = newGDIG_dB;

        end

        % IMPLICIT GAIN FREEZE
        % If state.locked == true:   % (that is a condition verified above)
        %
        % GRF and GDIG are NOT updated.
        %
        % This is the Gain Freeze.

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    END AGC RF - ADC - AGC DIG                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%