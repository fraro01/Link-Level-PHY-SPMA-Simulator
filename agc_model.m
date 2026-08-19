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

% takes in input the current time slot analyzed "rxIn" for the current 'k' 
% burst under alaysis inside the main for loop of the RX
function [rxOut, state] = agc_model(rxIn, Fs, state)

    %% PARAMETERS
    target_dBFS = -12; % decibels relative to full scale

    attackTime = 6e-6;
    releaseTime = 36e-6;

    rfGainMin_dB = 0;
    rfGainMax_dB = 60;

    digitalGainMin_dB = 0;
    digitalGainMax_dB = 24;

    adcNumBits = 10; % maximum number of available bits to represent values
    adcFullScale = 1; % full scale, hence range: [-1;+1], for clipping

    updatePeriod = 1e-6;
    samplesPerUpdate = max(1, round(updatePeriod * Fs));

    %% INITIAL STATE
    if isempty(state)
        state.GRF_dB = 0;
        state.GDIG_dB = 0;
    end

    rxOut = zeros(size(rxIn));

    %% ADC PARAMETERS
    % ADC maximum positive integer code.
    adcMaxCode = 2^(adcNumBits-1) - 1; % maximum representable number of 
    % ADC, in order to be symmetric to 0 --> -511,..,0,...,+511 (total 1023
    % not 1024).
    % Quantization step. distance between two consecutive quantized levels
    adcStep = adcFullScale / adcMaxCode;

    %% AGC LOOP
    for n = 1:samplesPerUpdate:length(rxIn) % start:pace:end
        % do not exceed the length of the signal
        idxEnd = min(n + samplesPerUpdate - 1, length(rxIn)); % indices 
        idx = n:idxEnd; % array of sample indices
        x = rxIn(idx); % related samples

        %% Measure input level
        rmsIn = sqrt(mean(abs(x).^2)); % Root Mean Square

        inputLevel_dBFS = 20*log10(max(rmsIn / adcFullScale, 1e-12));

        %% Required total gain
        desiredTotalGain_dB = target_dBFS - inputLevel_dBFS;

        %% RF AGC
        desiredRFgain_dB = desiredTotalGain_dB - state.GDIG_dB;

        if desiredRFgain_dB < state.GRF_dB
            % ATTACK: received level too high,
            % reduce gain quickly
            tau = attackTime;
        else
            % RELEASE: received level too low,
            % increase gain more slowly
            tau = releaseTime;
        end

        alpha = 1 - exp(-updatePeriod/tau);

        state.GRF_dB = state.GRF_dB + ...
            alpha*(desiredRFgain_dB - state.GRF_dB);

        state.GRF_dB = min( max(state.GRF_dB, rfGainMin_dB), rfGainMax_dB);

        %% Apply RF gain
        xRF = x * 10^(state.GRF_dB/20);

        %% ADC clipping
        % Separate I and Q components.
        xr = real(xRF);
        xi = imag(xRF);
        % Input clipping. [-1;+1], very clever with only one row of code
        xr = min(max(xr,-adcFullScale),adcFullScale);
        xi = min(max(xi,-adcFullScale),adcFullScale);

        %% ADC quantization
        % Uniform quantization. Uniform since the quantization step is 
        % equal. how many quantization steps corresponds the input value, 
        % round them to the nearest integer, then go back to the normalized 
        % value
        xr = adcStep * round(xr/adcStep);
        xi = adcStep * round(xi/adcStep);
        % Complex ADC output. after clipping and quantization
        xADC = xr + 1j*xi;

        %% DIGITAL AGC
        rmsADC = sqrt(mean(abs(xADC).^2));

        levelADC_dBFS = 20*log10(max(rmsADC/adcFullScale,1e-12));

        desiredDigitalGain_dB = target_dBFS - levelADC_dBFS;

        state.GDIG_dB = min( ...
            max(desiredDigitalGain_dB, digitalGainMin_dB), ...
            digitalGainMax_dB);

        %% Apply digital gain
        rxOut(idx) = xADC * 10^(state.GDIG_dB/20);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    END AGC RF - ADC - AGC DIG                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%