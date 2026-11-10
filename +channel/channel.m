%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           AVIONIC CHANNEL                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Channel pipeline:
% txSignal
%    -> TX power normalization
%    -> FSPL and Link Budget
%    -> Doppler
%    -> Rician multipath fading
%    -> oscillator CFO and Phase Offset
%    -> thermal AWGN
%    -> RX front-end gain

% The complex-baseband samples are interpreted as voltage amplitudes
% across a 1-ohm reference load.
% For a generic reference load R:
%     instantaneousPower_W = abs(x[n])^2 / R
% Since R = 1 ohm:
%     instantaneousPower_W = abs(x[n])^2
% This convention must be used consistently for both signal scaling (3.)
% and thermal-noise generation (9.)

function [blankingMask, ...
          rxSignal, ... % FINAL SIGNAL IN RX     
          rxSignalImpaired, ...
          rxSignalAtAntenna, ...
          rxSignalWithInterference, ...
          interferer_waveform, ...
          pathGains, ...
          pathDelays_s, ...
          delay_max, ...
          friendBurstGroundTruth, ...
          max_interferer_delay ...
         ] = channel(...
          txSignal, ... % signal in input
          numBursts, ...
          fMin, ...
          fMax, ...
          fCenter, ...
          carrierForBurst, ... % carrier frequencies per burst
          symbolRate, ...
          samplesPerSymbol, ...
          Fs, ...
          startTimes, ...
          symbolsPerBurstWithCP, ...
          filterSpan, ...
          rolloff, ...
          numFrequencyChannels, ...
          numInterferers, ...
          blankingPercentage, ...
          tunerDelay_s, ...
          dutyCycle, ...
          burstTxTime, ...
          syncReference, ...
          phaseReference, ...
          bitsPerBurst, ...
          numProtectedIndexBits, ...
          friendIndexCodeword, ...
          cpLength ...
          )

    %% 1. Avionic scenario parameters
    % ASSUMPTION!
    % Distance and velocities remain approximately constant during one packet.
    % we could implement a sort of time variance for future implementations.
    linkDistance_m = 50e3;          % TX-RX distance
    aircraftSpeed_mps = 686;        % aircraft relative speed magnitude [m/s]
    % Convention used here:
    % radialVelocity_mps > 0  --> approaching
    % radialVelocity_mps < 0  --> moving away
    
    %% 2. Link-budget parameters
    % Transmitter power
    txPower_dBm = 48; % [dBm]                formula: x dBm = 10^(x/10) mW; 
    % dBm (Decibels RELATIVE to 1 milliwatt)
    % Antenna gains
    txAntennaGain_dBi = 6; % [dBi] (Decibels RELATIVE to Isotropic antenna)
    rxAntennaGain_dBi = 6; % RX antenna gain [dBi]
    % Feeder, cable and connector losses
    txFeederLoss_dB = 0.5;  % TX cable/connector loss [dB]
    rxFeederLoss_dB = 0.5;  % RX cable/connector loss [dB]
    % Additional propagation and implementation losses
    otherPropagationLosses_dB = 1; % polarization, implementation, etc. [dB]
    % Receiver noise and front-end parameters:
    receiverNoiseFigure_dB = 2; % added by the receiver device, degrading SNR.
        
    %% 3. Normalize the transmitted waveform to the selected TX power
    % txSignal currently has arbitrary normalized MATLAB units.
    % Select only active transmitted samples (ommit empty tx intervals).
    txActiveSamples = abs(txSignal) > 1e-12;
    if ~any(txActiveSamples)
        error('The transmitted waveform contains no active samples.');
    end
    % Average power of the current normalized MATLAB waveform.
    txSignalPowerArbitrary = mean(abs(txSignal(txActiveSamples)).^2);
    % Convert transmitter power from dBm to watts. 
    txPower_W = 10^((txPower_dBm - 30)/10); % standard formula
    % FUNDAMENTAL SCALING: |txSignalPhysical|^2 is interpreted as instantaneous
    % power in watts, and its active-sample average is txPower_W.
    txSignalPhysical = txSignal .* sqrt(txPower_W / txSignalPowerArbitrary);
    % square root is necessary since we are modifying the amplitude
    
    %% 4. Free Space Path Loss and link budget
    speedOfLight_mps = physconst('LightSpeed'); % for the wavelength calculus
    % Sample-wise propagation attenuation.
    largeScaleAmplitude = zeros(size(txSignalPhysical));
    % one value for each burst.
    fsplPerBurst_dB = zeros(numBursts,1);
    % for the calculating the desired: 'rxFrontEndGain_dB'
    propagationGainPerBurst_dB = zeros(numBursts,1); 
    % different bursts are on different carriers leading to different FSPL
    for k = 1:numBursts %iterate over all the N bursts
        % Actual RF carrier used by burst k.
        carrierHz = carrierForBurst(k);
        % Wavelength corresponding to the current carrier.
        wavelength_m = speedOfLight_mps / carrierHz;
        % Free Space Path Loss in dB.
        fsplPerBurst_dB(k) = fspl(linkDistance_m, wavelength_m);
        % Propagation gain excluding txPower_dBm.
        % txPower_dBm is excluded because txSignalPhysical has already
        % been normalized to the actual transmitter power.
        % link budget calculation
        propagationGain_dB = ...
                txAntennaGain_dBi ...
                + rxAntennaGain_dBi ...
                - txFeederLoss_dB ...
                - rxFeederLoss_dB ...
                - otherPropagationLosses_dB ...
                - fsplPerBurst_dB(k); 
        % store it
        propagationGainPerBurst_dB(k) = propagationGain_dB;
        % Convert power gain in dB into an amplitude multiplier.
        propagationAmplitudeGain = 10^(propagationGain_dB/20);
        % First sample occupied by the pulse-shaped burst.
        burstStartSample = round(startTimes(k)*Fs) + 1;
        % Include useful waveform and RRC filter tail.
        burstLengthSamples = ...
            (symbolsPerBurstWithCP + filterSpan) * samplesPerSymbol;
        burstEndSample = min( ...
            burstStartSample + burstLengthSamples - 1, ...
            length(txSignalPhysical));
        % Assign the attenuation to samples belonging to this burst.
        largeScaleAmplitude(burstStartSample:burstEndSample) = ...
            propagationAmplitudeGain;
    end
    % Apply large-scale gains and losses. (expressed as amplitude, not power)
    txSignalAfterLinkBudget = txSignalPhysical .* largeScaleAmplitude;


    %% 5. Multipath parameters
    % for a Rice model:
    %   - If K = 0, the dominant path disappears, and the Rician model reduces 
    %     perfectly to a Rayleigh model.
    %   - As K → ∞, the signal resembles an uncorrupted, stable signal with 
    %     almost zero fading.
    % PathDelays values are expressed in seconds, using symbolRate.
    % Here:
    %   path 1 -> LOS, zero delay ASSUMPTION!
    % The absolute propagation delay is neglected in the proposed channel 
    % model. As the simulator targets link-level physical-layer performance, 
    % the receiver is assumed to be synchronized to the first arriving path. 
    % Consequently, only the excess delays between the multipath components are 
    % modeled, since they determine the channel impulse response and directly 
    % affect equalization and detection performance (relative delays).
    %   path 2 -> reflection delayed by N symbol;
    %   path 3 -> reflection delayed by N+X symbols.
    % delays are taken relative from the symbol rate and samples per symbol and
    % more importantly: are strictly related to the "Lest" parameter of the
    % receiver
    pathDelays_s = [0, ...
                    1/symbolRate, ...
                    2/symbolRate];
    delay_max = max(pathDelays_s); % for the Lest parameter in 'receiver.m'
    % Average relative power of the paths relative to the strongest path, 
    % LOS/echo paths [dB]
    averagePathGains_dB = [0, ... % 0 gain with respect to reference (itself is taken as reference indeed it makes sense zero)
                          -12, ... % with -12 dB --> the path has 6.31% of the power of the first path, or equivalently a 16th
                          -18]; % with -18 dB --> the path has 1.58% of the power of the first path
    % First path is Rician with K>0 (component specular/LOS is dominant), other
    % paths are Rayleigh (because KFactor is supplied as a vector), with a 
    % diffuse component without LOS. 
    % here we've set "KFactor" as a linear ratio, not in dB, e.g.: linear_K=X
    % --> the deterministic component is X times stronger (power) than the 
    % diffuse casual component.
    % the K-factor Rician determines the ratio between: power of the
    % deterministic component (LOS) and the diffuse casual component. Since a
    % Rician channel contains a deterministic component and a diffuse casual
    % component
    ricianKFactor = [30, ...
                     0, ... % the deterministic compoment LOS power is 0 --> Rayleigh channel
                     0]; % the deterministic compoment LOS power is 0 --> Rayleigh channel
    % because a Rician channel with K --> 0 approaches a Rayleigh channel
    % NOTE! actually the others two minor taps are reflected by the sea and
    % probably by the wing of the aircraft --> we could model them with a
    % predominant deterministic component with respect to the diffuse component
    % as well! --> Rice instead of Rayleigh (K>0).
    % Notice that the doppler shift is usually applied only for the Rice,
    % because the Rayleigh is constituted only by many diffuse components by
    % which attributing a single neat frequency shift is too difficult, hence
    % in this case the doppler is only accounted as doppler spread spectrum.
    
    %% 6. Doppler parameters
    % Corresponding wavelength.
    minimumWavelength_m = speedOfLight_mps / fMax;
    centralWavelength_m = speedOfLight_mps / fCenter;
    % Conservative maximum Doppler magnitude used to parametrize the
    % diffuse fading process of all channel paths, for doppler spread spectrum
    maximumDopplerShift_Hz = aircraftSpeed_mps / minimumWavelength_m;
    % Signed Doppler shift of the direct LOS component
    directPathDopplerShift_Hz = aircraftSpeed_mps / centralWavelength_m;
    
    %% 7. Rician multipath fading channel
    ricianChannel = comm.RicianChannel( ...
        'SampleRate', Fs, ... % how many samples per second the channel contains
        'PathDelays', pathDelays_s, ... % seconds
        'AveragePathGains', averagePathGains_dB, ... % average power of each path
        'NormalizePathGains', true, ... % normalize the total average power of all the paths of the multipath to 1 (otherwise is greater), see below
        'KFactor', ricianKFactor, ... % ratio between power of LOS component / power casual of diffuse component
        'DirectPathDopplerShift', [directPathDopplerShift_Hz*cosd(10) directPathDopplerShift_Hz*cosd(45) directPathDopplerShift_Hz*cosd(5)], ... % doppler of the LOS component of each Rician path; it is applied only for Rician paths
        'DirectPathInitialPhase', [0 pi/3 -pi/2], ... % initial phase of LOS component of each Rician path; they depend on the distance TX-RX. It represents the phase with which the specific path arrives at the RX when the simulation starts.
        'MaximumDopplerShift', maximumDopplerShift_Hz, ... % maximum Doppler magnitude of the diffuse fading components paths
        'PathGainsOutputPort', true, ... % returns also the instaneous gains of all the paths
        'RandomStream', 'Global stream'); % keeps generating casual fading components
    % Apply time-varying multipath channel.
    [rxSignalFaded, pathGains] = ricianChannel(txSignalAfterLinkBudget);
    % rxSignalFaded includes:
    % - LOS component;
    % - echoes;
    % - delays of multipath;
    % - temporal variation of fading;
    % - Doppler LOS;
    % - Doppler spectrum of the diffused components.
    % "pathGains" includes the complex instaneous coefficients of the paths

    % regarding: "'NormalizePathGains', true":
    % normalize the total average power of all the paths of the multipath to 1 
    % (otherwise is greater), hence 0 dB; [1+0,0631+0,01585]≈1,07895 --> 
    % P_1≈0,927;P_2≈0,0585;P_3≈0,0147; it avoids that the multipath introduces 
    % more power with respect to the one already determined by link-budget, 
    % since the channel can temporarily increase or decrease power due to 
    % interferences with paths. We do this because the link-budget already 
    % contains all attenuations, the multipath must not take into account 
    % others arbitrary powers. Hence distribute the total power throught the 
    % paths, but keep the average power equal to the one already determined 
    % by the link-budget. But notice that the istantaneous value is not fixed!
    
    %% 8. Oscillator CFO and phase offset
    % Hardware-related oscillator Carrier Frequency Offset (CFO) and Phase 
    % Impairments.
    % set the following Offset parameters according to the simulation
    % These are distinct from physical propagation Doppler.
    trueCarrierPhaseOffset_rad = deg2rad(1); % the Phase offset introduced
    trueCarrierFrequencyOffset_Hz = 250; % the CFO introduced in Hz
    sampleIndex = (0:length(rxSignalFaded)-1).'; % we start from 0 --> -1
    timeAxis = sampleIndex / Fs; % in seconds
    % formula used below: ϕ(t) = ϕ0 + 2πΔft | CFO makes rotating continuously 
    % the phase, Δf rotations per second; a complete rotation values 2π radians
    % hence: dϕ/dt ​= 2π⋅1000
    carrierImpairmentPhase = trueCarrierPhaseOffset_rad ...
                         + 2*pi * trueCarrierFrequencyOffset_Hz * timeAxis;
    % Impair the signal with the parameters defined above
    rxSignalImpaired = rxSignalFaded .* exp(1j*carrierImpairmentPhase);
    
    %% 9. Colored Noise from Interferers (external auxiliary script)
    % colored noise of the intererfers, comes from the auxiliar external script
    % NOTE! from the script: 'interferers.m' you will notice that the
    % interferers generate a waveform that is exactly compliant to the one
    % generated by: 'transmitter.m', with the difference that they come
    % from different stochastic distances in the space, with different
    % powers, different velocities ecc... --> hence it is nornaml not to
    % properly detect all the bursts sent from all the interferers, due to
    % the power in reception of their bursts received by 'receiver.m'
    [interferer_waveform, ...
    friendBurstGroundTruth, ...
    max_interferer_delay] = interferers( ...
                                numInterferers, ...
                                numBursts, ...
                                fMin, ...
                                fMax, ...
                                numFrequencyChannels, ...
                                rolloff, ...
                                filterSpan, ...
                                symbolRate, ...
                                samplesPerSymbol, ...
                                Fs, ...
                                tunerDelay_s, ...
                                dutyCycle, ...
                                syncReference, ...
                                phaseReference, ...
                                bitsPerBurst, ...
                                numProtectedIndexBits, ...
                                friendIndexCodeword, ...
                                cpLength ...
                                ); % calls the external script
    maxLength = max(length(interferer_waveform), length(rxSignalImpaired));
    rxSignalImpaired(end+1:maxLength, 1) = 0;
    interferer_waveform(end+1:maxLength, 1) = 0; 
    rxSignalWithInterference = rxSignalImpaired + interferer_waveform;
    
    %% 10. Thermal AWGN and receiver noise figure
    % the SNR as a consequence of: transmitted power, distance, antennas gains,
    % losses, bandwidth, temperature, noise figure.
    % Convert receiver noise figure from dB to linear noise factor.
    % Receiver-input noise power, standard formula:
    % Pn = k * T * B * F    [W]
    % where:
    % - kT : noise spectral density;
    % - kTB : total noise power ideal;
    % - F is the linear noise factor, it adds to the ideal total noise power  
    % the noise added in reception by RX device.
    % notice that B is the width of the frequencies window through which the RX
    % let pass the noise, (the greater --> the more noise we let get into RX).
    % here it is calculated by comm.ThermalNoise() taking advantage of its
    % relation with the sampling frequency Fs.
    
    % Physical convention used in this simulation:
    %   ReferenceLoad = 1 ohm
    % Therefore, for every complex-baseband sample x[n]:
    %   instantaneousPower_W = abs(x[n])^2 
    % since P = V^2/R
    % This is coherent with the previous TX normalization, where the active
    % waveform was scaled so that its mean abs(x)^2 equals txPower_W.
    % NoiseMethod = 'Noise figure':
    % the receiver noise figure characterizes the noise introduced by the
    % receiver electronics, assuming the standard reference temperature 290 K.
    % Add290KAntennaNoise = true, adds also: thermal noise received from an
    % equivalent 290 K antenna, hence it takes into account also the EM noise 
    % of the environment; other than the equivalent receiver circuit noise 
    % associated with the noise figure.
    % Consequently, the total input-referred noise power is equivalent to:
    %   Pnoise = k * 290 K * Fs * F 
    % for ReferenceLoad = 1 ohm.
    referenceLoad_Ohm = 1;
    thermalNoiseChannel = comm.ThermalNoise( ...
        'NoiseMethod', 'Noise figure', ...
        'NoiseFigure', receiverNoiseFigure_dB, ... % equivalent noise of the RX electronics
        'SampleRate', Fs, ...
        'ReferenceLoad', referenceLoad_Ohm, ... % how to interpret the numerical amplitude of the signal
        'Add290KAntennaNoise', true); % whether we want to introdue the thermal 
    % noise of the antenna at 290 K other than the one of the RX device.
    % Add antenna and receiver thermal noise at the receiver input.
    rxSignalAtAntenna = thermalNoiseChannel(rxSignalWithInterference);
    % rxSignalAtAntenna = rxSignalWithInterference;
    
    % DIAGNOSTIC check 10*log10(mean(abs(rxSignalAtAntenna).^2))
    %% 11. Receiver Analog Front-End FIXED Gain (Opposed to AGC, but
    % Complementary)
    % Analog gain before the normalized ADC model. LNA/IF gain before ADC [dB]
    % Apply a fixed analog front-end gain to bring the signal in reception
    % to make it compliant to our RX ADC.
    % Convert gain in dB into an amplitude multiplier.
    % so this block performs: LNA, mixer, IF amplifier, baseband amplifier
    % Apply the same gain to signal AND noise. So to bring the signal in a
    % level compatible for my ADC. It is a fixed front-end analog gain.
    % rxSignal = rxSignalAtAntenna * rxFrontEndAmplitudeGain;
    nominalPropagationGain_dB = mean(propagationGainPerBurst_dB);
    nominalRxPower_dBm = txPower_dBm + nominalPropagationGain_dB;
    adcFullScalePower_dBm = 30;   % 'adcFullScale=1', R=1 ohm
    targetFrontEnd_dBFS = -13;    % set: '-12' to avoid work for the AGC!
    % of course based on the: 'target_dBFS' value from: 'agc_model.m'
    targetFrontEndPower_dBm = adcFullScalePower_dBm + targetFrontEnd_dBFS;    
    rxFrontEndGain_dB = targetFrontEndPower_dBm - nominalRxPower_dBm;
    rxFrontEndAmplitudeGain = 10^(rxFrontEndGain_dB/20);
    rxSignal = rxSignalAtAntenna * rxFrontEndAmplitudeGain;
    % DIAGNOSTIC check 10*log10(mean(abs(rxSignal).^2))

    %% 12. Blanking Compensation when the RX is Transmitting
    burstDurationSamples = burstTxTime * Fs; % how many samples per burst
    % call the external auxiliary function for TxRxN Blanking Compensation    
    blankingMask = blanking_mask(rxSignal, ...
                                blankingPercentage, ...
                                burstDurationSamples);
    % Apply local-transmit blanking
    rxSignal = rxSignal .* blankingMask; % element-wise multiplication with 1/0
    
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        END AVIONIC CHANNEL                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%