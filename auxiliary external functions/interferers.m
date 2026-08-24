% This is a secondary auxiliary file for the simulator, hence there must 
% not be any code before the 'function' declaration.

%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
* 
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           INTERFERERS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function interferer_waveform = interferers(numInterferers, txSCFDEBursts)
    
    % initialize the final interferer_waveform
    interferer_waveform = 0;

    for i = 1:numInterferers
    
        %% Originate the Packet
        % retrieve the proper size of the final transmitted bursts inside the final
        % transmitted signal
        numBursts = size(txSCFDEBursts, 2);
        bitsPerBurst = 3 * size(txSCFDEBursts, 1);
        txBursts = randi([0 1], bitsPerBurst, numBursts);
        
        %% 8-PSK Modulation
        M = 8;
        bitsPerSymbol = log2(M); % number of bits per symbol transmitted
        modulatedBursts = zeros(size(txSCFDEBursts, 1), numBursts); % matrix
        for k = 1:numBursts
            bits_k = txBursts(:,k); % k column of the 8-psk padded burst
            % matrix with #bitsPerSymbol rows, and columns automatically set by [],
            % we then transpose with .' hence each row contains one symbol
            bitGroups = reshape(bits_k, bitsPerSymbol, []).';
            % converts each row from binary to decimal (pskmod requires decimal)
            symbolIndices = bi2de(bitGroups, 'left-msb');
            % assigns to the k column the result of the modulation
            modulatedBursts(:,k) = pskmod( ...
                symbolIndices, ...
                M, ...
                pi/M, ... % initial phase of the constellation (shifted by pi/M) for pahse and synchronization purposes in receiver, the decisional boundaries coincide with neat angles
                'gray'); % Gray coding, hence close symbols differ for only 1 bit
        end
        
        % NOTE! from now on should be the same as the TRANSMITTER block of the main
        % simulator! excpet that now we have switched the name: "txSCFDEBursts" 
        % with "modulatedBursts"
        
        %% Time-Frequency Hopping
        % NOTE! we are also distribuiting the burst over a pseudo-casual time
        % pattern.
        % symbol rate = 13 MHz; maximum 17 MHz, check if the one derived is
        % compliant to these numbers.
        % parameters
        numFrequencyChannels = 4;
        symbolsPerBurstWithCP = size(modulatedBursts,1);
        dutyCycle = 0.25;
        fMin = 1435e6;
        fMax = 1518e6;
        % Configuration of Root Raised Cosine TX Filter, we configure the "rolloff"
        % now, because is used to calculate the Symbol Rate from the formula below
        rolloff = 0.35; % how much extra band we use for the RRC filter
        % TLC formula for Symbol Rate: symbolRate = channelBandwidth/(1+rolloff)
        symbolRate = ((fMax - fMin)/numFrequencyChannels)/(1 + rolloff);
        % timing
        burstTxTime = symbolsPerBurstWithCP / symbolRate;     % seconds
        slotDuration = burstTxTime / dutyCycle;               % seconds
        % fundamental check of the logic: burstTxTime / slotDuration = dutyCycle
        % Frequency channels
        channelBandwidth = (fMax - fMin) / numFrequencyChannels;
        carrierFrequencies = fMin + channelBandwidth/2 + ...
            (0:numFrequencyChannels-1)*channelBandwidth; % list comprehension (0:i)
        % Time-frequency hopping pattern
        rng();   % if you want results reproducibility insert a number
        timeSlots = 1:numBursts; % vector column
        % Random offset inside each slot, without taking into account RRC tails in
        % which a burst could be transmitted even in the final edge of its timeSlot
        % maxOffset = slotDuration - burstTxTime;
        % The below definition of maxOffset better cope with the tails introduced 
        % by the Pulse Shaping Filters.
        filterSpan = 12; % RRC filter duration in symbols NOTE! the longer it lasts
        % over time, the sharper it is over the frequencies. 
        % I have defined it here, even though the RRC is implemented in the next
        % block. But it is essential to define it here in order to prevent the
        % transmission of a burst to create tails that exceed their own allocated 
        % time slot.
        shapedBurstTime = (symbolsPerBurstWithCP + filterSpan) / symbolRate; 
        % TX tuner delay: 1 microsecond
        tunerDelay_s = 1e-6;
        % Maximum random intra-slot offset.
        % Reserve tunerDelay_s at the beginning so that the complete
        % transmitted waveform remains inside its allocated time slot.
        % the tuner delay is absorbed inside the allocated time slot of the burst!!
        maxOffset = slotDuration - shapedBurstTime - tunerDelay_s;
        % Random transmission position inside the available part of the slot
        intraSlotOffset = rand(numBursts,1) * maxOffset; % uniform distribution!
        % array of values between 1 and #numFrequencyChannels of length #numBursts
        frequencyChannels = randi( ...
            [1 numFrequencyChannels], ...
            numBursts, ...
            1);
        carrierForBurst = carrierFrequencies(frequencyChannels).';
        startTimes = (timeSlots - 1).' * slotDuration + intraSlotOffset + ...
                    tunerDelay_s; % arrray
        endTimes = startTimes + burstTxTime; % array

        % ASSUMPTION!
        % startTimes represents the burst arrival time at the simulated RX,
        % not the physical transmission time at the interfering transmitter.
        % Therefore, the absolute propagation delay distance_i/c is not
        % explicitly modeled.
        % distance_i is still used for FSPL and channel geometry.
        % The Rician PathDelays represent only multipath excess delays
        % relative to the first arriving path.
        
        %% Final transmitted EQUIVALENT-BASEBAND signal with Pulse Shape Filter RRC
        % NOTE! Due to Nyquist and thanks to the fact that we are simulating in
        % baseband, we must use: samplesPerSymbol > 2*(fmax_occupied)/symbolRate
        % in which fmax_occupied = [max(relativeCarrierFrequencies) +
        % (relativeCarrierFrequencies(4) - relativeCarrierFrequencies(3))]
        % Atually, with RRC the one-sided band of a single channel is: 
        % Bsingle = symbolRate/2 * (1 + rolloff) --> hence fmax_occupied is less.
        % NOTE! see: https://en.wikipedia.org/wiki/Raised-cosine_filter for time
        % and frequency domain plots, for a better understanding of RRC. RRC makes
        % symbols extremely overlap over the time, but thanks to the first Nyquist
        % criterion for pulse shaping filters, exactly in the sampling points the
        % interferers are equal to zero.
        samplesPerSymbol = 8; % number of samples per symbol
        % sample rate --> how much finely we represent the waveform over the time.
        % symbol rate --> how quick the info changes.
        Fs = symbolRate * samplesPerSymbol; % sampling frequency, since we are 
        % constructing the signal from the symbol to the the continuous time domain 
        % we need more samples to represent one symbol.
        % NOTE! each RRC burst has null-to-null bandwidth: symbolRate*(1+rolloff)
        % IF we transmit simultaneously multiple bursts per time slot we could
        % experience ACI, since our channel spacing is roughly 13.83 MHz
        % filterSpan = 12; filter duration in symbols, NOTE! the longer it lasts
        % over time, the sharper it is over the frequencies. I have commented it
        % since it has been already defined above for the TF-Hopping, so that a
        % burst transmission cannot happen at the very end of its time slot
        txRRC = comm.RaisedCosineTransmitFilter( ...
            'Shape','Square root', ...
            'RolloffFactor',rolloff, ...
            'FilterSpanInSymbols',filterSpan, ...
            'OutputSamplesPerSymbol',samplesPerSymbol); % how many samples it 
        % outputs per symbol the RRC
        totalTime = max(endTimes); % final time for transmitting the whole packet 
        % (last burst)
        extraSymbols = filterSpan; % Extra zero symbols appended to flush the TX/RX 
        % pulse-shaping filters.
        tailSamples = filterSpan * samplesPerSymbol; % tail extention in samples 
        numSamplesTotal = ceil(totalTime * Fs) + tailSamples + 1; % +1 for safety 
        % reasons
        txSignal = zeros(numSamplesTotal,1);
        % we move to baseband, to avoid computational overheads of L-band (Nyquist)
        fCenter = mean(carrierFrequencies);
        % Frequencies relative to bandwidth center
        relativeCarrierFrequencies = carrierFrequencies - fCenter;
        for k = 1:numBursts
            reset(txRRC); % clean the filter memory per each burst
            burstSymbols = modulatedBursts(:,k); % from where we sample
            % Add zeros to flush the TX/RX filter memory response tail for last
            % symbol
            burstSymbolsPadded = [burstSymbols; zeros(extraSymbols,1)];
            % Pulse shaping
            shapedBurst = txRRC(burstSymbolsPadded);
            % startTimes(k) indicates the beginning of the shaped TX waveform,
            % not the sampling instant of the first useful symbol.
            startSample = round(startTimes(k)*Fs) + 1;  % initial sampling window
            % ending sampling window
            endSample = startSample + length(shapedBurst) - 1;
            % absolute temporal axes
            n = (startSample-1:endSample-1).';
            t = n/Fs; % in order to shift a signal we must multiply it by a complex 
            % time variant sinusoid
            fRel = relativeCarrierFrequencies(frequencyChannels(k));
            % originally the burst was center in 0 Hz, now we center it to fRel Hz
            % multiplying a complex sinusoid varying in time, is equal to shift its 
            % signal in the frequency
            burstUpconverted = shapedBurst .* exp(1j*2*pi*fRel*t);
            % populating the signal to send
            txSignal(startSample:endSample) = ...
                txSignal(startSample:endSample) + burstUpconverted;
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   CHANNEL GENERATION PER INTERFER                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1. Topology
    % x = a + (b-a)*rand(); generates a a number picked uniformly between:
    % [a; b]
    distance_i = 0 + (500000 - 0)*rand(); % [m]
    velocity_i = 0 + (1372 - 0)*rand();   % [m/s]
    angle_i = randi([0 180]);             % degree
    %% 2. Power and Generic Gain Settings
    % transmitter power
    txPower_dBm = randi([40 48]); % [dBm]
    propagationGeneralGain = 0 + (13 - 0)*rand(); % [dB]
    %% 3. Normalize the Interferer waveform to the selected TX Interf power
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
    %% 4. Free Space Path Loss and Link Budget
    speedOfLight_mps = physconst('LightSpeed'); % for the wavelength calculus
    % Sample-wise propagation attenuation.
    largeScaleAmplitude = zeros(size(txSignalPhysical));
    % one value for each burst.
    fsplPerBurst_dB = zeros(numBursts,1);
    % different bursts are on different carriers leading to different FSPL
    for k = 1:numBursts %iterate over all the N bursts
        % Actual RF carrier used by burst k.
        carrierHz = carrierForBurst(k);
        % Wavelength corresponding to the current carrier.
        wavelength_m = speedOfLight_mps / carrierHz;
        % Free Space Path Loss in dB.
        fsplPerBurst_dB(k) = fspl(distance_i, wavelength_m);
        % Propagation gain excluding txPower_dBm.
        % txPower_dBm is excluded because txSignalPhysical has already
        % been normalized to the actual transmitter power.
        % link budget calculation
        propagationGain_dB = propagationGeneralGain - fsplPerBurst_dB(k); 
        % Convert power gain in dB into an amplitude multiplier.
        propagationAmplitudeGain = 10^(propagationGain_dB/20);
        % First sample occupied by the pulse-shaped burst.
        burstStartSample = round(startTimes(k)*Fs) + 1;
        % Include useful waveform and RRC filter tail.
        burstLengthSamples = ...
            (symbolsPerBurstWithCP + extraSymbols) * samplesPerSymbol;
        burstEndSample = min( ...
            burstStartSample + burstLengthSamples - 1, ...
            length(txSignalPhysical));
        % Assign the attenuation to samples belonging to this burst.
        largeScaleAmplitude(burstStartSample:burstEndSample) = ...
            propagationAmplitudeGain;
    end
    % Apply large-scale gains and losses. (expressed as amplitude, not power)
    txSignalAfterLinkBudget = txSignalPhysical .* largeScaleAmplitude;
    %% 5. Doppler Parameters
    % Corresponding wavelength.
    minimumWavelength_m = speedOfLight_mps / fMax;
    centralWavelength_m = speedOfLight_mps / fCenter;
    % Conservative maximum Doppler magnitude used to parametrize the
    % diffuse fading process of all channel paths, for doppler spread spectrum
    maximumDopplerShift_Hz = velocity_i / minimumWavelength_m;
    % Signed Doppler shift of the direct LOS component
    directPathDopplerShift_Hz = velocity_i / centralWavelength_m;
    %% 6. Multipath & Fading
    pathDelays_s = [0, ...
                    1/(symbolRate) + (2/symbolRate - 1/(symbolRate))*rand(), ...
                    2/symbolRate + (3/symbolRate - 2/(symbolRate))*rand()];
    averagePathGains_dB = [0, ... % 0 gain with respect to reference (itself is taken as reference indeed it makes sense zero)
                           randi([-16 -8]), ... % with -12 dB --> the path has 6.31% of the power of the first path, or equivalently a 16th
                           randi([-22 -14])];
    ricianKFactor = [randi([25 35]), ...
                     0, ... % the deterministic compoment LOS power is 0 --> Rayleigh channel
                     0];    
    ricianChannel = comm.RicianChannel( ...
        'SampleRate', Fs, ... % how many samples per second the channel contains
        'PathDelays', pathDelays_s, ... % seconds
        'AveragePathGains', averagePathGains_dB, ... % average power of each path
        'NormalizePathGains', true, ... % normalize the total average power of all the paths of the multipath to 1 (otherwise is greater), see below
        'KFactor', ricianKFactor, ... % ratio between power of LOS component / power casual of diffuse component
        'DirectPathDopplerShift', [directPathDopplerShift_Hz*cosd(angle_i) directPathDopplerShift_Hz*cosd(45) directPathDopplerShift_Hz*cosd(5)], ... % doppler of the LOS component of each Rician path; it is applied only for Rician paths
        'DirectPathInitialPhase', [0 pi/3 -pi/2], ... % initial phase of LOS component of each Rician path; they depend on the distance TX-RX. It represents the phase with which the specific path arrives at the RX when the simulation starts.
        'MaximumDopplerShift', maximumDopplerShift_Hz, ... % maximum Doppler magnitude of the diffuse fading components paths
        'PathGainsOutputPort', false, ... % returns also the instaneous gains of all the paths
        'RandomStream', 'Global stream'); % keeps generating casual fading components
    % Apply time-varying multipath channel.
    interfererTXWaveformFaded = ricianChannel(txSignalAfterLinkBudget);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                  END CHANNEL GENERATION PER INTERFER                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % since 'txSignal' and 'interfer_waveform' might have different lengths
    % we bring them to the same length, with zero-padding
    if isscalar(interferer_waveform) && interferer_waveform == 0
    interferer_waveform = interfererTXWaveformFaded;
    else
        maxLength = max(length(interferer_waveform), ...
                length(interfererTXWaveformFaded));
        interferer_waveform(end+1:maxLength, 1) = 0;
        interfererTXWaveformFaded(end+1:maxLength, 1) = 0;
        interferer_waveform = interferer_waveform + interfererTXWaveformFaded;
    end

    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          END INTERFERERS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%