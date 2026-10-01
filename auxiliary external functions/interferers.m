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

function interferer_waveform = interferers( ...
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
                                friendIndexCodeword, ... % since we have FRIENDLY interferers
                                cpLength ...
                                )
    
    % initialize the final interferer_waveform
    interferer_waveform = 0;

    % NOTE! it is basically as close as possible to the 'transmitter.m'
    % script.
    for i = 1:numInterferers
    
        %% Originate the Packet   
        % Safety checks for the reserved FRIEND index.
        assert(length(friendIndexCodeword) == numProtectedIndexBits, ...
            'friendIndexCodeword must contain numProtectedIndexBits bits.');        
        assert(all(ismember(friendIndexCodeword, [0 1])), ...
            'friendIndexCodeword must be binary.');
        % Column-vector representation, consistent with the transmitter framing.
        friendIndexBits = friendIndexCodeword(:); % convert a row array to column       
        % Generate a random payload for every interferer burst.
        % bitsPerBurst is exactly the same payload length used by transmitter.m.
        interfererPayload = randi([0 1], ...
                                  bitsPerBurst, ...
                                  numBursts);
        %% Burst framing
        % EXACT SAME PHY STRUCTURE AS THE DESIRED TRANSMITTER:
        % sync | phase reference | protected index | payload
        % The ONLY intentional structural difference is that every interferer
        % uses the same reserved FRIEND index codeword as indexing bits.       
        framedBursts = zeros(bitsPerBurst + ...
                             length(syncReference) + ...
                             length(phaseReference) + ...
                             numProtectedIndexBits, ...
                             numBursts);
        % fill the above empty matrix
        for k = 1:numBursts        
            framedBursts(:,k) = [ ... % take all rows of cloumn k vertical concatenation of three vectors
                syncReference;         % 48 bit
                phaseReference;        % 48 bit
                friendIndexBits;       % 15 bit: friend burst identifier
                interfererPayload(:,k) % 144 bit coded, from column k
            ];        
        end
        
        %% 8-PSK Modulation
        M = 8;
        bitsPerSymbol = log2(M);        
        % Same padding rule used by transmitter.m.
        bitsPerFramedBurst = size(framedBursts, 1);        
        numPadBits8PSK = mod(bitsPerSymbol - ...
            mod(bitsPerFramedBurst, bitsPerSymbol), bitsPerSymbol);        
        % Add zero padding after the payload.
        framedBurstsPadded = [ ...
            framedBursts;
            zeros(numPadBits8PSK, numBursts)
        ];        
        bitsPerBurst8PSK = size(framedBurstsPadded, 1);        
        symbolsPerBurst = bitsPerBurst8PSK / bitsPerSymbol;        
        modulatedBursts = zeros(symbolsPerBurst, numBursts);        
        for k = 1:numBursts        
            bits_k = framedBurstsPadded(:,k);        
            bitGroups = reshape(bits_k, bitsPerSymbol, []).';
            symbolIndices = bi2de(bitGroups, 'left-msb');        
            modulatedBursts(:,k) = pskmod( ...
                                        symbolIndices, ...
                                        M, ...
                                        pi/M, ...
                                        'gray');        
        end

        %% SC-FDE Modulation, Cyclic Prefix
        txSCFDEBurstsInterferer = zeros(symbolsPerBurst + cpLength, numBursts);        
        for k = 1:numBursts        
            burstSymbols = modulatedBursts(:,k);        
            % Same cyclic-prefix construction used by transmitter.m.
            cyclicPrefix = burstSymbols(end-cpLength+1:end);        
            txSCFDEBurstsInterferer(:,k) = [ ...
                cyclicPrefix;
                burstSymbols
            ];        
        end
        
        %% Time-Frequency Hopping
        symbolsPerBurstWithCP = size(txSCFDEBurstsInterferer, 1);        
        burstTxTime = (symbolsPerBurstWithCP + filterSpan) / symbolRate;        
        channelBandwidth = (fMax - fMin) / numFrequencyChannels;        
        carrierFrequencies = fMin + channelBandwidth/2 + ...
            (0:numFrequencyChannels-1)*channelBandwidth;
        
        frequencyChannels = randi( ...
            [1 numFrequencyChannels], ...
            numBursts, ...
            1);
        carrierForBurst = carrierFrequencies(frequencyChannels).';
       
        %% Duty Cycle + Tuner Delay + Pseudo-Random Burst Timing
        if dutyCycle == 1
            tunerDelayPerBurst = zeros(numBursts,1); 
            startTimes = zeros(numBursts,1); 
            for k = 2:numBursts
                if frequencyChannels(k) ~= frequencyChannels(k-1)
                    tunerDelayPerBurst(k) = tunerDelay_s;   
                else
                    tunerDelayPerBurst(k) = 0;
                end
                % cumulative!
                startTimes(k) = startTimes(k-1) ...
                              + burstTxTime ...
                              + tunerDelayPerBurst(k);
            end
            randomListeningIntervals = 0;
            
        else % dutyCycle < 1

            % Per-burst tuner delay.
            % tunerDelayPerBurst(k) is the delay required BEFORE burst k if the
            % carrier changes with respect to burst k-1.
            tunerDelayPerBurst = zeros(numBursts,1);   
            for k = 2:numBursts    
                if frequencyChannels(k) ~= frequencyChannels(k-1)
                    tunerDelayPerBurst(k) = tunerDelay_s;
                else
                    tunerDelayPerBurst(k) = 0;
                end    
            end
            
            % Total time during which the node is effectively transmitting.
            totalTxTime = numBursts * burstTxTime;    
            % Total observation time required by the desired duty cycle.
            totalTimeDutyCycle = totalTxTime / dutyCycle;    
            % Total non-transmission time available.
            totalNonTxTime = totalTimeDutyCycle - totalTxTime;    
            % Part of the non-TX time necessarily occupied by tuner switching.
            totalTunerDelay = sum(tunerDelayPerBurst);
            
            % Check that the requested duty cycle is physically compatible with
            % the required tuner delays.
            if totalTunerDelay > totalNonTxTime
                totalNonTxTime = totalNonTxTime + (totalTunerDelay-totalNonTxTime);
                %error(['Requested dutyCycle is too high to accommodate the ' ...
                %       'required tuner delays.']);
            end
            
            % Remaining non-TX time that can be distributed as random listening.
            totalRandomListeningTime = totalNonTxTime - totalTunerDelay;
            
            % Randomly distribute the available listening time.
            % There are numBursts+1 possible listening intervals:
            % listen | B1 | listen+tuner | B2 | ... | listen+tuner | BN | listen
            randomWeights = rand(numBursts+1,1); % #numBursts+1 positive random 
            % numbers between (0;1); we take #numBursts+1 intervals since there is
            % a small RX interval before the transmission of the very first burst
            randomListeningIntervals = totalRandomListeningTime * ...
                randomWeights / sum(randomWeights); % normalization to 1 --> 
            % --> sum(randomListeningIntervals) = totalRandomListeningTime.
            % Construct the burst start times.
            startTimes = zeros(numBursts,1); % starttime TX cumulative
            % Initial random listening interval before the first burst.
            startTimes(1) = randomListeningIntervals(1);    
            for k = 2:numBursts
            
                startTimes(k) = ...
                    startTimes(k-1) ...                  % previous burst start
                    + burstTxTime ...                    % previous burst duration
                    + tunerDelayPerBurst(k) ...          % possible frequency switching
                    + randomListeningIntervals(k);       % random listening time
            
            end    
    
        end         
        
        % End time of each burst.
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
        txRRC = comm.RaisedCosineTransmitFilter( ...
            'Shape','Square root', ...
            'RolloffFactor',rolloff, ...
            'FilterSpanInSymbols',filterSpan, ...
            'OutputSamplesPerSymbol',samplesPerSymbol); % how many samples it 
        % outputs per symbol the RRC
        totalTime = max(endTimes) + randomListeningIntervals(end); % final time for transmitting the whole packet 
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
            burstSymbols = txSCFDEBurstsInterferer(:,k); % from where we sample
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