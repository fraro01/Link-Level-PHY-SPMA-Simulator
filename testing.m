clc;
clear;
close all;

%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
* 
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           TRANSMITTER                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Originate the Packet
% Set the length of the packet in bits
packetLength = 1331;
% Generate random bits (0 or 1)
packet = randi([0 1], packetLength, 1); %column vector not row for applying 
% turbo coding


%% Turbo Coding
% TODO! add trellis and numeiterations
% NOTE! study the turbo coding functionings
% Interleaver turbo
interleaverIndices = randperm(packetLength);
% Turbo encoder rate 1/3, (it is the default)
% Use as "InterleaverIndices" the one defined above
turboEnc = comm.TurboEncoder('InterleaverIndices', interleaverIndices);
% Coding
encodedPacket = turboEnc(packet);


%% Divide the Original Packet into Bursts
% ASSUMPTION! we assume the number of bits per packet fixed, hence packet
% length fixed
numBursts = 28;
% length of the turbo encoded packet
numEncodedBits = length(encodedPacket);
% "integer superior part" taken --> exceeding bits --> padding bits
bitsPerBurst = ceil(numEncodedBits / numBursts);
numPaddingBits = bitsPerBurst*numBursts - numEncodedBits;
% padding the encoded packet with #numPaddingBits zeros added at the end
% NOTE! the Rx MUST "unpad" the received packet
encodedPacketPadded = [encodedPacket; zeros(numPaddingBits,1)];
% (bitsPerBurst x numBursts) matrix keeping the order of bits
bursts = reshape(encodedPacketPadded, bitsPerBurst, numBursts);


%% 48 sync reference bits, 48 phase reference bits, 5 indexing bits
% the sync and the reference bits remain fixed for all bursts
rng(1); % fixed seed to keep fixed only the two random variables below
syncReference  = randi([0 1], 48, 1);   % 48 synchronization bits
phaseReference = randi([0 1], 48, 1);   % 48 phase reference bits
rng('shuffle'); % come back to random sequence for all rest of the script
% integer superior part of the log2 --> add minimum number of indexing bits
numIndexBits = ceil(log2(numBursts));   % 5 bits for indexing 28 bursts
% matrix of zeros
framedBursts = zeros(bitsPerBurst + 48 + 48 + numIndexBits, numBursts);
for k = 1:numBursts
    % decimal to binary, with numIndexBits as #bits of result, most
    % significant bit is put left (used in TLC), then we transpose it: .'
    burstIndexBits = de2bi(k-1, numIndexBits, 'left-msb').';
    framedBursts(:,k) = [ % take all rows of cloumn k
    % vertical concatenation of three vectors
        syncReference;        % 48 bit
        phaseReference;       % 48 bit
        burstIndexBits;       % 5 bit: burst identifier
        bursts(:,k)           % 144 bit coded, from column k
    ];
end


%% 8-PSK Modulation
M = 8;
bitsPerSymbol = log2(M); % number of bits per symbol transmitted
bitsPerFramedBurst = size(framedBursts, 1);
% padding in order to render each burst divisible by bitPerSymbol
numPadBits8PSK = mod(bitsPerSymbol - ...
    mod(bitsPerFramedBurst, bitsPerSymbol), bitsPerSymbol);
% add the above padding
framedBurstsPadded = [ ...
    framedBursts;
    zeros(numPadBits8PSK, numBursts)
];
bitsPerBurst8PSK = size(framedBurstsPadded, 1);
symbolsPerBurst = bitsPerBurst8PSK / bitsPerSymbol;
modulatedBursts = zeros(symbolsPerBurst, numBursts); % matrix
for k = 1:numBursts
    bits_k = framedBurstsPadded(:,k); % k column of the 8-psk padded burst
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


%% SC-FDE Modulation, Cyclic Prefix 
cpLength = 16; % TODO! this must be set acccording to the channel length
% adding the space for the cyclic prefix
txSCFDEBursts = zeros(symbolsPerBurst + cpLength, numBursts);
for k = 1:numBursts
    % already 8-PSK modulated bursts
    burstSymbols = modulatedBursts(:,k);
    % getting the cyclic prefix
    cyclicPrefix = burstSymbols(end-cpLength+1:end);
    % vertical concatenation of the two vectors: CP + modulatedBurst
    txSCFDEBursts(:,k) = [ ...
        cyclicPrefix;
        burstSymbols
    ];
end


%% Time-Frequency Hopping
% NOTE! we are also distribuiting the burst over a pseudo-casual time
% pattern.
% symbol rate = 13 MHz; maximum 17 MHz, use it as a hyperparameter
% parameters
numFrequencyChannels = 4;
symbolRate = 13e6; % 13MHz
symbolsPerBurstWithCP = size(txSCFDEBursts,1);
dutyCycle = 0.25; 
fMin = 1435e6;
fMax = 1518e6;
% timing
burstTxTime = symbolsPerBurstWithCP / symbolRate;     % seconds
slotDuration = burstTxTime / dutyCycle;               % seconds
totalPacketTime = numBursts * slotDuration;           % seconds
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
%  transmission of a burst to create tails that exceed their own allocated 
% time slot.
shapedBurstTime = ...
    (symbolsPerBurstWithCP + filterSpan) / symbolRate; 
maxOffset = slotDuration - shapedBurstTime;


intraSlotOffset = rand(numBursts,1) * maxOffset; % uniform distribution!
% array of values between 1 and #numFrequencyChannles of length #numBursts
frequencyChannels = randi( ...
    [1 numFrequencyChannels], ...
    numBursts, ...
    1);
carrierForBurst = carrierFrequencies(frequencyChannels).';
startTimes = (timeSlots - 1).' * slotDuration + intraSlotOffset; % arrray
endTimes = startTimes + burstTxTime; % array

%{
% display hopping table
hoppingTable = table( ...
    (1:numBursts).', ...
    timeSlots.', ...
    startTimes*1e6, ...
    endTimes*1e6, ...
    frequencyChannels, ...
    carrierForBurst/1e6, ...
    'VariableNames', ...
    {'BurstID','TimeSlot','Start_us','End_us', ...
     'FrequencyChannel','Carrier_MHz'});
disp(hoppingTable);
%}

% TODO! do we really need the following structure?
%{
% Build transmission structure
for k = 1:numBursts
    transmission(k).burstID = k; %#ok<SAGROW>
    transmission(k).timeSlot = timeSlots(k); %#ok<SAGROW>
    transmission(k).startTime_s = startTimes(k); %#ok<SAGROW>
    transmission(k).endTime_s = endTimes(k); %#ok<SAGROW>
    transmission(k).frequencyChannel = frequencyChannels(k); %#ok<SAGROW>
    transmission(k).carrierFrequency_Hz = carrierForBurst(k); %#ok<SAGROW>
    transmission(k).signal = txSCFDEBursts(:,k); %#ok<SAGROW>
end
%}

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
% Configuration of Root Raised Cosine TX Filter
rolloff = 0.35; % how much extra band we use
% NOTICE! each RRC burst has null-to-null bandwidth: symbolRate*(1+rolloff)
% IF we transmit simultaneously multiple bursts per time slot we could
% experience ACI, since our channel spacing is roughly 13.83 MHz --> set
% rolloff accordingly
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
totalFilterDelaySymbols = filterSpan;
totalTime = max(endTimes); % final time for transmitting the whole packet 
% (last burst)
symbolsPerBurstWithCP = size(txSCFDEBursts,1);
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
    burstSymbols = txSCFDEBursts(:,k); % from where we sample
    % Add zeros to flush the TX/RX filter memory
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          END TRANSMITTER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           SIMPLE CHANNEL                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% example of simple multipath channel, column vector, three taps,
% IT WORKS ON THE SAMPLES, NOT ON THE SYMBOLS, the channel sees samples
h = [0.9;       % LOS
     0.4+0.2j;  % echo 1
     0.2-0.1j]; % echo 2

% Apply channel response as a convolution without retroaction as a FIR 
% between the Tx signal and the channel.
rxSignal = filter(h, 1, txSignal);

% AWGN
SNRdB = 100; % signal to noise ratio in dB (P_signal/P_noise), but power 
% calculus are usually performed linearly...
SNRlin = 10^(SNRdB/10); % converts the above SNR from dB to linear value 
% "the signal is SNRlin times greater than the noise."

% Measure actual signal power only after channel effects
% Noise variance necessary to obtain the desired SNR set above
active = abs(rxSignal) > 1e-12; % it stores 1 and 0 (True/False)
signalPower = mean(abs(rxSignal(active)).^2); % calculates signal power
% only for non-zeros symbols, thanks to threshopld set above (1e-12)
noiseVar = signalPower / SNRlin; % NOTE! the noiseVar, based on the above 
% formulas is the P_noise, For zero-mean noise, average power equals 
% variance.
noise = sqrt(noiseVar/2) * ... % complex noise, divided by 2, because the variance gets divided between Re and Im
    (randn(size(rxSignal)) + 1j*randn(size(rxSignal))); % complex noise
% randn generates gaussian noise with mean 0 and variance 1

rxSignal = rxSignal + noise; % add noise to the attenuated received signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        END SIMPLE CHANNEL                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             RECEIVER                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Conversion in reception of EQUIVALENT-BASEBAND signal with RRC filter 
%% Realistic Burst Synchronization Using Preamble Correlation
%% Coarse AGC --> ADC --> Fine AGC
% The Rx knows:
%   - "Slot" duration
%   - Frequency hopping sequence
%   - Sync symbols
% The Rx DOES NOT know:
%   - startTimes(k) hence the start time of transmission of each burst
%   - intraSlotOffset(k) the offset of the bursts inside their own slot
%   - the exact sample of burst start
% numbers of complete symbols inside the syncReference
numSyncSymbols = floor(length(syncReference) / bitsPerSymbol);
% The first #numSyncSymbols of each burst are symbols for synchronization, 
% and are equal for all bursts
knownSyncSymbols = modulatedBursts(1:numSyncSymbols, 1);
% Normalization of the preamble, so to render the metric less dependent on
% the power received
knownSyncSymbolsNormalized = knownSyncSymbols / ...
    sqrt(sum(abs(knownSyncSymbols).^2));

% RRC filter RX
% NOTE!:
% We do not decimate so far, because we still don't know the correct
% temporal phase. Hence to output of this filter mantains #samplesPerSymbol
% samples per symbol
rxRRCsync = comm.RaisedCosineReceiveFilter( ... % pulse shaping filter RX
    'Shape', 'Square root', ...
    'RolloffFactor', rolloff, ...
    'FilterSpanInSymbols', filterSpan, ...
    'InputSamplesPerSymbol', samplesPerSymbol, ...
    'DecimationFactor', 1);
% where to store the syncohronized bursts
rxSCFDEBursts = zeros(symbolsPerBurstWithCP, numBursts);
% Set the threshold in between 0 e 1. The higher the amount of sync symbols 
% the higher could be set, and viceversa. e.g. with only 5 symbols --> 0.55
syncDetectionThreshold = 0.55;

% COARSE and FINE AGC configuration parameters (2.1. and 8.1.)
% Approximate desired power at the output of the Coarse AGC.
% A power equal to 0.25 correspondes to an RMS amplitude of sqrt(0.25)=0.5
% We do not use power=1, so to leave a certain margin to an eventual
% saturation of RX.
coarseTargetPower = 0.25;
% Coarse gain limits.
% the Coarse AGC can attenuate up to a factor of: 0.05
% or amplify up to a factor of: 30.
coarseMinGain = 0.05;
coarseMaxGain = 30;
% Fine gain limits.
% The Fine AGC has to perform only a small final correction, because the
% main part has been already corrected by the Coarse AGC.
fineMinGain = 0.5;
fineMaxGain = 2;

% NOTE! this is the main external cycle that iterates over all the bursts
for k = 1:numBursts
    reset(rxRRCsync); % clean the filter

    % 1. EXTRACTION OF THE ENTIRE TIME SLOT
    % We do not use startTimes(k).
    % the RX looks at the whole slot k.
    slotStartTime = (k-1) * slotDuration;
    slotEndTime   = k * slotDuration;
    slotStartSample = floor(slotStartTime * Fs) + 1;
    slotEndSample   = min(ceil(slotEndTime * Fs), length(rxSignal));
    rxSlotUpconverted = rxSignal(slotStartSample : slotEndSample);

    % 2. DOWNCONVERSION
    % RX knows the FH sequence. But doesn't know the exact time of start 
    % burst
    globalSampleIndices = (slotStartSample-1 : slotEndSample-1).';
    tSlot = globalSampleIndices / Fs; % Fs = [samples/second], do the math
    fRel = relativeCarrierFrequencies(frequencyChannels(k)); % FH scheme
    % Downconversion bring the entire time slot back to baseband
    % raw because we still have to appli Coarse AGC
    rxSlotBasebandRaw = rxSlotUpconverted .* exp(-1j * 2*pi * fRel *tSlot);

    % 2.1. SIMPLE COARSE AGC
    % NOTE! it is not causal since the code reads the whole slot, hence it
    % takes advantage of samples that will actually arrive in the future,
    % but a real recevier cannot do this in reality in real time!
    % instead it would update the gain sample by sample.
    % AGC Coarse is required in order to bring the received signal in an
    % amplitude interval compatible with the front-end of the RX, before
    % synchronization. It must guarantee that the analog signal does not
    % overcome the full scale of the ADC.
    % Put the coarse AGC after downconversion since it's easier in the
    % equivalent-baseband, because AGC depends on the modulus of the
    % signal, the frequency traslation does not modify amplitude.
    % Calculate the istantaneous power of all the samples of the slot.
    slotSamplePower = abs(rxSlotBasebandRaw).^2;
    % The slot contains:
    % - a part with the burst;
    % - a part mainly with noise.
    % In order to avoid that the long noise part lowers too much the 
    % estimate, we only consider samples with relatively high power.
    % the 75th percentile divides the samples as follow:
    % - the 75% has power lower than the threshold;
    % - the 25% has power higher than the threshold.
    % NOTE! adjust the percentile based on the Duty Cycle set.
    powerThreshold = prctile(slotSamplePower,75);
    % Select the samples probaly belonging to the burst.
    activeSamples = slotSamplePower >= powerThreshold; % array of: 1/0
    % Safety check.
    % Normally activeSamples contains at least 25% of samples.
    % If for some unforseen reasons no sample is selected, we use the whole 
    % slot to calculate the average.
    if any(activeSamples) % if at least there is a: 1 in the array --> True
        measuredCoarsePower = mean( ...
            slotSamplePower(activeSamples));
    else
        measuredCoarsePower = mean(slotSamplePower);
    end
    % Calculate the required gain to obtain the objective power.
    % Remember: power is the square of the amplitude --> sqrt()
    % If:
    % Pout = gain^2 * Pin
    % Then:
    % gain = sqrt(Pout / Pin)
    coarseGain = sqrt( ...
        coarseTargetPower / ...
        max(measuredCoarsePower,1e-12));
    % Limit the gain to the allowed interval.
    coarseGain = min( ...
        max(coarseGain,coarseMinGain), ...
        coarseMaxGain);
    % Apply the same gain to the whole slot.
    rxSlotBaseband = coarseGain * rxSlotBasebandRaw;

    % 2.1.1. Simple ADC Model
    adcNumBits = 10; % maximum number of available bits to represent values
    adcFullScale = 1; % full scale, hence range: [-1;+1], for clipping
    % ADC maximum positive integer code.
    adcMaxCode = 2^(adcNumBits-1) - 1; % maximum representable number of 
    % ADC, in order to be symmetric to 0 --> -511,..,0,...,+511 (total 1023
    % not 1024).
    % Quantization step. distance between two consecutive quantized levels
    adcStep = adcFullScale / adcMaxCode;    
    % Separate I and Q components.
    rxADCReal = real(rxSlotBaseband);
    rxADCImag = imag(rxSlotBaseband);    
    % Input clipping. [-1;+1], very clever with only one row of code
    rxADCReal = min(max(rxADCReal,-adcFullScale),adcFullScale);
    rxADCImag = min(max(rxADCImag,-adcFullScale),adcFullScale);    
    % Uniform quantization. Uniform since the quantization step is equal.
    % how many quantization steps corresponds the input value, round them
    % to the nearest integer, then go back to the normalized value
    rxADCReal = adcStep * round(rxADCReal/adcStep);
    rxADCImag = adcStep * round(rxADCImag/adcStep);    
    % Complex ADC output. after clipping and quantization
    rxSlotBasebandADC = rxADCReal + 1j*rxADCImag;    

    % 3. RCC FILTER WITHOUT DECIMATION
    % The exit still keeps the #samplesPerSymbol samples per symbol
    rxSlotMatched = rxRRCsync(rxSlotBasebandADC);

    % 4. CORRELATION METRIC CALCULUS
    % For each possible sample q, we take:
    % q,
    % q + samplesPerSymbol,
    % q + 2*samplesPerSymbol,
    % ...
    % So to obtain a possible sequence of received symbols
    % We confront it with: knownSyncSymbols.
    numMatchedSamples = length(rxSlotMatched);
    % the first sync symbol must have as header before itslef the whole CP,
    % hence we go after the CP for position correlation.
    firstCandidate = cpLength * samplesPerSymbol + 1;
    % After the first sync symbol there must be available:
    %   - All sync symbols;
    %   - The whole rest part of the burst.
    % In other words we are extracting the entire burst from its "slot".
    lastCandidate = numMatchedSamples ...
        - (symbolsPerBurstWithCP - cpLength - 1) * samplesPerSymbol; % here 
    % we are taking the POSSIBLE LAST POSITION of the burst inside its slot 
    % in terms of samples, starting directly AFTER the CP.
    numCandidates = lastCandidate - firstCandidate + 1; % number of 
    % candidated POSITIONS for starting the preamble (in samples)
    syncMetric = zeros(numCandidates, 1);
    for candidateIdx = 1:numCandidates
        q = firstCandidate + candidateIdx - 1;
        % possible samples of sync symbols 
        syncSampleIndices = q + ...
            (0:numSyncSymbols-1).' * samplesPerSymbol;
        % Possible samples candidated
        receivedSyncCandidate = ...
            rxSlotMatched(syncSampleIndices);
        % Local energy of the received sequence
        receivedEnergy = sum( ...
            abs(receivedSyncCandidate).^2);
        % if the energy is above this threshold we proceed, otherwise no
        if receivedEnergy > 1e-12
            % normalize it, to compare with the normalized known sequence
            receivedSyncNormalized = ...
                receivedSyncCandidate / sqrt(receivedEnergy);
            % Complex correlation normalized.
            % the absolute value, renders the metric insensible to a phase
            % rotation in common.
            correlationValue = ...
                knownSyncSymbolsNormalized' * receivedSyncNormalized; % the 
            % knownSyncSymbolsNormalized' makes the conjugated transpose
            syncMetric(candidateIdx) = abs(correlationValue)^2; % in order 
            % to cope with phase rotation of the channel, the modulus
            % reimans unchanged while the angle could change
        end
    end

    % 5. LOOK FOR THE MAXIMUM CORRELATION
    % take the maximum and its position
    [bestMetric, bestCandidateIdx] = max(syncMetric);
    % position of the first sync symbol in the output of the adapted filter
    firstSyncSample = firstCandidate + bestCandidateIdx - 1;

    % 6. DETECTION DECISION
    if bestMetric < syncDetectionThreshold
        warning( ...
            ['Burst %d not detected: synchronization metric ' ...
             '= %.3f.'], ...
            k, ...
            bestMetric);
        continue;
    end

    % 7. FIND WHERE THE CP START
    % The preamble is the first element after the CP
    firstCPSample = firstSyncSample - cpLength * samplesPerSymbol;

    % 8. EXTRACTION OF THE SYMBOLS OF THE BURST
    % now we know both:
    %   - tepmoral position;
    %   - correct sampling phase.
    % We sample a value each #samplePerSymbols samples
    burstSampleIndices = firstCPSample + ...
        (0:symbolsPerBurstWithCP-1).' * samplesPerSymbol;

    if burstSampleIndices(1) < 1 || ... % || (or)
        burstSampleIndices(end) > numMatchedSamples
        warning( ...
            'Burst %d detected, but not entirely contained in the slot.', ...
            k);
        continue;
    end

    rxBurstSymbols = ...
    rxSlotMatched(burstSampleIndices);

    % 8.1. SIMPLE FINE AGC 
    % AGC Fine is required in order to finalize the coarse AGC, over the
    % precise found synchronization symbols.
    % Fine AGC must be implemented once we have understood which samples
    % belong to the burst, the timing must be known.
    % The symbols known by the receiver are the ones generated by:
    % - 48 synchronization reference bits;
    % - 48 phase reference bits.
    % Hence in total we have:
    % 48 + 48 = 96 knwon bit.
    % With 8-PSK each symbol contains 3 bits:
    % 96 / 3 = 32 known symbols.
    numFinePilotSymbols = floor( ...
        (length(syncReference) + length(phaseReference)) ...
        / bitsPerSymbol);    
    % The transmitted pilots are the first symbols of the burst after the 
    % CP. Since Sync and Phase reference are equal fot all bursts,
    % We can extract them from the first column of: modulatedBursts.
    knownFinePilotSymbols = ...
        modulatedBursts(1:numFinePilotSymbols,1);    
    % rxBurstSymbols contains at the beginning the cyclic prefix.
    % So the first pilot symbol is in position:
    % cpLength + 1
    finePilotStart = cpLength + 1;    
    % And the last pilot symbol known.
    finePilotEnd = ...
        finePilotStart + numFinePilotSymbols - 1;    
    % Extract the pilots received.
    receivedFinePilotSymbols = ...
        rxBurstSymbols(finePilotStart:finePilotEnd);    
    % Average power of the received pilots.
    measuredFinePilotPower = mean( ...
        abs(receivedFinePilotSymbols).^2);  
    % Average power of the transmitted pilots.
    % for a normal 8-PSK modulation, should be around 1.
    % In any case it is better to directly calculate it.
    targetFinePilotPower = mean( ...
        abs(knownFinePilotSymbols).^2);
    % Calculate the fine gain.
    fineGain = sqrt( ...
        targetFinePilotPower / ...
        max(measuredFinePilotPower,1e-12));
    % Limit the correction
    fineGain = min( ...
        max(fineGain,fineMinGain), ...
        fineMaxGain);
    % apply the same gain to the whole burst.
    rxBurstSymbols = fineGain * rxBurstSymbols;
    % Save the burst after Coarse and Fine AGC.
    rxSCFDEBursts(:,k) = rxBurstSymbols;
end

%testing
%{
errMax = max(abs(rxSCFDEBursts(:) - txSCFDEBursts(:)));
errRMS = sqrt(mean(abs(rxSCFDEBursts(:) - txSCFDEBursts(:)).^2));
fprintf('Max error: %.3e\n', errMax);
fprintf('RMS error: %.3e\n', errRMS);
%}


%% PILOT-BASED Channel Estimation using Sync + Phase reference bits and 
%% SC-FDE Demodulation with ZF Equalization
N = symbolsPerBurst;
Ncp = cpLength;
% Known pilot bits available in each burst
knownPilotBits = [syncReference; phaseReference];
numPilotSymbols = floor(length(knownPilotBits)/bitsPerSymbol);
knownPilotSymbols = modulatedBursts(1:numPilotSymbols, 1);
numPilotSymbols = length(knownPilotSymbols);

% NOTE! the pilot tones (symbols) must be >> number of taps (Lest). 
% Currenlty the way it is does not really satisfy this. there are two 
% possible ways to operate:
% 1. Use dedicated Pilot Bursts as used in the current LinkLevelPHY.m file
% 2. Augment the length of the sync and reference bits, at least covering a
% total of 16-32 symbols.
% by doing so you will notice a strong decrease in numErrors8PSK even with
% higher AWGN.
% Length ESTimate
Lest = 3; % length of the discrete impulsive response we want to estimate. 
% it represents #Lest taps at symbols level, hence we estimate the channel 
% memory length in symbols, hence it is totally different from what we did 
% in the Channel implementation before the RRC RX. Indeed we are doing this
% after the RRC RX

rxEqualizedBursts = zeros(symbolsPerBurst, numBursts); % equalized bursts
H_est_bursts = zeros(N, numBursts); % channel estimation in frequency
h_est_bursts = zeros(Lest, numBursts); % channel estimation in time
for k = 1:numBursts % cycle over all bursts
    rxBurstCP = rxSCFDEBursts(:,k);
    % Remove CP
    rxBurst = rxBurstCP(cpLength+1 : cpLength+symbolsPerBurst);
    % Received known pilot portion
    rxPilot = rxBurst(1:numPilotSymbols);
    % Build Least Squares matrix:
    % rxPilot(n) ≈ h(1)*x(n) + h(2)*x(n-1) + ... means that we are
    % modelling the received symbols as convoluted with the channel
    numRows = numPilotSymbols - Lest + 1; % how many equations we construct 
    % using the pilots in order to estimate the taps (Lest)
    A = zeros(numRows, Lest); % "A" matrix containing the known pilots for 
    % channel estimate
    y = zeros(numRows, 1); % vector of the corresponding received samples

    row = 1; % initialize the index of the current row of matrix A
    for n = Lest:numPilotSymbols % to build a row we need x(n), x(n-1), 
        % x(n-2). this is the reason why we start from Lest
        for ell = 0:Lest-1 % cycle over the delays of the channel
            A(row, ell+1) = knownPilotSymbols(n-ell); % "A" contains the 
            % known transmitted pilots put in a convoluted form
        end
        y(row) = rxPilot(n); % pilot symbol received at time n
        row = row + 1;
    end

    % Least Squares channel estimate, since we have more equations than
    % unknown variables.
    h_hat = A \ y; % solves the system: A * h_hat ≈ y, hence we get the 
    % estimates of the channel from our knonw pilots
    h_est_bursts(:,k) = h_hat; % saves the times response estimates by 
    % burst k
    % Convert time-domain channel estimate to frequency domain
    h_hat_padded = [h_hat; zeros(N-Lest,1)];
    H_hat = fft(h_hat_padded, N); % channel frequency response
    H_est_bursts(:,k) = H_hat; % saves the frequency response estimates by 
    % burst k
    % SC-FDE Equalization
    Y = fft(rxBurst, N);
    % ZF with protection against very small H
    epsFloor = 1e-6;
    W = conj(H_hat) ./ (abs(H_hat).^2 + epsFloor);
    Xhat = W .* Y; % equalization
    rxEqualizedBursts(:,k) = ifft(Xhat, N);
end


%% 8-PSK Demodulation
rxModulatedBursts = rxEqualizedBursts; 
% due to channel effects
rxFramedBurstsPadded = zeros(symbolsPerBurst*bitsPerSymbol, numBursts);

for k = 1:numBursts
    rxSymbols_k = rxModulatedBursts(:,k);
    rxSymbolIndices = pskdemod( ...
        rxSymbols_k, ...
        M, ...
        pi/M, ...
        'gray');
    rxBitGroups = de2bi(rxSymbolIndices, bitsPerSymbol, 'left-msb');
    rxFramedBurstsPadded(:,k) = reshape(rxBitGroups.', [], 1);
end

% Remove 8-PSK padding
rxFramedBursts = rxFramedBurstsPadded( ...
    1:end-numPadBits8PSK, ...
    :);
% check for the errors Tx - Rx
%{%}
numErrors8PSK = sum(rxFramedBursts(:) ~= framedBursts(:))


%% Turbo decoder
% Random burst loss simulation
numLostBursts = 0;                     % choose manually
lostBurstIdx = randperm(numBursts, numLostBursts);
% Recover encoded packet from bursts
headerLength = length(syncReference) + ...
               length(phaseReference) + ...
               numIndexBits;
rxBurstsPayload = rxFramedBursts(headerLength+1:end,:); % remove headers
rxEncodedPacketPadded = reshape(rxBurstsPayload,[],1); % stack into column
rxEncodedPacket = rxEncodedPacketPadded(1:numEncodedBits); % leave padding

% Ideal Channel: No noise
% The decoder does not want bit 0/1, instead soft values.
% Typical convention:
% bit 0 -> -1
% bit 1 -> +1
rxSoft = 2*double(rxEncodedPacket) - 1;

payloadBitsPerBurst = size(rxBurstsPayload,1);
for b = lostBurstIdx % iterate over erased bursts
    firstBit = (b-1)*payloadBitsPerBurst + 1; % take also the first burst, 
    % hence 0, furthermore MATLAB starts from 1 not 0 
    lastBit  = min(b*payloadBitsPerBurst, numEncodedBits); % handle also 
    % the case of the last burst in which there was padding (minimum len)
    rxSoft(firstBit:lastBit) = 0; % set these bits to zero
end

% numIterations = 6; the default is 6
turboDec = comm.TurboDecoder( ...
    'InterleaverIndices', interleaverIndices);
    %'NumIterations', numIterations);
decodedPacket = turboDec(rxSoft);

% Just to compare
%{%}
numFinalErrors = sum(packet ~= decodedPacket)