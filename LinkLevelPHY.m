clc;
clear;
close all;

%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
* 
%}

%% Originate the Packet
% Set the length of the packet in bits
packetLength = 1331;
% Generate random bits (0 or 1)
packet = randi([0 1], packetLength, 1); %column vector not row for applying 
% turbo coding

%% Turbo Coding
% Interleaver turbo
interleaverIndices = randperm(packetLength);
% Turbo encoder rate 1/3, (it is the standard)
% Use as "InterleaverIndices" the one defined above
turboEnc = comm.TurboEncoder('InterleaverIndices', interleaverIndices);
% Coding
encodedPacket = turboEnc(packet);

% TODO! add trellis and numeiterations
% NOTE! study the turbo coding functionings
%{
%% Ideal Channel: No noise
%The decoder does not want bit 0/1, instead soft values.
% Typical convention:
% bit 0 -> -1
% bit 1 -> +1
rxSoft = 2*encodedPacket - 1;

%% Turbo decoder
%numIterations = 6;
turboDec = comm.TurboDecoder( ...
    'InterleaverIndices', interleaverIndices);
    %'NumIterations', numIterations);

decodedPacket = turboDec(rxSoft);

% Just to compare
numErrors = sum(packet ~= decodedPacket);
%}

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
% (144 x numBursts) matrix keeping the order of bits
bursts = reshape(encodedPacketPadded, bitsPerBurst, numBursts);

%% 16 sync reference bits, 16 phase reference bits, 5 indexing bits 
% the sync and the reference bits remain fixed for all bursts
syncReference  = randi([0 1], 16, 1);   % 16 synchronization bits
phaseReference = randi([0 1], 16, 1);   % 16 phase reference bits
% integer superior part of the log2 --> add minimum number of indexing bits
numIndexBits = ceil(log2(numBursts));   % 5 bits for indexing bursts
% matrix of zeros
framedBursts = zeros(bitsPerBurst + 16 + 16 + numIndexBits, numBursts);
for k = 1:numBursts
    % decimal to binary, with numIndexBits as #bits of result, most
    % significant bit is put left (used in TLC), then we transpose it: .'
    burstIndexBits = de2bi(k-1, numIndexBits, 'left-msb').';
    framedBursts(:,k) = [ % take all rows of cloumn k
    % vertical concatenation of three vectors
        syncReference;        % 16 bit
        phaseReference;       % 16 bit
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

%{
%% 8-PSK Demodulation
rxModulatedBursts = modulatedBursts; % ...ideally they should be different 
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
numErrors8PSK = sum(rxFramedBursts(:) ~= framedBursts(:));
%}

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

%{
%% SC-FDE Demodulation, Basic Channel + MMSE
% example of simple multipath channel, column vector, three taps
h = [0.9;      % LOS
    0.4+0.2j;  % echo 1
    0.2-0.1j]; % echo 2

% noise
SNRdB = 20; % signal to noise ratio in dB
SNRlin = 10^(SNRdB/10); % converts the above SNR from dB to linear value
noiseVar = 1/SNRlin; % noise variance

% Channel + Noise
rxSCFDEBursts = zeros(size(txSCFDEBursts));
for k = 1:numBursts
    txBurst = txSCFDEBursts(:,k);
    rxBurst = filter(h, 1, txBurst); % channel response as a convolution, 
    % without retroaction as a FIR between the Tx signal and the channel.
    noise = sqrt(noiseVar/2) * ... % complex noise, divided by 2, because the variance gets divided between Re and Im
        (randn(size(rxBurst)) + 1j*randn(size(rxBurst))); % complex noise
    rxSCFDEBursts(:,k) = rxBurst + noise;
end

% frequency response of the channel, padded so to get the same number of
% bins in output desired for the element wise operations per symbolPerBurst
H = fft(h, symbolsPerBurst);

% MMSE equalizer
% NOTE!  in this first implementation we do not estimate the channel!
% ./ element-wise division, abs(H).^2 squared modulus of H --> channel
% power for each frequency
Wmmse = conj(H) ./ (abs(H).^2 + noiseVar);

% Remove CP + equalization MMSE
rxEqualizedBursts = zeros(symbolsPerBurst, numBursts);
for k = 1:numBursts
    rxBurstCP = rxSCFDEBursts(:,k);
    % Remove cyclic prefix
    rxBurst = rxBurstCP(cpLength+1 : cpLength+symbolsPerBurst);
    % FFT
    Y = fft(rxBurst, symbolsPerBurst);
    % Frequency-domain MMSE equalization, .* element-wise multiplication
    Xhat = Wmmse .* Y; % burt estimation in the frequency domain
    % IFFT: go back to the symbol 8-PSK
    rxEqualizedBursts(:,k) = ifft(Xhat, symbolsPerBurst);
end
%}

%% Time-Frequency Hopping
% NOTE! we are also distribuiting the burst over a pseudo-casual time
% pattern.
% symbol rate = 13 MHz; maximum 17 MHz, use it as a hyperparameter
% parameters
numFrequencyChannels = 6;
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
% Random offset inside each slot
maxOffset = slotDuration - burstTxTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
This block is required in order to better cope with the tails introduced by
the Pulse Shaping Filters. But since with a simple simulation I found out
that actually the results do not really change, I have decided to comment
it and keep the simple older implementation, in which a burst could be
transmitted even in the final edge of its timeSlot.

shapedBurstTime = ...
    (77 + 12) / 13e6;
maxOffset = slotDuration - shapedBurstTime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

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

%% Final transmitted EQUIVALENT-BASEBAND signal with Pulse Shape Filter RRC
% NOTE! Due to Nyquist and thanks to the fact that we are simulating in
% baseband, we must use: samplesPerSymbol > 2*(fmax_occupied)/symbolRate
% in which fmax_occupied = [max(relativeCarrierFrequencies) +
% (relativeCarrierFrequencies(6) - relativeCarrierFrequencies(5))]
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
Fs = symbolRate * samplesPerSymbol; % sampling frequency
% Configuration of Root Raised Cosine TX Filter
rolloff = 0.35; % how much extra band we use
% NOTICE! each RRC burst has null-to-null bandwidth: symbolRate*(1+rolloff)
% IF we transmit simultaneously multiple bursts per time slot we could
% experience ACI, since our channel spacing is roughly 13.83 MHz --> set
% rolloff accordingly
filterSpan = 12; % filter duration in symbols, NOTE! the longer it lasts
% over time, the sharper it is over the frequencies.
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

%{%}
%% Conversion in reception of EQUIVALENT-BASEBAND signal with RRC filter
% pulse shaping filter Rx
rxRRC = comm.RaisedCosineReceiveFilter( ...
    'Shape','Square root', ...
    'RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',filterSpan, ...
    'InputSamplesPerSymbol',samplesPerSymbol, ...
    'DecimationFactor',samplesPerSymbol); % after filtering keep one sample 
    % every #samlpesPerSymbol

rxSignal = txSignal; % ASSUMPTION! ideal channel
rxSCFDEBursts = zeros(symbolsPerBurstWithCP, numBursts);
for k = 1:numBursts
    reset(rxRRC); % clean the filter
    % Same time window used at transmitter
    startSample = round(startTimes(k)*Fs) + 1;
    % Notice that this is different from Tx
    samplesPerShapedBurst = ...
    (symbolsPerBurstWithCP + extraSymbols) * samplesPerSymbol;
    endSample = startSample + samplesPerShapedBurst - 1;
    % Extract received burst samples
    rxBurstUpconverted = rxSignal(startSample:endSample);
    % Absolute temporal axis
    n = (startSample-1:endSample-1).';
    t = n/Fs;
    % Relative carrier frequency used for this burst
    fRel = relativeCarrierFrequencies(frequencyChannels(k));
    % Downconversion: bring burst back to baseband
    rxBurstBaseband = rxBurstUpconverted .* exp(-1j*2*pi*fRel*t);
    % Matched filter + downsampling (from samples we go back to symbols)
    % The matched filtering is declined to the waveform not to preamble!
    rxBurstSymbolsPadded = rxRRC(rxBurstBaseband);
    % Remove total TX+RX filter delay
    rxBurstSymbols = rxBurstSymbolsPadded( ...
        totalFilterDelaySymbols + 1 : ...
        totalFilterDelaySymbols + symbolsPerBurstWithCP);
    % Store recovered SC-FDE burst
    rxSCFDEBursts(:,k) = rxBurstSymbols;
end

% testing
errMax = max(abs(rxSCFDEBursts(:) - txSCFDEBursts(:)));
errRMS = sqrt(mean(abs(rxSCFDEBursts(:) - txSCFDEBursts(:)).^2));
fprintf('Max error: %.3e\n', errMax);
fprintf('RMS error: %.3e\n', errRMS);


