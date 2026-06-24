clc;
clear;
close all;

%{
* ctrl + f --> TODO! for coding 
* ctrl + f --> NOTE! for something over which we must pay attention
* ctrl + f --> ASSUMPTION! for modelling assumptions
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
% the sync and the reference bits remain fixed for each burst
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
        burstIndexBits;       % 5 bit: identificativo burst
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
symbolsPerBurstWithCP = 77;
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
rng(1234);   % results reproducibility
timeSlots = 1:numBursts; % vector column
% Random offset inside each slot
maxOffset = slotDuration - burstTxTime;
intraSlotOffset = rand(numBursts,1) * maxOffset;
% array of values between 1 and #numFrequencyChannles of length #numBursts
frequencyChannels = randi( ...
    [1 numFrequencyChannels], ...
    numBursts, ...
    1);
carrierForBurst = carrierFrequencies(frequencyChannels).';
startTimes = (timeSlots - 1).' * slotDuration + intraSlotOffset; % arrray
endTimes = startTimes + burstTxTime; % array

%{%}
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

% TODO! GENERA IL SEGNALE VERO E PROPRIO!
