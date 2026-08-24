clc;
clear;
close all;

%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
* 
%}

% NOTE! Make sure to have exactly the same name of the folders and MATLAB
% scripts as the ones I have described

% Add the two directories paths
addpath('auxiliary external functions') % for auxiliary functions
addpath('utils') % for testing functions
% transmitter, channel and receiver folders have already the '+' before the
% name, acting as a package/namespace

%% TRANSMITTER PARAMATERS
% Set the length of the packet in bits
packetLength = 1331; % bits
dutyCycle = 0.25; % percentage TX/RX
numBursts = 28;
% the sync and the reference bits remain fixed for all bursts
rng(1); % fixed seed to keep fixed only the two random variables below
syncReference  = randi([0 1], 48, 1);   % 48 synchronization bits
phaseReference = randi([0 1], 48, 1);   % 48 phase reference bits
rng('shuffle'); % come back to random sequence for all rest of the script
% Frequencies & RRC Parameters
fMin = 1435e6;
fMax = 1518e6;
numFrequencyChannels = 4;
rolloff = 0.35; % how much extra band we use for the RRC filter
filterSpan = 12; % RRC filter duration in symbols NOTE! the longer it lasts
% over time, the sharper it is over the frequencies.
% TLC formula for Symbol Rate: symbolRate = channelBandwidth/(1+rolloff)
symbolRate = ((fMax - fMin)/numFrequencyChannels)/(1 + rolloff);
samplesPerSymbol = 8; % number of samples per symbol
Fs = symbolRate * samplesPerSymbol; % sampling frequency, since we are 
% constructing the signal from the symbol to the the continuous time domain 
% we need more samples to represent one symbol.
cpLength = 10; % TODO! Must exceed the maximum channel delay in symbols
% adding the space for the cyclic prefix

%% CHANNEL PARAMETERS
% INTERFERERS
numInterferers = 0; % set the number of interferer nodes, NOTE! it 
% resembles an ALOHA network!
% BLANKING COMPENSATION
% use percentage with respect to the original value to tweak the number of 
% transmitted bursts and the respective occupied time
fractionTxBursts = round(numBursts * 1.0); % number of TX bursts while 
% Receiving, the number of TX bursts will be in any case smaller than 
% 'numBursts',
ex_fractionBurstOccupiedTime = 1.0; % time occupied by a TX burst, will be 
% multiplied by burstTxTime

%% SIMULATION PARAMETERS
% amount of packet transmissions we want to simulate
packets = 30;
% for storing values
finalEncodedPacketErroredBits = 0;
finalOriginalPacketErroredBits = 0;
lostBurst = 0;
errRMS_runs = zeros(1, packets);
errMax_runs = zeros(1, packets);
packetCompletelyReconstructed = false(1, packets);

%% SIMULATION LOOP
for i = 1:packets

    %% TRANSMITTER
    [packet, ... % original packet with random bits
     encodedPacket, ... % encoded packet above
     numEncodedBits,  ... % number of encoded bits due to Turbo Coding
     interleaverIndices, ... % for turbo decoder
     M, ... % modulation type
     bitsPerSymbol, ... 
     numPadBits8PSK, ... % padding bits for 8PSK
     symbolsPerBurst, ...
     symbolsPerBurstWithCP, ...
     modulatedBursts, ... % after 8PSK
     txSCFDEBursts, ... % after SCFDE
     frequencyChannels,...
     fCenter, ...
     relativeCarrierFrequencies, ... % equivalent baseband
     carrierForBurst, ... 
     burstTxTime, ...
     slotDuration, ...
     startTimes, ... % per burst
     txSignal ... % FINAL TRANSMITTED SIGNAL
    ] = transmitter.transmitter(...
     packetLength, ...
     dutyCycle, ... % percentage: TX/RX
     numBursts, ... % how many bursts a packet is divided into
     syncReference, ... % fixed preamble above
     phaseReference, ... % fixed preamble above
     fMin, ... % lower boundary 
     fMax, ... % upper boundary 
     numFrequencyChannels, ... 
     rolloff, ... % RRC configuration
     filterSpan, ... % RRC configuration
     symbolRate, ...
     samplesPerSymbol, ...
     Fs, ... % sampling frequency
     cpLength); % Cyclic Prefix for SCFDE
      
    %% CHANNEL
    % since we get 'burstTxTime' just from above
    fractionBurstOccupiedTime = ex_fractionBurstOccupiedTime * burstTxTime;
    [blankingMask, ... % blanking mask TxRxN
     rxSignal, ... % FINAL SIGNAL IN RX
    ] = channel.channel(...
     txSignal, ... % signal in input
     numBursts, ... 
     fMax, ...
     fCenter, ...
     carrierForBurst, ... % carrier frequencies per burst
     symbolRate, ...
     samplesPerSymbol, ...
     Fs, ...
     slotDuration, ...
     startTimes, ...
     symbolsPerBurstWithCP, ...
     txSCFDEBursts,  ...
     filterSpan, ...
     numInterferers, ...
     fractionTxBursts, ...
     fractionBurstOccupiedTime);
    
    %% RECEIVER
    [errMax, ...
     errRMS, ...
     totalErroredBitsFromEncodedPacket, ... 
     totalErroredBitsFromOriginalPacket, ...
     burstDetected, ...
    ] = receiver.receiver( ...
        rxSignal, ... % input signal
        blankingMask, ...
        numBursts, ...
        syncReference, ...
        phaseReference, ...
        cpLength, ...
        packet, ...
        encodedPacket, ...
        numEncodedBits, ...
        interleaverIndices, ...
        M, ...
        bitsPerSymbol, ...
        numPadBits8PSK, ...
        symbolsPerBurst, ...
        symbolsPerBurstWithCP, ...
        modulatedBursts, ...
        txSCFDEBursts, ...
        frequencyChannels, ...
        relativeCarrierFrequencies, ...    
        rolloff, ...
        filterSpan, ...
        symbolRate, ...
        samplesPerSymbol, ...
        Fs, ...    
        slotDuration);

    %% ACCUMULATE METRICS
    if totalErroredBitsFromOriginalPacket == 0
        packetCompletelyReconstructed(i) = true;
    end
    errRMS_runs(i) = errRMS;
    errMax_runs(i) = errMax;
    finalEncodedPacketErroredBits = finalEncodedPacketErroredBits + totalErroredBitsFromEncodedPacket;
    finalOriginalPacketErroredBits = finalOriginalPacketErroredBits + totalErroredBitsFromOriginalPacket;
    lostBurst = lostBurst + sum(~burstDetected);
    
end

fprintf('\n--- REPORT ---\n'); % title
fprintf('Max error over %d runs; after synchronizations, before equalization: %.3e\n', packets, max(errMax_runs));
fprintf('Mean RMS error over %d runs; after synchronizations, before equalization: %.3e\n', packets, mean(errRMS_runs));
fprintf('STD RMS error over %d runs; after synchronizations, before equalization: %.3e\n', packets, std(errRMS_runs));
fprintf('---------------------\n'); % div
fprintf('Hard errors inferred from LLR sign: %d\n', ...
        finalEncodedPacketErroredBits);
fprintf('Total final errored bits from original packets: %d\n', ...
        finalOriginalPacketErroredBits);
fprintf('---------------------\n'); % div
fprintf('Total bursts erased: %d\n', lostBurst);
fprintf('Total packets completely reconstructed: %d\n', sum(packetCompletelyReconstructed));
fprintf('---------------------\n'); % div
fprintf('Bit Error Rate before decoding %.2f%%\n', 100 * (finalEncodedPacketErroredBits/ (packets*numEncodedBits)));
fprintf('Bit Error Rate after decoding %.2f%%\n', 100 * (finalOriginalPacketErroredBits / (packets*packetLength)));
fprintf('Burst Loss Rate %.2f%%\n', 100 * (lostBurst / (packets * numBursts)));
fprintf('Packet Error Rate %.2f%%\n', 100 * (sum(~packetCompletelyReconstructed)/(packets)));
fprintf('------------------------------------------\n'); % report end