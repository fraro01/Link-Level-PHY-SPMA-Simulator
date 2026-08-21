%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           TRANSMITTER                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [txSignal, numBursts, carrierForBurst, startTimes, Fs, ...
    symbolsPerBurstWithCP, extraSymbols, samplesPerSymbol, symbolRate, ...
    fMax, fCenter, txSCFDEBursts, burstTxTime, slotDuration, ...
    syncReference, phaseReference, bitsPerSymbol, modulatedBursts, ...
    rolloff, filterSpan, cpLength, frequencyChannels, ...
    relativeCarrierFrequencies, symbolsPerBurst, M, numPadBits8PSK, ...
    numEncodedBits, encodedPacket, interleaverIndices, packet, ...
    dutyCycle] = transmitter()

    %% Originate the Packet
    % Set the length of the packet in bits
    packetLength = 1331;
    % Generate random bits (0 or 1)
    packet = randi([0 1], packetLength, 1); %column vector not row for applying 
    % turbo coding
    
    
    %% Turbo Coding
    % TODO! add trellis and numeiterations
    % NOTE! study the turbo coding functionings
    % Constituent convolutional code: rate 1/3
    trellis = poly2trellis(4, [13 15 17], 13);
    % Interleaver turbo
    interleaverIndices = randperm(packetLength);
    % Turbo encoder rate 1/3, (it is the default)
    % Use as "InterleaverIndices" the one defined above
    turboEnc = comm.TurboEncoder(... %'TrellisStructure', trellis, ...
                    'InterleaverIndices', interleaverIndices);
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
    framedBursts = zeros(bitsPerBurst + length(syncReference) + ...
                            length(phaseReference), ... % + numIndexBits, ... 
                            numBursts);
    for k = 1:numBursts
        % decimal to binary, with numIndexBits as #bits of result, most
        % significant bit is put left (used in TLC), then we transpose it: .'
      % burstIndexBits = de2bi(k-1, numIndexBits, 'left-msb').';
        framedBursts(:,k) = [ % take all rows of cloumn k
        % vertical concatenation of three vectors
            syncReference;        % 48 bit
            phaseReference;       % 48 bit
          % burstIndexBits;       % 5 bit: burst identifier
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
    cpLength = 16; % TODO! Must exceed the maximum channel delay in symbols
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
    % symbol rate = 13 MHz; maximum 17 MHz, check if the one derived is
    % compliant to these numbers.
    % parameters
    numFrequencyChannels = 4;
    symbolsPerBurstWithCP = size(txSCFDEBursts,1);
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
    
    
    %% Final transmitted EQUIVALENT-BASEBAND signal with Pulse Shape Filter RRC
    % NOTE! Due to Nyquist and thanks to the fact that we are simulating in
    % baseband, we must use: samplesPerSymbol > 2*(fmax_occupied)/symbolRate
    % in which fmax_occupied = [max(relativeCarrierFrequencies) +
    % (relativeCarrierFrequencies(4) - relativeCarrierFrequencies(3))]
    % Atually, with RRC the one-sided band of a single channel is: 
    % Bsingle = symbolRate/2 * (1 + rolloff) --> hence fmax_occupied is less.
    % TODO! in any case adjust the samplesPerSymbol accordingly.
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
    % experience ACI, since our channel spacing is roughly 13.83 MHz --> set
    % rolloff accordingly <-- TODO!
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

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          END TRANSMITTER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%