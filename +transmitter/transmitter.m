%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           TRANSMITTER                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [packet, ... % original packet with random bits
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
          startTimes, ... % per burst
          txSignal, ... % FINAL TRANSMITTED SIGNAL
          numProtectedIndexBits, ...
          burstIndexCodebook, ...
          friendIndexCodeword, ...
          bitsPerBurst ...
         ] = transmitter(...
          packetLength, ...
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
          cpLength, ... % Cyclic Prefix for SCFDE
          tunerDelay_s, ... % tuner delay for carriers
          dutyCycle, ...
          burstIndexCodebook, ...
          friendIndexCodeword ...
          )
 

    %% Originate the Packet
    % Generate random bits (0 or 1)
    packet = randi([0 1], packetLength, 1); %column vector not row for applying 
    % turbo coding    
    
    %% Turbo Coding
    % I could add trellis and numeiterations
    % Constituent convolutional code: rate 1/3
    % trellis = poly2trellis(4, [13 15 17], 13);
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
    % length of the turbo encoded packet
    numEncodedBits = length(encodedPacket);
    % "integer superior part" taken --> exceeding bits --> padding bits
    bitsPerBurst = ceil(numEncodedBits / numBursts);
    numPaddingBits = bitsPerBurst*numBursts - numEncodedBits;
    % padding the encoded packet with #numPaddingBits zeros added at the end
    % the Rx MUST "unpad" the received packet
    encodedPacketPadded = [encodedPacket; zeros(numPaddingBits,1)];
    % (bitsPerBurst x numBursts) matrix keeping the order of bits
    % bursts = reshape(encodedPacketPadded, bitsPerBurst, numBursts);
    bursts = reshape(encodedPacketPadded, numBursts, bitsPerBurst).'; % inter-burst interleaving

    %TODO! ERASEME, and uncomment the row just above
    % INTER-BURST INTERLEAVING INDICES
    % The interleaver is a fixed pseudo-random permutation that must be IDENTICAL
    % at the TX and at the RX, without exchanging any variable between them.
    % This is achieved by generating it from a hard-coded seed:
    %   1. rng(7,'twister') sets the Mersenne Twister generator to seed 7 and
    %      returns in 'prevState' the state the generator had before the call,
    %      i.e. the exact point reached by the simulation stream (seeded in
    %      main.m). It acts as a bookmark.
    %   2. randperm(N) is then deterministic: every call made right after
    %      rng(7,'twister') returns the very same permutation of 1:N, so TX and
    %      RX obtain the same indices.
    %   3. rng(prevState) restores the simulation stream, so that the rest of the
    %      simulation (new packet bits, burst timings, channel, noise,
    %      interferers) proceeds exactly as if these lines had never been
    %      executed. Without this step the random stream would stay stuck on
    %      seed 7 and every packet would see the same "random" values.
    % NOTE! The seed (7), the generator type ('twister') and N must be the same
    % at TX and RX, otherwise the de-interleaving returns scrambled LLRs.
    %{
    % inter-burst interleaving: random permutation with a fixed seed, known
    % to both TX and RX
    prevState = rng(7, 'twister');               % save simulation state, set interleaver seed
    interburstIdx = randperm(numBursts*bitsPerBurst);
    rng(prevState);                              % restore simulation state
    xInt = encodedPacketPadded(interburstIdx);   % encodedPacketPadded is already a column
    bursts = reshape(xInt, bitsPerBurst, numBursts);   % one column = one burst
    %}

    
    % essential to extract the proper length
    numProtectedIndexBits = size(burstIndexCodebook, 2);
    
    % matrix of zeros
    % sync | phase reference | burst index | payload
    framedBursts = zeros(bitsPerBurst + ...
                         length(syncReference) + ...
                         length(phaseReference) + ...
                         numProtectedIndexBits, ... 
                         numBursts);
    for k = 1:numBursts
        % Row k of the codebook represents burst identifier k-1.
        protectedBurstIndexBits = burstIndexCodebook(k,:).';
        % frame structure
        framedBursts(:,k) = [ ... % take all rows of cloumn k vertical concatenation of three vectors
            syncReference;           % 48 bit
            phaseReference;          % 48 bit
            protectedBurstIndexBits; % 15 bit: burst identifier
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
    % symbol rate = 13 MHz; maximum 17 MHz, check if the one derived is
    % compliant to these numbers.
    % parameters
    symbolsPerBurstWithCP = size(txSCFDEBursts,1);
    
    % Configuration of Root Raised Cosine TX Filter, we configure the "rolloff"
    % now, because is used to calculate the Symbol Rate from the formula below
    % timing
    % NOTE! below we also include the variable: 'filterSpan' we get the 
    % real physical time for transmitting a burst, because we take into
    % account also the tails of the RRC
    burstTxTime = (symbolsPerBurstWithCP + filterSpan)/ symbolRate; % seconds
    % Duty-cycle validity check.
    if dutyCycle <= 0 || dutyCycle > 1
        error('dutyCycle must be in the interval (0,1].');
    end
    % Frequency channels
    channelBandwidth = (fMax - fMin) / numFrequencyChannels;
    carrierFrequencies = fMin + channelBandwidth/2 + ...
        (0:numFrequencyChannels-1)*channelBandwidth; % list comprehension (0:i)        
    % array of values between 1 and #numFrequencyChannels of length #numBursts
    frequencyChannels = randi( ...
        [1 numFrequencyChannels], ...
        numBursts, ...
        1);
    carrierForBurst = carrierFrequencies(frequencyChannels).';

    %% Duty Cycle + Tuner Delay + Pseudo-Random Burst Timing
    if dutyCycle == 1

        % start times defintion and tuner delay if two consecutive bursts are 
        % sent on different carriers
        tunerDelayPerBurst = zeros(numBursts,1); % Per-burst tuner delay
        startTimes = zeros(numBursts,1); % start TX cumulative
        for k = 2:numBursts
            if frequencyChannels(k) ~= frequencyChannels(k-1)
                tunerDelayPerBurst(k) = tunerDelay_s;   % 1 us
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
    
    %% Final transmitted EQUIVALENT-BASEBAND signal with Pulse Shape Filter RRC
    % NOTE! Due to Nyquist and thanks to the fact that we are simulating in
    % baseband, we must use: samplesPerSymbol > 2*(fmax_occupied)/symbolRate
    % in which fmax_occupied = [max(relativeCarrierFrequencies) +
    % (relativeCarrierFrequencies(4) - relativeCarrierFrequencies(3))]
    % Atually, with RRC the one-sided band of a single channel is: 
    % Bsingle = symbolRate/2 * (1 + rolloff) --> hence fmax_occupied is less.
    % in any case make sure the samplesPerSymbol is set AT LEAST accordingly
    % NOTE! see: https://en.wikipedia.org/wiki/Raised-cosine_filter for time
    % and frequency domain plots, for a better understanding of RRC. RRC makes
    % symbols extremely overlap over the time, but thanks to the first Nyquist
    % criterion for pulse shaping filters, exactly in the sampling points the
    % interferers are equal to zero.
    
    % sample rate --> how much finely we represent the waveform over the time.
    % symbol rate --> how quick the info changes.
    
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
    % Complete observation interval:
    % includes also the final random listening interval.
    totalTime = endTimes(end) + randomListeningIntervals(end);
    % (last burst)
    tailSamples = filterSpan * samplesPerSymbol; % tail extention in samples 
    numSamplesTotal = ceil(totalTime * Fs) + tailSamples + 1; % +1 for safety 
    % reasons
    txSignal = zeros(numSamplesTotal,1);
    % we move to baseband, to avoid computational overheads of L-band (Nyquist)
    fCenter = mean(carrierFrequencies);
    % Frequencies relative to bandwidth center
    relativeCarrierFrequencies = carrierFrequencies - fCenter;
    
    %{
    The Root Raised Cosine filter is a memory FIR filter because the 
    response associated with each symbol spans multiple symbol intervals. 
    Since each burst in the simulator is pulse-shaped independently, before 
    each burst I reset the filter's internal state using reset (txRRC), 
    preventing the previous burst from influencing the next. Furthermore,
    after the last useful symbols, I add zero symbols to continue advancing 
    the filter and allow the residual transient response to fully 
    dissipate. These zeros do not represent transmitted information  
    symbols, but serve exclusively to flush the filter's tail.

    The precise shape of the initial transient (due to the start from the 
    null state) depends on the first symbols of the burst, while that of 
    the final transient depends on the last symbols. Since these sequences 
    are generally different, the two tails will not be mathematically 
    mirrored.

    When thinking about RRC always differentiate between Macroscopic vision
    --> burst; Microscopic vision --> symbols OR samples;
    Filter_Length = filterSpan*samplesPerSymbol + 1; the '+1' represents
    the center of the filter response.
    %}

    for k = 1:numBursts
        reset(txRRC); % clean the filter memory per each burst, otherwise 
        % MATLAB interprets two consecutive bursts as a continuous 
        % sequence, and not being independent; since the final state of the
        % first burst influences the the start of the second one.
        % RESET resets the state before the start of the new burst.
        burstSymbols = txSCFDEBursts(:,k); % from where we sample
        % Append zero-valued symbols to flush the TX RRC filter state and
        % preserve the complete pulse-shaping tail of the burst. Extra zero 
        % symbols appended to flush the TX/RX pulse-shaping filters. In 
        % order to copmletely exit from the filter the end of the current 
        % burst.
        burstSymbolsPadded = [burstSymbols; zeros(filterSpan,1)];     
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