%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             RECEIVER                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errMax, ...
          errRMS, ...
          totalErroredBitsFromEncodedPacket, ... 
          totalErroredBitsFromOriginalPacket, ...
          burstDetected, ...
          numPilotSymbols, ...
          rxEqualizedBursts, ...
          H_est_bursts, ... 
          vettore, ...             % diagnostic
          detections, ...          % diagnostic
          rxMatchedByChannel, ...  % diagnostic
          intraBurstDiagnostic, ... % diagnostic
          reportDiagnosticDDSCE, ...
          thresholdDDSCEBER ... % for DD-SCE report
         ] = receiver( rxSignal, ... % input signal
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
                                    relativeCarrierFrequencies, ...    
                                    rolloff, ...
                                    filterSpan, ...
                                    symbolRate, ...
                                    samplesPerSymbol, ...
                                    Fs, ...
                                    show_metrics_plots, ...
                                    i, ...
                                    numProtectedIndexBits, ...
                                    burstIndexCodebook, ... % Burst-index codebook configuration, RX knows the codebook, but NOT which burst will arrive when.
                                    delay_max, ...
                                    channelEstMargin, ...
                                    friendIndexCodeword ...
                                    )

    %% Conversion in reception of EQUIVALENT-BASEBAND signal with RRC filter 
    %% Realistic Burst Time Synchronization Using Preamble Correlation
    %% Finite States Machine: RF AGC --> ADC --> DIGITAL AGC
    %% Carrier Synchronization (Frequency Synch. and Phase Synch.)
    % The RX knows:
    %   - Set of possible frequency channels
    %   - Synchronization symbols
    %   - PHY waveform parameters
    %
    % The RX DOES NOT know:
    %   - Frequency hopping sequence
    %   - carrier frequency used by each burst
    %   - Burst transmission times
    %   - Burst arrival times
    %   - Exact sample of burst start
    %
    % Therefore burst acquisition is performed jointly in:
    %   1. FREQUENCY: scan all possible carrier frequencies
    %   2. TIME: search the known preamble over the complete received waveform
    
    % Configuration of knownSyncSymbols for Time Synchronization(4.CORRELATION)
    % numbers of complete symbols inside the syncReference
    numSyncSymbols = floor((length(syncReference) + length(phaseReference)) ... 
                                        / bitsPerSymbol);
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
        
    % ADAPTIVE DETECTION THRESHOLD
    % the 3 hyperparameters/thresholds below: 'backgroundWindowLength', 
    % 'guardLength', 'thresholdK', 'minUsableSyncFraction'; are validated 
    % through Monte Carlo campaigns, find the results and the method used 
    % in the Excel file: 'Montecarlo campaigns for thresholds validation'
    % Number of correlation samples used to estimate the local correlation 
    % background.
    backgroundWindowLength = 500; % validated with Monte Carlo
    % Number of correlation samples excluded around the candidate under test.
    % This avoids contaminating the background estimate with the useful peak.
    % if 'guardLength = 0 --> exclude ONLY the CUT (Candidate Under Test)'
    guardLength = numSyncSymbols * ( 32 / 100 ); % validated with Monte Carlo
    % Number of standard deviations above the local correlation background.
    thresholdK = 14; % validated with Monte Carlo
    % Minimum usable fraction of synchronization symbols after blanking.
    % Blanking Compensation variable threshold
    minUsableSyncFraction = 0.3;
    % this parameter sets the minimum fraction of sync symbols not-blanked by  
    % which we can still try to recover the burst
    % (usable_sync_symbols/total_sync_symbols).

    % Positions of the known carrier-pilot symbols inside rxBurstSymbols.
    % rxBurstSymbols has the structure:
    % CP | sync | phase reference | index | payload
    carrierPilotIndicesInBurst = cpLength + (1:numSyncSymbols);

    % TIME-FREQUENCY BURST ACQUISITION
    % Number of possible RF channels known by the receiver.
    % the RX knows ONLY the set of possible frequencies,
    % NOT the hopping sequence.
    numPossibleFrequencyChannels = length(relativeCarrierFrequencies);
    
    % Complete received timeline.
    numRxSamples = length(rxSignal); % whole received packet in samples
    globalSampleIndices = (0:numRxSamples-1).';
    tRx = globalSampleIndices / Fs; % array of time evolution of the packet
    
    % RRC group delay expressed in samples
    % for the blanking mask slot: since the RRC filter introduces a delay 
    % of filterSpan/2​⋅samplesPerSymbol we must delay also the blanking mask
    % when we correlate the candidate so to find the proper instants
    % blanked, so keeping into account the tails for the slot!
    rxRRCGroupDelaySamples = filterSpan * samplesPerSymbol / 2;
    
    % Store the matched-filter output for every possible frequency channel.
    % We use a cell array because later, once a detection is found, we need
    % to return to the corresponding channel waveform and extract the burst.
    % the cell array allows to store data of different dimensions and types
    % together
    rxMatchedByChannel = cell(numPossibleFrequencyChannels,1); % it will be 
    % of length: #numChannels x #numRxSamples
    % Synchronization quantities for each scanned frequency for the
    % 'MetricsPlots'
    syncMetricByChannel = cell(numPossibleFrequencyChannels,1);
    syncThresholdByChannel = cell(numPossibleFrequencyChannels,1);
    
    % Structure containing all possible burst detections found during the
    % time-frequency search.
    detections = struct( ...
        'firstSyncSample', {}, ...
        'channelIndex', {}, ...
        'frequency', {}, ...
        'score', {}, ...
        'metricIndex', {});
    
    % Diagnostic output of the AGC.
    vettore = [];

    %% 1. FREQUENCY SCANNING OVER THE COMPLETE RECEIVED WAVEFORM
    for ch = 1:numPossibleFrequencyChannels
    
        % Candidate relative carrier frequency.
        % The receiver is simply testing one of the possible RF channels.
        fRel = relativeCarrierFrequencies(ch);

        %% 1.1 DOWNCONVERT THE WHOLE RECEIVED WAVEFORM   
        % Downconversion bring the entire time slot back to baseband
        % raw because we still have to apply AGC
        rxChannelBasebandRaw = rxSignal .* exp(-1j * 2*pi * fRel * tRx);
    
        %% 1.2 AGC + ADC    
        % Each frequency channel scan must start from the same receiver state.
        % Do NOT carry the AGC state from one frequency channel to another.
        % since we are simulating a BANK/RAKE demodulator, with 4 parallel 
        % demodulators.
        % FINITE STATE MACHINES, retrieves the external auxiliary function: 
        % 'agc_model' computing: 
        % - 3.1. RF AGC
            % RF AGC is required in order to bring the received signal in an
            % amplitude interval compatible with the front-end of the RX, 
            % before synchronization. It must guarantee that the analog signal 
            % does not overcome the full scale of the ADC.
            % Put the RF AGC after downconversion since it's easier in the
            % equivalent-baseband, because AGC depends on the modulus of the
            % signal, the frequency traslation does not modify amplitude.
        % - 3.2. ADC (CLIPPING + QUANTIZATION) 
        % - 3.3. DIGITAL AGC 
            % AGC DIG is required in order to finalize the RF AGC, over the
            % precise found synchronization symbols.
        [rxChannelBasebandADC, ...
         agcDiagnostic] = ...
            agc_model(rxChannelBasebandRaw, Fs);
    
        % Optional diagnostic storage.
        vettore = [vettore, agcDiagnostic]; %#ok<AGROW> % append with the initialized 
        % one before the for loop over the channels.
    
        %% 1.3 RRC MATCHED FILTER    
        % The same System object is reused for independent channel hypotheses.
        % Reset its internal memory before processing another full waveform.
        reset(rxRRCsync);    
        % The exit still keeps the #samplesPerSymbol samples per symbol 
        % (RRC without decimation)
        rxChannelMatched = rxRRCsync(rxChannelBasebandADC);
    
        % Save it because, if we later detect a burst on this channel,
        % we will need this waveform for burst extraction.
        rxMatchedByChannel{ch} = rxChannelMatched; % store the RRC RX 
        % output
    
    %% 2. TIME SYNCHRONIZATION OVER THE COMPLETE TIMELINE
        numMatchedSamples = length(rxChannelMatched); % total amount of samples
    
        % Earliest possible position of the first synchronization symbol.
        % There must be enough space before it for the cyclic prefix.
        firstCandidate = cpLength * samplesPerSymbol + 1;    
        % Latest possible synchronization position for which the complete
        % burst can still be extracted from the received waveform.
        lastCandidate = numMatchedSamples ...
            - (symbolsPerBurstWithCP - cpLength - 1) * samplesPerSymbol;    
        numCandidates = lastCandidate - firstCandidate + 1; % number of 
        % candidated POSITIONS for starting the preamble (in samples)
    
        % Safety check. (if there are no candidates)
        if numCandidates <= 0
            continue; % skip directly to the next channel
        end
        
        % candidated POSITIONS for starting the preamble (in samples)
        syncMetric = zeros(numCandidates,1);
        % for blanking compensation
        numBlankedSyncSymbols = zeros(numCandidates,1);
        % variable threshold for blanking compensation per candidate
        syncThresholdCandidate = zeros(numCandidates,1);
        % Fraction of usable synchronization symbols for each candidate.
        usableFractionCandidate = zeros(numCandidates,1);
        
        % For each possible sample q, we take:
        % q,
        % q + samplesPerSymbol,
        % q + 2*samplesPerSymbol,
        % ...
        % So to obtain a possible sequence of received symbols
        % We confront it with: knownSyncSymbols.
        for candidateIdx = 1:numCandidates
            % Candidate sample corresponding to the first sync symbol.
            % hence we go after the CP for position correlation.
            q = firstCandidate + candidateIdx - 1;
            % Candidate synchronization-symbol SAMPLING POSITIONS.
            syncSampleIndices = q + (0:numSyncSymbols-1).' * samplesPerSymbol;
            % Extract possible synchronization sequence.
            receivedSyncCandidate = rxChannelMatched(syncSampleIndices);
            % Candidate energy.
            receivedEnergy = sum(abs(receivedSyncCandidate).^2);
            % if the energy is above this threshold we proceed, otherwise no
            if receivedEnergy > 1e-12    
                % Normalize candidate sequence. So, to compare with the 
                % normalized known sequence
                receivedSyncNormalized = ...
                receivedSyncCandidate / sqrt(receivedEnergy);
                % Normalized complex correlation with known preamble.
                correlationValue = ...
                knownSyncSymbolsNormalized' * receivedSyncNormalized; % the 
                % knownSyncSymbolsNormalized' makes the conjugated transpose.
                % Phase-insensitive synchronization metric.
                syncMetric(candidateIdx) = abs(correlationValue)^2; % in order 
                % to cope with phase rotation of the channel, the modulus
                % reimans unchanged while the angle could change
            end
            
            %% BLANKING COMPENSATION
            % rxChannelMatched is referred directly to the complete RX
            % timeline, and is delayed by the RX RRC group delay, so that
            % delay must be removed before checking the blanking state.
            % We only compensate the RRC group delay.
            globalBlankingSampleIndices = syncSampleIndices - rxRRCGroupDelaySamples;    
            % Make sure indices are inside the global blanking mask.
            validMaskIndices = globalBlankingSampleIndices >= 1 & ...
                globalBlankingSampleIndices <= length(blankingMask);    
            % Default assumption outside modeled waveform:
            % Receiver available. --> ones
            candidateBlankingStates = ones(numSyncSymbols,1);    
            % Read the actual local-TX blanking state from the GLOBAL mask.
            % This is the key step for determining which of the 32 candidate sync symbols
            % are actually blanked. candidateBlankingStates is initially set to all ones.
            % Only the candidates that satisfy the blanking-mask boundaries, identified by
            % validMaskIndices, are updated. Their values are set to the corresponding
            % entries of the global blankingMask, using globalBlankingSampleIndices to map
            % each candidate sync symbol to its position in the global mask.
            candidateBlankingStates(validMaskIndices) = ...
               blankingMask(globalBlankingSampleIndices(validMaskIndices));    
            % Number of blanked synchronization symbols.
            numBlankedSyncSymbols(candidateIdx) = ...
                sum(candidateBlankingStates == 0);    
            % Fraction of usable synchronization symbols, for blanking
            % compensation
            numUsableSyncSymbols = ...
                numSyncSymbols - numBlankedSyncSymbols(candidateIdx); 
            % varibale threshold for blanking compensation
            usableFraction = numUsableSyncSymbols / numSyncSymbols;
            % for lowering the detection threshold according to the number
            % of symbols blanked.
            usableFractionCandidate(candidateIdx) = usableFraction;
            % important check-point if too many symbols are blanked
            if usableFraction < minUsableSyncFraction    
                % Reject candidate, candidate cannot be used because too 
                % much of the synchronization sequence has been blanked.
                syncThresholdCandidate(candidateIdx) = Inf; % with 'inf' the 
                % candidate is not valid
            %{    
            else % calculate the variable threshold of Blanking Compens
                % Blanking-compensated threshold.
                syncThresholdCandidate(candidateIdx) = ...
                    syncDetectionThreshold * usableFraction;   
            %}
            end
        end      

        % ADAPTIVE CORRELATION THRESHOLD
        % CA-CFAR Cell-Averaging Constant False Alarm Rate
        halfWindow = floor(backgroundWindowLength / 2);
        halfGuard  = floor(guardLength / 2);        
        for candidateIdx = 1:numCandidates        
            % Skip candidates already rejected because of excessive blanking.
            if isinf(syncThresholdCandidate(candidateIdx))
                continue;
            end        
            % Limits of the local correlation window.
            windowStart = max(1, candidateIdx - halfWindow);
            windowEnd   = min(numCandidates, candidateIdx + halfWindow);        
            % Indices belonging to the local window.
            backgroundIndices = (windowStart:windowEnd).';        
            % Guard region around the candidate under test.
            guardStart = max(1, candidateIdx - halfGuard);
            guardEnd   = min(numCandidates, candidateIdx + halfGuard);        
            % Remove the guard region from the background samples.
            backgroundIndices( ...
                backgroundIndices >= guardStart & ...
                backgroundIndices <= guardEnd) = []; % remove elements from vector     
            % Safety check.
            if isempty(backgroundIndices)        
                syncThresholdCandidate(candidateIdx) = Inf;
                continue;        
            end        
            % Correlation samples representing the local background.
            backgroundMetric = syncMetric(backgroundIndices);
            % Estimate local correlation background.
            muBackground = mean(backgroundMetric);        
            % Estimate local correlation-background dispersion.
            sigmaBackground = std(backgroundMetric);
            % Adaptive detection threshold based on local correlation background.
            backgroundThreshold = muBackground + thresholdK * sigmaBackground;            
            % Blanking compensation.
            usableFraction = usableFractionCandidate(candidateIdx); 
            % Adaptive detection threshold blanknig compensated
            syncThresholdCandidate(candidateIdx) = backgroundThreshold * usableFraction;     
        end
        
    %% 3. NORMALIZED DETECTION SCORE
        % Since every candidate has its own blanking-compensated threshold,
        % compare candidates through the normalized detection score:
        %       detectionScore(q) = syncMetric(q) / threshold(q)
        % Hence:
        %       score > 1  --> candidate passes its threshold
        %       score = 1  --> candidate exactly at threshold
        %       score < 1  --> candidate fails detection    
        detectionScore = zeros(numCandidates,1);  
        % extract only the valid candidates (not infinite)
        validCandidates = isfinite(syncThresholdCandidate);    
        detectionScore(validCandidates) = ...
            syncMetric(validCandidates) ./ ... % element-wise, since it is directly with the arrays
            syncThresholdCandidate(validCandidates);

        % for the 'metricsPlots'
        syncMetricByChannel{ch} = syncMetric; % detectionScore;
        syncThresholdByChannel{ch} = syncThresholdCandidate; % Blanking-compensated CA-CFAR

        
    %% 4. FIND ALL PREAMBLE PEAKS ON THIS FREQUENCY CHANNEL
        
        % NOTE! this would make sense ONLY if we do not care about other
        % transmissions! outside the link TX-RX; but since there are also the
        % FRIENDS (firenldy interferers) we must temporarly erase this part!
        % Two detections belonging to the same physical burst must not be
        % interpreted as two separate bursts.
        % A burst cannot realistically start again a few samples after itself.
        % Use approximately 80% of one burst duration as minimum separation.
       %  minPeakDistanceSamples = ...
       %     max(1, round(0.8 * symbolsPerBurstWithCP * samplesPerSymbol));
    
        % we find all the possible peaks that overcome the threshold in the 
        % current channel, for sure they could be more than 28 even just 
        % for A channel! But then in section: 'BURST INDEX DECODING AND 
        % REORDERING', in the '8-PSK DEMODULATION' block, we will keep only 
        % the bursts with the proper indexing numbers.
        [peakScores, peakCandidateIndices] = ...
            findpeaks(detectionScore, ...
                      'MinPeakHeight', 1 ... % score >= 1 --> threshold overcome     
                      ); 
                % 'MinPeakDistance', minPeakDistanceSamples ... % avoid interpreting different close samples of the same maximum as different bursts
                

        % store the detection, the ordering now does not make any sense NOW
        for p = 1:length(peakCandidateIndices) % how many peaks were found
            candidateIdx = peakCandidateIndices(p);    
            % Convert metric index into matched-filter sample index.
            firstSyncSample = firstCandidate + candidateIdx - 1; 
            newDetection.firstSyncSample = firstSyncSample;    
            newDetection.channelIndex = ch;    
            newDetection.frequency = fRel;    
            newDetection.score = peakScores(p);    
            newDetection.metricIndex = candidateIdx;    
            detections(end+1) = newDetection; %#ok<AGROW>
        end

    end

%{ 
    % NOTE! this would make sense ONLY if we do not care about other
    % transmissions! outside the link TX-RX; but since there are also the
    % FRIENDS (firenldy interferers) we must temporarly erase this part!
    %% 5. REMOVE DUPLICATE TIME DETECTIONS ACROSS FREQUENCY CHANNELS
    % due to leakage, CFO, non-ideal filters...
    % basically: same time, different frequencies --> keep the strongest
    % correlation
    if ~isempty(detections) % if there are some detections   
        % First process strongest detections.
      % [sortedValues, indices], hence we only take their indices ranging in: [1, #detections]
        [~, scoreOrder] = sort([detections.score], 'descend'); % ordered from Max to Min   
        % indeed the structure of the detections ordered by score 
        detectionsByScore = detections(scoreOrder); 
        % create an empty struct array with the same structure (fields) of 
        % 'detections'  
        detectionsClean = detections([]);  
        % Two detections within the interval below are considered the same 
        % burst.
        duplicateToleranceSamples = 2 * samplesPerSymbol;

        for d = 1:length(detectionsByScore)    
            currentDetection = detectionsByScore(d);    
            isDuplicate = false;    
            % at first external iteration is empty, the it expands
            for a = 1:length(detectionsClean)     
                if abs( ...
                    currentDetection.firstSyncSample ...
                    - detectionsClean(a).firstSyncSample) ...
                    <= duplicateToleranceSamples

                    isDuplicate = true;
                    % it does not make any sense to compare with other 
                    % detections, we have already found out that is a 
                    % duplicate, hence --> exit the for cycle
                    break; 
                end
            end
    
            % Since detections are inspected from strongest to weakest,
            % the strongest frequency hypothesis is automatically retained.
            if ~isDuplicate
                % add the detection at the end of the strcut array ('end' 
                % returns the index of the last element)
                detectionsClean(end+1) = currentDetection; %#ok<AGROW>
            end    
        end    
        detections = detectionsClean;  % replace the 'detections' structure  
    end
%}

    %% 6. ORDER DETECTED BURSTS BY ARRIVAL TIME
    if ~isempty(detections)
      % [values, indices]
        [~, timeOrder] = sort([detections.firstSyncSample], 'ascend'); % from min to max
        detections = detections(timeOrder); % ordered in time
    end 

    %% 7. PROCESS ALL DETECTED BURSTS
    % this cycle iterates over all detected bursts
    % NOTE! The receiver does NOT know how many bursts are present.
    % The number of acquisition candidates is determined exclusively
    % by the synchronization detector.
    numDetectedBursts = length(detections);

    % where to store the syncohronized bursts
    rxSCFDEBursts = zeros(symbolsPerBurstWithCP, numDetectedBursts);
    % Storage for the Carrier Frequency Synchronization Estimates
    estimatedCarrierPhaseOffset_rad = nan(numDetectedBursts, 1);
    estimatedCarrierFrequencyOffset_Hz = nan(numDetectedBursts, 1);
    burstDetectedAcquisition = false(numDetectedBursts,1);

    for k = 1:numDetectedBursts
        % k now means:
        % "k-th burst detected in chronological order"
        currentDetection = detections(k);
        % Frequency channel autonomously discovered by the RX.
        detectedChannel = currentDetection.channelIndex;    
        % detectedFrequency = currentDetection.frequency;    
        % First synchronization symbol position.
        firstSyncSample = currentDetection.firstSyncSample;    
        % Retrieve matched-filter output corresponding to the frequency
        % autonomously selected during acquisition.
        rxChannelMatched = rxMatchedByChannel{detectedChannel};    
        numMatchedSamples = length(rxChannelMatched);

        %% 8. FIND CYCLIC PREFIX START    
        % The preamble is the first element after the CP
        firstCPSample = firstSyncSample - cpLength * samplesPerSymbol;
        % now we know both:
            %   - temporal position;
            %   - correct sampling phase.
            % We sample a value each #samplePerSymbols samples
        burstSampleIndices = firstCPSample + (0:symbolsPerBurstWithCP-1).' ...
            * samplesPerSymbol;
        if burstSampleIndices(1) < 1 || ... % || (or)
           burstSampleIndices(end) > numMatchedSamples    
            % Detection exists, but complete burst cannot be extracted.
            burstDetectedAcquisition(k) = false; % burst marked as erased
            continue; % skip to the next burst of the for cycle
        end
        
        % Burst Erasure,If the condition above is overcome --> the burst is not 
        % erased
        burstDetectedAcquisition(k) = true;
    
        rxBurstSymbols = rxChannelMatched(burstSampleIndices);
    
    
    %% 9. FEED-FORWARD CARRIER SYNCHRONIZATION (ACQUISITION)
    % NOTE!: in literature there are plenty of articles concerning the CFO
    % estimation in SC-FDE and its importance, e.g.: 
    % "Are SC-FDE systems robust to CFO?, by: Ying Chen et al."
    % Furthermore in my architecture, it is of core importantce, since a
    % bad CFO estimation might have a catastrphic domino effect for the
    % subsequent architecture blocks, since they are in sequence (channel 
    % estimation, equalization)!
        %% 9.1. EXTRACT THE KNOWN CARRIER PILOTS
        receivedCarrierPilotSymbols = ... % take only the known pilots!
            rxBurstSymbols(carrierPilotIndicesInBurst);
        
        %% 9.2. REMOVE THE KNOWN PILOT MODULATION
        % keep in mind the error we must correct: ϕ(t) = ϕ0 + 2πΔft
        % If:
        % receivedPilot(m) ≈ knownPilot(m) * exp(j*(2πΔft+ϕ0​))
        % then:
        % receivedPilot(m) * conj(knownPilot(m)) ≈ exp(j*(2πΔft+ϕ0​))
        % because 8-PSK symbols have unit magnitude!
        carrierPilotProducts = ...
            receivedCarrierPilotSymbols .* conj(knownSyncSymbols); % so
        % NOTE! carrierPilotProducts contains only the complex exponential: 
        % exp(j*(2πΔft + ϕ0​)) only in the case of a flat channel! In case of
        % multipath, like ours, it also contains: current pilots, other pilots
        % and symbols, ISI, transient of RRC --> hence the phase we feed to the
        % polyfit below also contains: linear phase due to CFO, phase/distortion
        % due to multipath, noise. Neighbouring symbols contribute to each
        % received symbol
        % hence it works as long as: the LOS dominates, other taps are small,
        % pilots are long, CFO generates a slope enough evident, timing is
        % accurate, SNR is high, channel is constant through the burst
        
        %% 9.3. MEASURE AND UNWRAP THE PILOT PHASE
        wrappedPilotPhase = angle(carrierPilotProducts);
        % angle() returns values between -pi and +pi.
        % unwrap() removes artificial jumps of approximately 2*pi, converting
        % the phase evolution into an approximately continuous line.
        unwrappedPilotPhase = unwrap(wrappedPilotPhase);
        
        %% 9.4. ESTIMATE PHASE INTERCEPT AND PHASE SLOPE
        % The first known pilot has relative symbol index zero.
        carrierPilotRelativeIndex = (0:numSyncSymbols-1).';
        % Linear least-squares fit:
        % keep in mind the error we must correct: ϕ(t) = 2πΔft + ϕ0​
        % phase(m) ≈ slope*m + intercept
        % p(1) = slope in radians per symbol
        % p(2) = intercept in radians
        % polyfit(x,y,u): 
        % x --> indipendent variable;
        % y --> dependent variable;
        % u is the degree of the polynome (1 is the linear)
        % then it returns the slope and the intercept
        phaseFit = polyfit(carrierPilotRelativeIndex, unwrappedPilotPhase, 1);
        estimatedPhaseSlope_radPerSymbol = phaseFit(1);
        estimatedPhaseIntercept_rad = phaseFit(2);
        % ϕ(t) = 2πΔft + ϕ0 ...that for us it tranlsates to:
        % unwrappedPilotPhase = ...
        % estimatedPhaseSlope_radPerSymbol * carrierPilotRelativeIndex ...
        % + estimatedPhaseIntercept_rad                
    
        %% 9.5. CONVERT PHASE SLOPE TO CFO
        % phase_variation_per_symbol = [2π * ​Δf] / Symbol_Rate
        % thus: Δf = [phase_variation_per_symbol * Symbol_Rate] / 2*π​
        estimatedCarrierFrequencyOffset_Hz(k) = ...
            estimatedPhaseSlope_radPerSymbol * symbolRate / (2*pi);
        % Wrap the initial phase estimate to the interval [-pi, +pi].
        estimatedCarrierPhaseOffset_rad(k) = ...
            angle(exp(1j*estimatedPhaseIntercept_rad));                
        
        %% 9.6. CORRECT CP and useful burst symbols
        % Symbol index relative to the first synchronization/pilot symbol:
        % CP symbols: -cpLength, ..., -1
        % First synchronization: 0
        % Last useful symbol: symbolsPerBurst-1
        burstRelativeSymbolIndex = (-cpLength:symbolsPerBurst-1).'; % it 
        % start from negative since it also includes the CP, but the
        % reference point (0) coincides with the first sample of the
        % sync reference symbols (as we did in TX). Negative values are
        % inserted since the angular phase rotation allows them, it is
        % an elegant consequence of the conversions: linear<->angular
        assert( ...
            length(burstRelativeSymbolIndex) == length(rxBurstSymbols), ...
       'Carrier-correction vector and received burst have different lengths.');
        % estimate of ϕ(t) = ϕ0 + 2πΔft
        estimatedCarrierPhaseEvolution = ...
            estimatedCarrierPhaseOffset_rad(k) ...
            + 2*pi ...
            * estimatedCarrierFrequencyOffset_Hz(k) ...
            * burstRelativeSymbolIndex ...
            / symbolRate;
        % correct the burst symbols with the opposite (-1j) of the impairments 
        % that were introduced by the channel
        rxBurstSymbolsCarrierCorrected = ...
                rxBurstSymbols .* exp(-1j*estimatedCarrierPhaseEvolution);
    
        % finally at the very end now save the corrected received burst
        rxSCFDEBursts(:,k) = rxBurstSymbolsCarrierCorrected;

        % Matche Filter and Carrier Estimation PLOTS
        if show_metrics_plots == true && i == 1 && k == 1

            % Synchronization quantities corresponding exactly to the
            % frequency channel on which this burst was detected.
            detectedSyncMetric = syncMetricByChannel{detectedChannel};        
            detectedSyncThreshold = syncThresholdByChannel{detectedChannel};        
            % Index of this detection inside the synchronization metric.
            detectedPeakIndex = currentDetection.metricIndex;
        
            MetricPlots.blindMatchedFilter_plot(k, ...
                                                detectedSyncMetric, ...
                                                detectedPeakIndex, ...
                                                detectedSyncThreshold);
        
            MetricPlots.blindFrequencyOffset_plot(k, ...
                                            phaseFit, ...
                                            carrierPilotRelativeIndex, ...
                                            unwrappedPilotPhase);
        end

    end
    

    %% PILOT-BASED Channel Estimation using Sync + Phase reference bits and 
    %% SC-FDE Demodulation with MMSE Equalization
    % NOTE! In literature there are specific sequences that take advantage 
    % of the autocorrelation. They try to set columns of A such that
    % correspondent delayed versions of the training are as much as
    % possible othogonal. In other words: a delayed copy of the training
    % must be very different from the original training.
    % Usuallly in order to find them we could impost a typical optimization
    % problem of Operational Research, with the training sequences and the
    % estimator. --> read the script: "generateOptimizedPreamble.m"

    % NOTE! this is a standard symbol-spaced channel estimation and 
    % equalization. Not fractionally-spaced.
    % Known pilot bits available in each burst
    knownPilotBits = [syncReference; phaseReference];
    numPilotSymbols = floor(length(knownPilotBits)/bitsPerSymbol);
    knownPilotSymbols = modulatedBursts(1:numPilotSymbols, 1);
    numPilotSymbols = length(knownPilotSymbols); % over-write 
    % 'numPilotSymbols' to make sure of its length
    
    % NOTE! the pilot tones (symbols) must be >> number of taps (Lest). 
    % there are two possible ways to operate:
    % 1. Use dedicated Pilot Bursts;
    % 2. Augment the length of the sync and reference bits, at least covering a
    % total of 16-32 symbols.
    % by doing so you will notice a strong decrease in numErrors8PSK even with
    % higher AWGN.
    % 'Lest' stands for --> Length ESTimate
    symbolPeriod = 1 / symbolRate; % duration in seconds of a symbol
    Lest = ceil(delay_max/symbolPeriod) + 1 + channelEstMargin; % length of the 
    % discrete impulsive response we want to estimate. It is the number of 
    % coefficients/taps of the DISCRETE EQUIVALENT CHANNEL the RX assumes 
    % that must estimate. It represents #Lest taps at symbols level, hence 
    % we estimate the channel memory length in symbols, hence it is totally 
    % different from what we did in the Channel implementation before the 
    % RRC RX. Indeed we are doing this after the RRC RX
    
    rxEqualizedBurstsDetected = zeros(symbolsPerBurst, numDetectedBursts); % equalized bursts
    H_est_burstsDetected = zeros(symbolsPerBurst, numDetectedBursts); % channel estimation in frequency
    h_est_burstsDetected = zeros(Lest, numDetectedBursts); % channel estimation in time
    % estimated noise variance after MMSE equalization, for 8-PSK LLR demod
    noiseVarAfterEq = zeros(numDetectedBursts,1);
    for k = 1:numDetectedBursts % cycle over all bursts
        % if the current burst was not found (erased), we directly skip its
        % channel estimation/MMSE
        if ~burstDetectedAcquisition(k)
            continue;
        end
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
        % Pilot residual after channel reconstruction, to estimate White Noise,
        % by taking advantage of the pilots, and just before the equalizer
        % indeed! for the MMSE denominator we need: (σ_n)^2/(σ_x)^2
        pilotResidual = y - A*h_hat;
        % Estimate effective noise + unmodelled interference variance
        % If: (σ_n)^2 --> 0, MMSE converges to ZF
        noiseVarHat = mean(abs(pilotResidual).^2); % if the signal has mean 
        % value equals to zero (E(X)=0) --> Power_X = variance(X); since:
        % Power_X = variance(X) + |mean(X)|^2. Generally speaking AWGN is
        % modelled as: N(0, sigma_n^2), hence it has mean value = 0.
        noiseVarHat = max(noiseVarHat, 1e-12); % floor for safety
        % Symbol power (again for MMSE)
        Px = mean(abs(knownPilotSymbols).^2); % since the constellation we 
        % have used has unitary energy --> Px(8-PSK_symbol) = 1;
        h_est_burstsDetected(:,k) = h_hat; % saves the times response estimates by 
        % burst k
        % Convert time-domain channel estimate to frequency domain
        h_hat_padded = [h_hat; zeros(symbolsPerBurst-Lest,1)]; % padding for FFT, we want 
        % a total of #symbolsPerBurst frequency bins in the FFT, the reason
        % is explained below in the row of the: 'fft(rxBurst, N)'
        H_hat = fft(h_hat_padded, symbolsPerBurst); % channel frequency response
        H_est_burstsDetected(:,k) = H_hat; % saves the frequency response estimates by 
        % burst k
        % SC-FDE Equalization
        Y = fft(rxBurst, symbolsPerBurst); % for equalization frequency-by-frequency we 
        % need a channel estimate with the same identical frequency bins.
        % MMSE equalizaer
        W = conj(H_hat) ./ (abs(H_hat).^2 + noiseVarHat/Px);
        Xhat = W .* Y; % equalization
        % Back to time/symbol domain
        rxBurstEqualized = ifft(Xhat, symbolsPerBurst);
        
        % POST-MMSE residual gain estimation using pilots (for potential 
        % ampitude gain and phase rotations on the burst)
        % Extract equalized pilots
        rxPilotEqualized = rxBurstEqualized(1:numPilotSymbols);    
        % Estimate residual complex gain:
        % rxPilotEqualized ~= gHat * knownPilotSymbols + residual
        % Least Sqaures estimate of gHat, (notice we cannot simplify since 
        % are vectors and matrices operations)
        gHat = (knownPilotSymbols' * rxPilotEqualized) / ...
               (knownPilotSymbols' * knownPilotSymbols);    
        % Numerical safety if too low gain
        if abs(gHat) < 1e-8
            gHat = 1;
        end
        % Normalize the whole burst    
        rxBurstEqualized = rxBurstEqualized / gHat;
        % Store normalized equalized burst
        rxEqualizedBurstsDetected(:,k) = rxBurstEqualized;   
        % Effective post-MMSE variance estimation    
        % Equalized and normalized pilots
        rxPilotEqualized = rxBurstEqualized(1:numPilotSymbols);    
        % Effective residual actually seen by the soft demapper (8-PSK demod)
        % variance from ground-truth!
        postEqResidual = rxPilotEqualized - knownPilotSymbols;    
        % Effective post-equalizatoin error variance (is a power since it is 
        % the square modulus of a complex number representing symbols)
        noiseVarAfterEq(k) = ...
        sum(abs(postEqResidual).^2) / max(numPilotSymbols-1, 1);  % Np-1: 
        % estimating the complex residual gain gHat from the same pilots
        % consumes one complex degree of freedom (equivalently, two real 
        % DoFs). The variance is evaluated after normalization by gHat, 
        % i.e. on the same symbol scale seen by the soft demapper.
        % Numerical floor, to avoid a too low variance that could explose LLR
        noiseVarAfterEq(k) = max(noiseVarAfterEq(k), 1e-8);
        
        % NOTE! based on the 'noiseVarAfterEq', (we can inspect its values
        % from the diagnostic structure below: 'intraBurstDiagnostic'), we
        % can see that usually if the total number of detections made by 
        % the CA-CFAR contains some false alarms, the number of false
        % alarms is usually equal to the number of outliers in the 
        % 'noiseVarAfterEq' array. --> we can implement a burst erasure
        % mechanism here! But only after having properly calibrated the
        % hyperparameters of the CA-CFAR algorithm:
        % ('backgroundWindowLength', 'guardLength', 'thresholdK').
        % use the median value for this, as long as the number of false
        % detections is less than half the total detections made
    
    end
    
    % Metric plot outputted only for the first iteration
    if show_metrics_plots == true && i == 1
    MetricPlots.constellation_plots(rxSCFDEBursts, ...
                                    rxEqualizedBurstsDetected, ...
                                    cpLength, ...
                                    symbolsPerBurst, ...
                                    numPilotSymbols, ...
                                    burstDetectedAcquisition);
    end

    
    %% 8-PSK Demodulation
    % take the equalized bursts
    rxModulatedBursts = rxEqualizedBurstsDetected;
    % Temporary storage indexed by DETECTION order.
    % At this stage column k does NOT necessarily correspond to the 
    % chronologically transmitted burst k. The true burst position will be 
    % recovered from 'burstIndexBits'.
    rxFramedBurstsPaddedLLR_detected = ...
        zeros(symbolsPerBurst * bitsPerSymbol, numDetectedBursts);
    
    % Parameters for decision-directed intra-burst erasure
    % The two hyperparameters/thresholds below are validated through Monte 
    % Carlo campaigns, find the results and the method used in the Excel 
    % file: 'Montecarlo campaigns for thresholds validation'
    localWindowSymbols = 12; % validated with Monte Carlo
    gammaErasure = 70;   % validated with Monte Carlo
    HeaderLength = length(syncReference) + ...
                          length(phaseReference) + ...
                          numProtectedIndexBits;
    % First LLR belonging to coded payload
    firstPayloadLLR = HeaderLength + 1;
    % Decision Directed Intra-burst erasure DIAGNOSTIC storage
    intraBurstDiagnostic = struct('noiseVarAfterEq', [], ...
                                  'sigmaDD', [], ...
                                  'reliabilityRatio', [], ...
                                  'erasedWindow', [], ...
                                  'firstSymbol', [], ...
                                  'lastSymbol', [], ...
                                  'firstPayloadLLR', [], ... % first LLR position
                                  'lastPayloadLLR', [], ... % last LLR position
                                  'berGroundTruth', [], ... % for DD-SCE BER diagnostic
                                  'numBitErrors', [], ...   % for DD-SCE BER diagnostic
                                  'numBitsWindow', []);     % for DD-SCE BER diagnostic
    
    % 8-PSK demodulation cycle
    for k = 1:numDetectedBursts
        % if the current burst was not found (erased), we directly skip its
        % demodulation and we go to the next burst.
        if ~burstDetectedAcquisition(k)
            % Erased burst --> leave all LLRs equal to zero, initializated in:
            % 'rxFramedBurstsPaddedLLR_detected'; sync failure⇒burst erasure⇒LLR=0⇒
            % ⇒Turbo decoder receives no information from that burst.
            continue;
        end
        rxSymbols_k = rxModulatedBursts(:,k);
        % estimated noise variance at the MMSE equalizer output
        noiseVariance_k = noiseVarAfterEq(k); 
        % MATLAB pskdemod LLR convention:
        %   positive LLR -> bit 0 more likely
        %   negative LLR -> bit 1 more likely
        burstLLR = pskdemod(rxSymbols_k, ...
                            M, ...
                            pi/M, ...
                            'gray', ...
                            'OutputType','approxllr', ... % instead of 'llr' is simpler and more robust, it uses max-log approximation instead of the sum of the likelihoods
                            'NoiseVariance', noiseVariance_k);
        
        % DECISION-DIRECTED
        % Hard 8-PSK decisions
        decidedIndices = pskdemod(rxSymbols_k, ...
                                  M, ...
                                  pi/M, ...
                                  'gray');
        % Reconstruct ideal constellation symbols
        decidedSymbols = pskmod(decidedIndices, ...
                                M, ...
                                pi/M, ...
                                'gray');    
        % Symbol-wise post-equalization DD residual
        ddResidual = rxSymbols_k - decidedSymbols;        
        numSymbols = length(rxSymbols_k); % #symbolsPerBurst        
        % DIAGNOSTIC
        % exclude the preamble for the residual calculus
        firstPayloadSymbol = floor(HeaderLength / bitsPerSymbol) + 1;
        numPayloadSymbolsDD = max(0, numSymbols - firstPayloadSymbol + 1);
        numLocalWindows = ceil(numPayloadSymbolsDD / localWindowSymbols);
        % Balanced partition of payload symbols across DD-SCE windows
        baseWindowSymbols = floor(numPayloadSymbolsDD / numLocalWindows);
        extraSymbols = mod(numPayloadSymbolsDD, numLocalWindows);        
        % The first extraSymbols windows contain one additional symbol
        windowSizes = baseWindowSymbols * ones(numLocalWindows, 1);
        windowSizes(1:extraSymbols) = windowSizes(1:extraSymbols) + 1;

        sigmaDD_windows = nan(numLocalWindows,1);
        reliabilityRatio_windows = nan(numLocalWindows,1);
        erasedWindow = false(numLocalWindows,1);        
        firstSymbol_windows = nan(numLocalWindows,1);
        lastSymbol_windows = nan(numLocalWindows,1);        
        firstPayloadLLR_windows = nan(numLocalWindows,1);
        lastPayloadLLR_windows = nan(numLocalWindows,1);        
        firstSymbol = firstPayloadSymbol;
        % DECISION-DIRECTED INTRABURST ERASURE
        % refers to a signal processing and communication technique where a 
        % receiver uses its own internal hard or soft decisions (its "best 
        % guesses" of the transmitted data symbols) as a substitute for 
        % known training sequences or pilot symbols to adapt filters, track 
        % channel states, or correct errors.
        % for cycle over the INTRA-BURST WINDOWS
        for windowIdx = 1:numLocalWindows
            lastSymbol = firstSymbol + windowSizes(windowIdx) - 1;
            symbolIndices = firstSymbol:lastSymbol;  
            % Local DD residual power
            sigmaDD_m = mean(abs(ddResidual(symbolIndices)).^2);    
            % Normalize against the GLOBAL post-EQ pilot estimate
            reliabilityRatio = sigmaDD_m / max(noiseVariance_k, 1e-12);            
            % LLR range corresponding to the current SYMBOL window
            firstLLR_window = (firstSymbol - 1)*bitsPerSymbol + 1;
            lastLLR_window  = lastSymbol*bitsPerSymbol;            
            % Select ONLY the portion belonging to the coded payload.
            % This protects:
            %   - sync reference bits
            %   - phase reference bits
            %   - protected burst-index bits
            firstLLR_toErase = max(firstLLR_window, firstPayloadLLR);
            lastLLR_toErase  = lastLLR_window;            
            isErased = false; % for diagnostic
            % if reliabilityRatio > gammaErasure && firstLLR_toErase <= lastLLR_toErase            
            %    burstLLR(firstLLR_toErase:lastLLR_toErase) = 0;
            %    isErased = true; % for diagnostic
            % end   
            %{%}
            idx = firstLLR_toErase:lastLLR_toErase;
            if ~isempty(idx)
                if reliabilityRatio > gammaErasure && firstLLR_toErase <= lastLLR_toErase
                    burstLLR(idx) = 0;                                   % hopeless window: erase
                    isErased = true;
                elseif reliabilityRatio > 1 && firstLLR_toErase <= lastLLR_toErase
                    burstLLR(idx) = burstLLR(idx) / reliabilityRatio;    % noisier window: lower confidence
                end
            end
            
            % for diagnostic
            sigmaDD_windows(windowIdx) = sigmaDD_m;
            reliabilityRatio_windows(windowIdx) = reliabilityRatio;
            erasedWindow(windowIdx) = isErased;        
            firstSymbol_windows(windowIdx) = firstSymbol;
            lastSymbol_windows(windowIdx) = lastSymbol;            
            % for diagnostic purposes
            if firstLLR_toErase <= lastLLR_toErase
                firstPayloadLLR_windows(windowIdx) = firstLLR_toErase;
                lastPayloadLLR_windows(windowIdx) = lastLLR_toErase;
            end

             % Starting symbol of the next balanced window
            firstSymbol = lastSymbol + 1;
        end
        
        % Diagnostic storage
        intraBurstDiagnostic(k).noiseVarAfterEq = noiseVariance_k;
        intraBurstDiagnostic(k).sigmaDD = sigmaDD_windows;
        intraBurstDiagnostic(k).reliabilityRatio = reliabilityRatio_windows;
        intraBurstDiagnostic(k).erasedWindow = erasedWindow;        
        intraBurstDiagnostic(k).firstSymbol = firstSymbol_windows;
        intraBurstDiagnostic(k).lastSymbol = lastSymbol_windows;        
        intraBurstDiagnostic(k).firstPayloadLLR = firstPayloadLLR_windows;
        intraBurstDiagnostic(k).lastPayloadLLR = lastPayloadLLR_windows;

        % Store according to DETECTION order for the moment.
        rxFramedBurstsPaddedLLR_detected(:,k) = burstLLR; % a soft value per each bit
    end


    %% BURST INDEX DECODING AND REORDERING
    % Length of the known part preceding the burst identifier.
    knownHeaderLength = length(syncReference) + length(phaseReference);    
    % Positions of all CODEBOOK-CODED burst-index bits.
    indexBitPositions = knownHeaderLength + (1:numProtectedIndexBits);
    % Final LLR matrix.
    % IMPORTANT:
    % columns of this matrix will correspond to the ACTUAL chronologically 
    % transmitted burst identifiers, not to detection order.
    rxFramedBurstsPaddedLLR = zeros(symbolsPerBurst * bitsPerSymbol, numBursts);
    % Final burst detection vector indexed by the true burst identifier.
    % it is simply the chronologically reordered version
    burstDetectedReordered = false(numBursts,1);

    headerLength = length(syncReference) + ...
                   length(phaseReference) + ...
                   numProtectedIndexBits;

    % Since the metrics will compare 'rxEqualizedBursts' with 'modulatedBursts'
    % we should reorder some of them with the following 3 arrays that will
    % substitute their original related ones once reordered after the following cycle
    rxEqualizedBurstsReordered = zeros(symbolsPerBurst, numBursts);   % for reordered metrics
    rxSCFDEBurstsReordered = zeros(symbolsPerBurstWithCP, numBursts); % for reordered metrics
    H_est_burstsReordered = zeros(symbolsPerBurst, numBursts);        % for reordered metrics
    % DD-SCE diagnostic reordered by transmitted burst index, structure
    % for storing each one of the 28 bursts. Leave it like this, it will be
    % overwritten (so far we have the first detection repeated #numBursts
    % times).
    % the reason is to take advantage of the struct fields already written!
    intraBurstDiagnosticReordered = repmat(intraBurstDiagnostic(1), numBursts, 1);

    % whether at low SNR two detections declare the same ID, we take 
    % advantage of the LLR reliability:
    bestIndexReliability = -inf(numBursts,1);
    % main cycle for reordering the bursts
    for detectionIdx = 1:numDetectedBursts
        % Skip missing detections.
        if ~burstDetectedAcquisition(detectionIdx)
            continue; % skip to the following detection
            % Erased burst --> leave all LLRs equal to zero
        end    
        % LLRs associated with this detected burst.
        currentBurstLLR = rxFramedBurstsPaddedLLR_detected(:,detectionIdx);    
        % Extract the LLRs associated with the coded burst identifier.
        protectedIndexLLR = currentBurstLLR(indexBitPositions);
        % Soft metric for every possible burst identifier.
        % MATLAB convention:
        %   positive LLR -> bit 0 more likely
        %   negative LLR -> bit 1 more likely
        % Therefore:
        %   code bit 0 -> +LLR contribution
        %   code bit 1 -> -LLR contribution
        indexScores = zeros(numBursts,1);
        for candidateIdx = 1:numBursts
            candidateCodeword = burstIndexCodebook(candidateIdx,:).';        
            bipolarCodeword = 1 - 2*candidateCodeword;     
            % comparison with the all the codewords and the index extracted 
            % from the burst
            indexScores(candidateIdx) = sum(bipolarCodeword .* protectedIndexLLR);        
        end        

        % FRIEND index hypothesis
        % Compute the FRIEND score using exactly the same soft metric
        % used for all desired burst identifiers.
        friendCodeword = friendIndexCodeword(:); % transpose it to column vector
        friendBipolarCodeword = 1 - 2*friendCodeword;
        friendIndexScore = sum(friendBipolarCodeword .* protectedIndexLLR);
        % Joint DESIRED / FRIEND decision
        % Best desired hypothesis.
        [bestDesiredScore, decodedBurstIndex] = max(indexScores);
        % Compare the FRIEND hypothesis against the best desired one.
        
        if friendIndexScore > bestDesiredScore
            % ---------------------------------------------------------
            % FRIENDLY INTERFERER
            % ---------------------------------------------------------
            % The acquisition is valid because FRIEND intentionally uses
            % the same preamble as the desired transmitter.
            % However, its protected index identifies it as FRIEND.
            % Therefore it must NOT enter the desired-burst reordering.
            %{
            % DIAGNOSTIC Debug information for desired bursts.
            friendMargin = friendIndexScore - bestDesiredScore;
            fprintf( ...
                ['Detection %d: FRIEND | ' ...
                 'friendScore = %.3f | ' ...
                 'bestDesired = %.3f | ' ...
                 'margin = %.3f\n'], ...
                detectionIdx, ...
                friendIndexScore, ...
                bestDesiredScore, ...
                friendMargin);
            %}
            continue; % skip to the following detection, NOTE! this is very 
            % important, because it means that we do not insert FRIENDS
            % into the arrays for the decoding below, for the final blocks.
        end

        % Reliability of the DESIRED decision
        % The reliability must include FRIEND among the competing
        % hypotheses.
        % Since bestDesiredScore is the winning hypothesis here:
        % (friendIndexScore < bestDesiredScore)
        % find its strongest competitor among:
        %   - all other desired burst IDs
        %   - FRIEND
        competingDesiredScores = indexScores;
        competingDesiredScores(decodedBurstIndex) = -inf; % set the best 
        % score to -inf, so to extract the second best one in the below line:
        secondBestDesiredScore = max(competingDesiredScores);
        strongestCompetitorScore = max(secondBestDesiredScore, friendIndexScore);
        % 'indexRealiability' is expressed as a MARGIN!
        indexReliability = bestDesiredScore - strongestCompetitorScore;
        
        %{
        % DIAGNOSTIC Debug information for desired bursts.
        fprintf( ...
            ['Detection %d: DESIRED %d | ' ...
             'desiredScore = %.3f | ' ...
             'friendScore = %.3f | ' ...
             'reliability = %.3f\n'], ...
            detectionIdx, ...
            decodedBurstIndex, ...
            bestDesiredScore, ...
            friendIndexScore, ...
            indexReliability);
        %}

        % NOTE! POSSIBLE BURST ERASURE CASE!
        % The hyperparameter/threshold below is validated through Monte 
        % Carlo campaigns, find the results and the method used in the 
        % Excel file: 'Montecarlo campaigns for thresholds validation'
        if indexReliability < 4
            continue;
        end

        % If another detected burst already claimed this identifier,
        % retain the one whose index bits are more reliable.
        if indexReliability <= bestIndexReliability(decodedBurstIndex)        
            continue; % if we are not sure enough we avoid inserting this 
            % detection, and we skip directly to the following one
        end
        bestIndexReliability(decodedBurstIndex) = indexReliability;                

        % Put this received frame in its TRUE transmitted position.
        rxFramedBurstsPaddedLLR(:,decodedBurstIndex) = currentBurstLLR;    
        burstDetectedReordered(decodedBurstIndex) = true;
        rxEqualizedBurstsReordered(:,decodedBurstIndex) = ... 
        rxEqualizedBurstsDetected(:,detectionIdx);  % for reordered metrics
        rxSCFDEBurstsReordered(:,decodedBurstIndex) = ...    
        rxSCFDEBursts(:,detectionIdx); % for reordered metrics        
        H_est_burstsReordered(:,decodedBurstIndex) = ...      
            H_est_burstsDetected(:,detectionIdx); % for reordered metrics
        % Reorder DD-SCE diagnostic
        intraBurstDiagnosticReordered(decodedBurstIndex) = ...
            intraBurstDiagnostic(detectionIdx);
    end
    
    % From this point onward 'burstDetected' is indexed according to the
    % transmitted burst identifier rather than acquisition/detection order.
    burstDetected = burstDetectedReordered; 
    rxEqualizedBursts = rxEqualizedBurstsReordered; % for reordered metrics
    rxSCFDEBursts = rxSCFDEBurstsReordered;         % for reordered metrics
    H_est_bursts = H_est_burstsReordered;           % for reordered metrics
    intraBurstDiagnostic = intraBurstDiagnosticReordered; % for DD-SCE BER



    %% SCE GROUND-TRUTH BER DIAGNOSTIC  
    % NOTE!
    % LIMITATIONS OF THE DD-SCE BER DIAGNOSTIC:
    % 1. Burst association:
    %    BER assumes correct matching between received and transmitted
    %    bursts. Incorrect burst IDs or FRIEND misclassification may
    %    lead to misleading BER values and TP/FP/FN/TN statistics.
    %
    % 2. Hard-decision BER:
    %    BER is computed from hard-bit decisions on equalized symbols,
    %    while the receiver uses soft LLRs. Therefore, this metric
    %    detects windows containing bit errors but does not assess
    %    whether erasing their LLRs benefits the Turbo decoder.
    for k = 1:numBursts    
        % Only consider detected and accepted bursts
        if ~burstDetected(k)
            continue;
        end    
        % True transmitted bits (including header and padding)
        txBits = pskdemod(modulatedBursts(:,k), ...
                          M, pi/M, 'gray', 'OutputType', 'bit');    
        % HARD DECISIONS before SCE erasure
        rxBits = pskdemod(rxEqualizedBursts(:,k), ...
                          M, pi/M, 'gray', 'OutputType', 'bit');    
        % One BER value per SCE window
        % how many windows adopted by DD-SCE
        nWindows = numel(intraBurstDiagnostic(k).erasedWindow);
        % for storing final data
        berWindows = nan(nWindows,1);
        numErrorsWindows = zeros(nWindows,1);
        numBitsWindows = zeros(nWindows,1);    
        % Last valid bit, excluding 8-PSK padding
        lastValidBit = numel(txBits) - numPadBits8PSK;    
        for m = 1:nWindows
            firstBit = intraBurstDiagnostic(k).firstPayloadLLR(m);
            lastBit = intraBurstDiagnostic(k).lastPayloadLLR(m);    
            % Exclude padding
            lastBit = min(lastBit, lastValidBit);    
            if isnan(firstBit) || isnan(lastBit) || firstBit > lastBit
                continue;
            end    
            % Actual bit errors in this window
            nErrors = sum(txBits(firstBit:lastBit) ~= ...
                          rxBits(firstBit:lastBit));    
            nBits = lastBit - firstBit + 1;    
            numErrorsWindows(m) = nErrors;
            numBitsWindows(m) = nBits;
            berWindows(m) = nErrors / nBits;
        end
        % Store diagnostic results
        intraBurstDiagnostic(k).berGroundTruth = berWindows;
        intraBurstDiagnostic(k).numBitErrors = numErrorsWindows;
        intraBurstDiagnostic(k).numBitsWindow = numBitsWindows;    
    end
    % SCE DECISION ACCURACY DIAGNOSTIC    
    TP = 0;
    FP = 0;
    FN = 0;
    TN = 0;    
    for k = 1:numBursts    
        if ~burstDetected(k)
            continue;
        end    
        ber = intraBurstDiagnostic(k).berGroundTruth;
        erased = intraBurstDiagnostic(k).erasedWindow;    
        % Ignore windows without a valid BER
        valid = ~isnan(ber); 
        thresholdDDSCEBER = 0.174;
        degraded = ber > thresholdDDSCEBER; % ber > 0.5 means that the majority of the bits 
        % in the given window are errored. 
        % For a Turbo Coder with rate 1/3, the theoretical limit should be
        % around 17.4% of raw BER for still being able to completely
        % reconstruct the original packet.
        % sum(& & &) --> how many elements satisfy simultaneously all three 
        % conditions
        TP = TP + sum(valid & degraded & erased);
        FP = FP + sum(valid & ~degraded & erased);
        FN = FN + sum(valid & degraded & ~erased);
        TN = TN + sum(valid & ~degraded & ~erased);    
    end    
    % Classification metrics
    precisionSCE = TP / max(TP + FP, 1);
    recallSCE = TP / max(TP + FN, 1);  
    % diagnostic Structure for storing results
    reportDiagnosticDDSCE.TP = TP;
    reportDiagnosticDDSCE.FP = FP;
    reportDiagnosticDDSCE.FN = FN;
    reportDiagnosticDDSCE.TN = TN;
    reportDiagnosticDDSCE.precisionSCE = 100*precisionSCE;
    reportDiagnosticDDSCE.recallSCE = 100*recallSCE;

 

    % Remove 8-PSK padding
    rxFramedBurstsLLR = rxFramedBurstsPaddedLLR(1:end-numPadBits8PSK, :);
    
    
    %% Turbo decoder
    % NOTE! Turbo decoder can correct a lack of information; but in 
    % general, it cannot autonomously correct a wrong LLRs mapping!
    % Recover encoded packet from bursts
    
    %remove headers
    rxBurstsPayloadLLR = rxFramedBurstsLLR(headerLength+1:end,:);
    % stack into column
    % rxEncodedPacketPaddedLLR = reshape(rxBurstsPayloadLLR,[],1);
    rxEncodedPacketPaddedLLR = reshape(rxBurstsPayloadLLR.', [], 1); % de-interleaving

    % leave padding
    rxEncodedPacketLLR = rxEncodedPacketPaddedLLR(1:numEncodedBits);
               
    % numIterations = 6; % the default is 6
    turboDec = comm.TurboDecoder( ... % 'TrellisStructure', trellis, ...
        'InterleaverIndices', interleaverIndices);
        %'NumIterations', numIterations);
    % convention used by turboDecoder MATLAB is with the reversed values (-)
    decodedPacket = turboDec(-rxEncodedPacketLLR); 

    
    %% DIAGNOSTIC CHECKPOINTS
    errMax = max(abs(rxSCFDEBursts(:) - txSCFDEBursts(:)));
    errRMS = sqrt(mean(abs(rxSCFDEBursts(:) - txSCFDEBursts(:)).^2));    
    % just to compare the whole packet from Tx to Rx
    rxEncodedPacketHardFromLLR = rxEncodedPacketLLR < 0;
    % convert LLR sign back to hard decisions
    totalErroredBitsFromEncodedPacket = sum(rxEncodedPacketHardFromLLR ~= encodedPacket);
    totalErroredBitsFromOriginalPacket = sum(packet ~= decodedPacket);   
    %{
    fprintf('Max error after synchronizations, before equalization: %.3e\n', errMax);
    fprintf('RMS error after synchronizations, before equalization: %.3e\n', errRMS);
    fprintf('Hard errors inferred from LLR sign: %d\n', ...
        totalErroredBitsFromEncodedPacket);
    fprintf('Total final errored bits from original packet: %d\n', ...
        totalErroredBitsFromOriginalPacket);
    fprintf('Total bursts erased: %d\n', ...
        sum(~burstDetected));
    %}
%burstDetectedAcquisition
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            END RECEIVER                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%