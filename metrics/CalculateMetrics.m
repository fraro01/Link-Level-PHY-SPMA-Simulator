
%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
* 
%}



% related-functions collector ad a Static Class
classdef CalculateMetrics

    methods (Static) % class with static methods, hence a container of 
                     % related functions, without requiring an object,
                     % without properties
        
        function [sinr_linear, ... % SINR linear
                  EbN0_linear, ... % E_b/N_0 linear
                  evmPerBurst_percent ... % EVM percentage
                  ] = calculate_metrics( ...
                  Fs, ...
                  rxSignalImpaired, ...
                  rxSignalAtAntenna, ...
                  rxSignalWithInterference, ...
                  interferer_waveform, ...
                  packetLength, ...
                  numBursts, ...
                  burstDetected, ...
                  modulatedBursts, ...
                  numPilotSymbols, ...
                  rxEqualizedBursts)                     
            %% SINR CALCULUS 
            %% calculated in the 'channel' before the 'receiver'
            % since the SINR metric calculates the interference and noise
            % that affect the desired signal, they must be calculated
            % according to the desired signal, hence to its active mask,
            % used as a common mask.
            activeMask = abs(rxSignalImpaired) > 1e-15;
            noiseWaveform = rxSignalAtAntenna - rxSignalWithInterference;
            signalPower_W = mean(abs(rxSignalImpaired(activeMask)).^2);
            interferencePower_W = mean( ...
                abs(interferer_waveform(activeMask)).^2);
            noisePower_W = mean(abs(noiseWaveform(activeMask)).^2);
            sinr_linear = signalPower_W / (interferencePower_W + ...
                noisePower_W);        
            %% E_b/N_0
            %% calculated in the 'channel' before the 'receiver'
            packetEnergy_J = sum(abs(rxSignalImpaired).^2) / Fs;
            Eb_J = packetEnergy_J / packetLength; % divided per INFO BITS!
            noisePower_W = mean(abs(noiseWaveform).^2);
            N0_W_per_Hz = noisePower_W / Fs;
            EbN0_linear = Eb_J / N0_W_per_Hz;        
            %% ERROR VECTOR MAGNITUDE CALCULATED ON THE PAYLOAD
            %% calculated in the 'receiver' after SC-FDE, before 8PSK demod
            evmPerBurst_percent = nan(numBursts,1);
            % for cycle over iterated over all bursts
            for k = 1:numBursts
                % if the current burst was not detected, set the values to 
                % NaN
                if ~burstDetected(k)
                    evmPerBurst_percent(k) = NaN;
                    continue;
                end
                % take the two measures        
                referencePayload = ...
                    modulatedBursts(numPilotSymbols+1:end,k); % reference, 
                % they are the symbols that directly exit from the 8PSK mod
                % excluded the preamble, since it is the "dirtiest" part in
                % reception, used for better estimating the following info
                % bits
                receivedPayload  = ...                    
                    rxEqualizedBursts(numPilotSymbols+1:end,k); % estimated
                % symbols, cleaned with equalization, right before the 8PSK
                % demodulation.
                % calculate the error
                errorVector = receivedPayload - referencePayload;
                % RMS 
                evmRMS = sqrt( ...
                   sum(abs(errorVector).^2)/sum(abs(referencePayload).^2));        
                % percentage
                evmPerBurst_percent(k) = 100 * evmRMS; % in percentage
            end        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % CA-CFAR FALSE DETECTIONS STATISTICS
        function result = calculate_false_detections( ...
                                                detections, ...
                                                desiredGroundTruth, ...
                                                friendGroundTruth, ...
                                                fMin, ...
                                                fMax, ...
                                                numFrequencyChannels, ...
                                                Fs, ...
                                                filterSpan, ...
                                                cpLength, ...
                                                samplesPerSymbol, ...
                                                minCpLength, ...
                                                max_interferer_delay)
        
            % 1. PARAMETERS
            max_symbol_delay = max(minCpLength, max_interferer_delay);
            % Initial temporal matching tolerance (in symbols)
            toleranceSymbols = 3;      
            assert(toleranceSymbols >= max_symbol_delay, ...
            ['tolerance symbols too short: toleranceSymbols = %d symbols, ' ...
             'but at least %d symbols are required by the maximum channel ' ...
             'delay of both TX and interferers.'], ...
            toleranceSymbols, max_symbol_delay);
            toleranceSamples = toleranceSymbols * samplesPerSymbol;        
            % TX + RX RRC group delay + cyclic prefix
            syncOffsetSamples = (filterSpan + cpLength) * samplesPerSymbol;
        
            % 2. BUILD COMPLETE GROUND TRUTH        
            % Combine DESIRED and FRIEND bursts
            allStartTimes = [ ...
                desiredGroundTruth.startTimes(:); ...
                friendGroundTruth.startTimes(:)];        
            allCarrierFrequencies = [ ...
                desiredGroundTruth.carrierForBurst(:); ...
                friendGroundTruth.carrierForBurst(:)];        
            numDesired = numel(desiredGroundTruth.startTimes); % 'numel' counts the number of elements
            numPhysicalBursts = numel(allStartTimes);        
            % 1 = DESIRED, 2 = FRIEND
            burstType = [... # vertical concatenation of the following two arrays
                ones(numDesired,1); ... % array of ones of equal length of the desired
                2*ones(numPhysicalBursts-numDesired,1)]; % array of twos of equal length of the friends
        
            % 3. EXPECTED SYNCHRONIZATION POSITIONS
            expectedSyncSamples = ...
                round(allStartTimes * Fs) + 1 + syncOffsetSamples;        
            % Convert absolute carrier frequencies into channel indices
            channelBandwidth = ...
                (fMax - fMin) / numFrequencyChannels;        
            carrierFrequencies = ...
                fMin + channelBandwidth/2 + ...
                (0:numFrequencyChannels-1)*channelBandwidth;        
            expectedChannels = zeros(numPhysicalBursts,1);
            for b = 1:numPhysicalBursts   
               %[value,   position     ]
                [~, expectedChannels(b)] = min( ...
                    abs(carrierFrequencies - allCarrierFrequencies(b)));
            end
        
            % 4. BUILD POSSIBLE DETECTION-BURST ASSOCIATIONS        
            numDetections = numel(detections);
            % Columns:
            % [detectionIndex, physicalBurstIndex, timingError]
            candidateMatches = zeros(0,3);
            for d = 1:numDetections        
                detectionSample = detections(d).firstSyncSample;
                detectionChannel = detections(d).channelIndex;        
                for b = 1:numPhysicalBursts        
                    % Frequency channel must match, if not --> skip to
                    % following detection
                    if detectionChannel ~= expectedChannels(b)
                        continue;
                    end        
                    timingError = abs( ...
                        detectionSample - expectedSyncSamples(b));     
                    % if they are close enough --> mark it as a detection
                    if timingError <= toleranceSamples 
                        % add a new row with the current detection
                        candidateMatches(end+1,:) = [ ...
                            d, b, timingError]; %#ok<AGROW>        
                    end
                end
            end
        
            % 5. ONE-TO-ONE MATCHING   
            % NOTE! it is a matching greedy search: if there are close
            % bursts in the same channel, a detection could be classified
            % as a false, even if it could have been true for another burst.
            % each physical burst can be associated only to one detection
            % and viceversa.
            % Prefer associations with the smallest timing error
            if ~isempty(candidateMatches)        
                candidateMatches = sortrows(candidateMatches,3);        
            end        
            detectionMatched = false(numDetections,1);
            burstMatched = false(numPhysicalBursts,1);
            trueDesired = 0;
            trueFriend = 0;        
            for m = 1:size(candidateMatches,1) % cycle over the number of rows       
                d = candidateMatches(m,1); % detectionIndex
                b = candidateMatches(m,2); % physicalBurstIndex
        
                % A detection and a physical burst can be used only once!!!
                % IF AT LEAST ONE OF THEM IS TRUE --> SKIP TO THE NEXT ROW
                if detectionMatched(d) || burstMatched(b) 
                    continue;
                end
                detectionMatched(d) = true;
                burstMatched(b) = true;        
                if burstType(b) == 1
                    trueDesired = trueDesired + 1;
                else
                    trueFriend = trueFriend + 1;
                end
            end
        
            % 6. FINAL METRICS
            result.totalDetections = numDetections;
            result.trueDetectionsDesired = trueDesired;
            result.trueDetectionsFriend = trueFriend;        
            result.trueDetections = trueDesired + trueFriend;        
            result.falseDetections = numDetections - result.trueDetections;        
            % False Discovery Proportion
            if numDetections > 0        
                result.FDP = result.falseDetections / numDetections;        
            else        
                result.FDP = 0;        
            end
        
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % For displaying values in the final Report
        function metrics_calculated = final_values(metrics, ...
                numEncodedBits, packetLength, numBursts)
            % taking advantage of the MATLAB 'struct'
            metrics_calculated.errMax = max(metrics.errMax_runs);
            metrics_calculated.errRMS = mean(metrics.errRMS_runs);
            metrics_calculated.errSTD_RMS = std(metrics.errRMS_runs);
            metrics_calculated.finalEncodedBitErrors = ...
                metrics.finalEncodedPacketErroredBits;
            metrics_calculated.finalOriginalBitErrors = ...
                metrics.finalOriginalPacketErroredBits;
            metrics_calculated.SINR = 10*log10(mean(metrics.SINR_runs));
            metrics_calculated.EVM_pct = mean(metrics.EVM_runs, 'all', ...
                'omitnan');
            metrics_calculated.EVM_dB = 20*log10(mean(metrics.EVM_runs, ...
                'all', 'omitnan')/100);
            metrics_calculated.EbN0_linear = ...
                mean(metrics.EbN0_linear_runs);
            metrics_calculated.EbN0_dB = ...
                10*log10(mean(metrics.EbN0_linear_runs));
            metrics_calculated.lostBurst = metrics.lostBurst;
            metrics_calculated.ReconstructedPacket = ...
                sum(metrics.packetCompletelyReconstructed);
            metrics_calculated.BER_before_dec = 100 * ...
                (metrics.finalEncodedPacketErroredBits/ ...
                (metrics.packets*numEncodedBits));
            metrics_calculated.BER_after_dec = 100 * ...
                (metrics.finalOriginalPacketErroredBits / ...
                (metrics.packets*packetLength));
            metrics_calculated.BLR = 100 * (metrics.lostBurst / ...
                (metrics.packets * numBursts));
            metrics_calculated.PER = 100 * ...
                (sum(~metrics.packetCompletelyReconstructed) / ...
                (metrics.packets));
        end

    end % methods end
end