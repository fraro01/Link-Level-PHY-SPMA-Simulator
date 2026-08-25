
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
                    modulatedBursts(numPilotSymbols+1:end,k); % reference
                receivedPayload  = ...
                    rxEqualizedBursts(numPilotSymbols+1:end,k); % estimate
                % calculate the error
                errorVector = receivedPayload - referencePayload;
                % RMS 
                evmRMS = sqrt( ...
                   sum(abs(errorVector).^2)/sum(abs(referencePayload).^2));        
                % percentage
                evmPerBurst_percent(k) = 100 * evmRMS; % in percentage
            end        
        end
        
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