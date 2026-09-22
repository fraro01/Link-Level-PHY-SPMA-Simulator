
%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
*
%}



% related-functions collector ad a Static Class
classdef MetricPlots

    methods (Static) % class with static methods, hence a container of 
                     % related functions, without requiring an object,
                     % without properties

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Frame Structure PLOT
        function frame_plot(bitsPerBurst, ...
                            cpLength, ...
                            syncReference, ...
                            phaseReference, ...
                            numProtectedIndexBits, ...
                            M)            
            % Plot complete binary frame with illustrative Cyclic Prefix
            % Visualization convention:
            % 0 -> -1
            % 1 -> +1
            % 1. Build the different frame sections
            % Illustrative CP: copy the last #cpLength information bits
            randomBurst = randi([0 1], bitsPerBurst, 1);
            cpBits = randomBurst(end-(log2(M)*cpLength)+1:end);
            % Preamble
            preambleBits = [syncReference; ...
                            phaseReference];
            % randomly generated indexing
            random_index = randi([0 1], numProtectedIndexBits, 1);
            % Complete frame for visualization
            frameBits = [cpBits; ...
                         preambleBits; ...
                         random_index; ...
                         randomBurst];
            % Convert:
            % 0 -> -1
            % 1 -> +1
            framePlot = 2*frameBits - 1;
            % Bit indices
            n = 0:length(framePlot)-1; % starts from 0
            % 2. Define the three regions
            numCPBits       = length(cpBits);
            numPreambleBits = length(preambleBits);
            idxCP = 1:numCPBits;
            idxPreamble = numCPBits+1 : numCPBits+numPreambleBits;
            idxIndex = numCPBits+numPreambleBits+1 : numCPBits+numPreambleBits+numProtectedIndexBits;
            idxPayload = numCPBits+numPreambleBits+numProtectedIndexBits+1 : length(framePlot);
            % 3. Plot
            figure;
            hold on; % draw multiple things on the same plot
            % Cyclic Prefix
            % 'stem()' allows to visualize discrete data, furthermore 
            % automatically sets distinct colors
            stem(n(idxCP), ... % 'n' starts from 0
                 framePlot(idxCP), ...
                 'filled', ... % markers at the top of the sticks are filled instead of "empty circles"
                 'LineWidth', 1.2);
            % Preamble
            stem(n(idxPreamble), ...
                 framePlot(idxPreamble), ...
                 'filled', ...
                 'LineWidth', 1);
            % Index
            stem(n(idxIndex), ...
                 framePlot(idxIndex), ...
                 'filled', ...
                 'LineWidth', 1);
            % Information packet
            stem(n(idxPayload), ...
                 framePlot(idxPayload), ...
                 'filled', ...
                 'LineWidth', 0.7);
            % 4. Region boundaries
            xline(numCPBits - 0.5, ... % position in the x axes, '-0.5' sets the line exactly in between the two semantically different bits
                  '--', ... %  '--' dashed line
                  'CP / Preamble', ...
                  'LineWidth', 1.2);
            xline(numCPBits + numPreambleBits - 0.5, ...
                  '--', ...
                  'Preamble / Index', ...
                  'LineWidth', 1.2);
            xline(numCPBits + numPreambleBits + numProtectedIndexBits - 0.5, ...
                  '--', ...
                  'Index / Payload', ...
                  'LineWidth', 1.2);
            % 5. Figure formatting
            grid on;
            ylim([-1.2 1.2]); % visualization boundaries
            yticks([-1 0 1]); % sets the visualization ticks of the y axes
            xlim([-25 length(framePlot)+25]); % visualization boundary of x 
            % axes
            xlabel('Bit index n');
            ylabel('Binary value');
            title('Binary Burst Frame Structure');
            legend('Cyclic Prefix', ...
                   'Preamble', ...
                   'Burst Index', ...
                   'Information bits', ...
                   'Location', 'best'); % set the location of the legend 
            % autonomously in the best position
            hold off;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Channel Delay Taps PLOT
        function channelTaps_plot(pathDelays_s, ...
                                  pathGains)               
            % 'pathGains' is a matrix: #numberOfSamples X #numberOfPaths
            % Convert physical path delays to microseconds
            pathDelays_us = pathDelays_s * 1e6; % multiplied by 1e6, for seconds
            % Choose one time instant of the channel
            % Using approximately the sample in the middle of the simulated 
            % waveform, hence centre of the channel
            timeSampleToPlot = round(size(pathGains,1)/2); % take the number of 
            % rows of pathGains, hence the number of samples and then take its
            % middle, round it to the closest integer
            % Complex channel coefficients at that instant
            instantaneousTaps = pathGains(timeSampleToPlot,:); % take all columns 
            % of row 'timeSampleToPlot'    
            % Tap magnitudes
            instantaneousTapMagnitude = abs(instantaneousTaps); % since it is a 
            % complex value  
            % Instantaneous Channel Taps
            figure;    
            stem(pathDelays_us, ...
                 instantaneousTapMagnitude, ...
                 'filled', ...
                 'LineWidth', 1.3, ...
                 'MarkerSize', 7);    
            grid on; % activate the background grid on the plot
            xlabel('Excess delay [\mus]');
            ylabel('Tap magnitude |h_l|');    
            title(sprintf( ...
                'Instantaneous Multipath Channel Taps - Sample %d', ...
                timeSampleToPlot));    
            xlim([min(pathDelays_us)-0.1 ...
                  max(pathDelays_us)+0.1]);
            
            % Power Delay Profile    
            % Average power of each path over the complete simulated waveform
            pdpLinear = mean(abs(pathGains).^2, 1); % take the average over the first dimension, hence the rows --> time
            % Normalize with respect to strongest path
            pdpNormalized = pdpLinear / max(pdpLinear);    
            figure;    
            stem(pathDelays_us, ...
                 pdpNormalized, ...
                 'filled', ...
                 'LineWidth', 1.3, ...
                 'MarkerSize', 7);    
            grid on;    
            xlabel('Excess delay [\mus]');
            ylabel('Normalized average power');    
            title('Power Delay Profile (PDP) of Channel Taps');    
            ylim([0 1.1]);    
            xlim([min(pathDelays_us)-0.1 ...
                  max(pathDelays_us)+0.1]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% True vs Estimated Channel PLOT
        function channel_estimate(H_est_bursts, ...
                                  startTimes, ...
                                  Fs, ...
                                  numPilotSymbols, ...
                                  samplesPerSymbol, ...
                                  pathGains, ...
                                  symbolsPerBurst, ...
                                  burstDetected)
            if burstDetected(1) % do it only for the first burst IF is det  
                % 1. Estimated channel in frequency domain
                Hest = H_est_bursts(:, 1);    
                % 2. True channel
                % Take the channel coefficients
                burstStartSample = round(startTimes(1)*Fs) + 1;
                pilotStartSample = burstStartSample;
                pilotEndSample = pilotStartSample + ...
                    numPilotSymbols*samplesPerSymbol - 1;
                % Safety
                pilotStartSample = max(1,pilotStartSample);
                pilotEndSample = min(pilotEndSample,size(pathGains,1));                
                % Average complex channel coefficient of each path
                hTrue = mean( ... % mean over the rows (time)
                    pathGains(pilotStartSample:pilotEndSample,:), 1).';                
                % Frequency response
                hTruePadded = [hTrue; ... % 0-padding for FFT
                               zeros(symbolsPerBurst-length(hTrue),1)];                
                Htrue = fft(hTruePadded,symbolsPerBurst);        
                % 3. Normalize both responses    
                HtrueNorm = abs(Htrue) / max(abs(Htrue));    
                HestNorm = abs(Hest) / max(abs(Hest));    
                % 4. FFT-bin axis    
                fftBin = 0:symbolsPerBurst-1;    
                % 5. Plot
                figure;    
                plot(fftBin, ...
                     HtrueNorm, ...
                     'LineWidth', 1.5);    
                hold on;    
                plot(fftBin, ...
                     HestNorm, ...
                     '--', ...
                     'LineWidth', 1.5);    
                grid on;    
                xlabel('FFT bin index');
                ylabel('Normalized magnitude');    
                title(sprintf( ...
                 'Physical Multipath Reference vs Symbol Rate LS-Estimated Channel Response of First Burst'));    
                legend('Physical channel reference', ...
                       'LS channel estimate', ...
                       'Location', 'best');    
                ylim([0 1.1]);    
                hold off;    
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Preamble Matched-Filter Correlation Metric PLOT
        function blindMatchedFilter_plot(detectionIdx, ...
                                         syncMetric, ...
                                         peakIndex, ...
                                         syncThreshold)
            figure;
        
            % Small window around detected correlation peak
            halfWindow = 250;        
            firstIdx = max(1, peakIndex - halfWindow);
            lastIdx  = min(length(syncMetric), peakIndex + halfWindow);
        
            plot(firstIdx:lastIdx, ...
                 syncMetric(firstIdx:lastIdx), ...
                 'LineWidth', 1.2);
        
            hold on;
        
            plot(firstIdx:lastIdx, ...
                 syncThreshold(firstIdx:lastIdx), ...
                 '--', ...
                 'LineWidth', 1.2);
        
            plot(peakIndex, ...
                 syncMetric(peakIndex), ...
                 'o', ...
                 'MarkerSize', 8, ...
                 'LineWidth', 1.5);
        
            grid on;
        
            xlabel('Matched-filter sample index');
            ylabel('Synchronization metric');
        
            title(sprintf('Blind Synchronization - Detection %d', ...
                detectionIdx));

            legend('Synchronization metric', ...
                    'Blanking-compensated threshold', ...
                    'Detected peak', ...
                    'Location', 'best');
        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Frequency Offset Estimation PLOT
        function blindFrequencyOffset_plot(detectionIdx, ...
                                   phaseFit, ...
                                   carrierPilotRelativeIndex, ...
                                   unwrappedPilotPhase)

            figure;
        
            plot(carrierPilotRelativeIndex, ...
                 unwrappedPilotPhase, ...
                 'o', ...
                 'MarkerSize', 6);
        
            hold on;
        
            fittedPhase = polyval( ...
                phaseFit, ...
                carrierPilotRelativeIndex);
        
            plot(carrierPilotRelativeIndex, ...
                 fittedPhase, ...
                 'LineWidth', 1.5);
        
            grid on;
        
            xlabel('Pilot symbol index');
            ylabel('Unwrapped phase [rad]');
        
            legend( ...
                'Measured pilot phase', ...
                'Linear phase fit', ...
                'Location', 'best');
        
            title(sprintf( ...
                'Carrier Frequency Offset Estimation - Detection %d', ...
                detectionIdx));
        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Constellation Before and After Equalization PLOT
        function constellation_plots(rxSCFDEBursts, ...
                                     rxEqualizedBursts, ...
                                     cpLength, ...
                                     symbolsPerBurst, ...
                                     numPilotSymbols, ...
                                     burstDetected)
            if burstDetected(1)
                % BEFORE EQUALIZATION
                rxBeforeEq = rxSCFDEBursts(cpLength+1 : ...
                            cpLength+symbolsPerBurst, 1);    
                % Remove known pilot portion
                rxBeforeEqPayload = rxBeforeEq(numPilotSymbols+1 : end);
                scatterplot(rxBeforeEqPayload);    
                title(sprintf(['8-PSK Constellation BEFORE Equalization of' ...
                            ' First Burst']));        
                % AFTER EQUALIZATION   
                rxAfterEqPayload = rxEqualizedBursts(numPilotSymbols+1:end, 1);    
                scatterplot(rxAfterEqPayload);    
                title(sprintf(['8-PSK Constellation AFTER MMSE ' ...
                    'Equalization of First Burst']));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end