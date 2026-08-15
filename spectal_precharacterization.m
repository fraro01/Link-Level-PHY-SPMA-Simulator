
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    SPECTRAL PRE-CHARACTERIZATION                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use an Equivalent Baseband representation of the SINGLE waveform and
% exclude the Time-Frequency Hopping. Since the RF traslation does not 
% change the waveform spectrum shape, and the time-frequency hopping scheme
% would only replicate the waveform N times around the chosen N carriers.

% IMPORTANT: in order to properly execute this script, at first, run the
% MATLAB link-level simulator, so to retrieve its workspace.

clc;
% whether you want to compare multiple figure/plot
%close all;



%% Symbol rate for a nominal 13 MHz channel
%{%}
% rolloff = 1; % [0;1]
% filterSpan = 15; % in symbols
channelBandwidth = 13e6;
symbolRate = channelBandwidth / (1 + rolloff);
Fs = symbolRate * samplesPerSymbol;


%% 1. Ideal baseband waveform
% same RRC filter
txRRCspectral = comm.RaisedCosineTransmitFilter( ...
    'Shape','Square root', ...
    'RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',filterSpan, ...
    'OutputSamplesPerSymbol',samplesPerSymbol);
% container for wich each cell can  contain data of different types and 
idealShapedBursts = cell(numBursts,1); % dimensions
for k = 1:numBursts % iterate over all 28 bursts
    reset(txRRCspectral); % clean the filter memory per each burst
    burstSymbols = [ ...
        txSCFDEBursts(:,k); % take the burst after the SC-FDE has added CP
        zeros(filterSpan,1) % Add zeros to flush the TX/RX filter memory response tail
    ];
    % retrieve the cell
    idealShapedBursts{k} = txRRCspectral(burstSymbols); % pulse shaping
end
% concatenate vertically the content
txWaveformIdeal = vertcat(idealShapedBursts{:});
% NOTE! INTERESTING TO CHECK!
% txWaveformIdeal = interferer_waveform;
% txWaveformIdeal = blankingMask;

%% 2. Power Spectral Density 
% we use PWELCH to estimate PSD since FFT could be quite noisy due to 
% variance; PWELCH divides the long signal into smaller blocks partially 
% overlapped, each block is multiplied by the Hannning window so to avoid 
% border distortions, we calculate the FFT over each block and take the 
% final average of all results per frequencies power over all blocks.
% segment length used by pwelch() to estimate PSD, Welch window cannot be 
% longer than the signal itself.
segmentLength = min(4096, length(txWaveformIdeal));
[psdLinear, frequencyAxis_Hz] = pwelch( ...
            txWaveformIdeal, ... % INPUT SIGNAL to calculate PSD
            hann(segmentLength,'periodic'), ... % Hann window to calculate PSD to smooth the side lobes avoids spectral leakages, periodic since we use it for FFT
            floor(segmentLength/2), ... % overlap between two consecutive window segments (50%) in samples
            16384, ... % number of FFT points
            Fs, ... % sampling frequency 
            'centered', ... % center the spectrum in 0 Hz since we are in baseband
            'psd'); % we get as output the PSD
% fundamental compromise of PWELCH():
% * frequency resolution: longer segments --> narrower peaks and closer
% frequencies more distinguishable;
% * estimate variance: shorter segments --> more segments to average -->
% smoother and more stable PSD.
% normalize the PSD with respect to the maximum (carrier) and convert to
% dBc
psdNormalized_dBc = 10*log10(psdLinear/max(psdLinear));

%% 3. Occupied Bandwidth: 99% and 99.9%
% cumsum() --> cumulative sum for an array
cumulativePowerNormalized = cumsum(psdLinear) / sum(psdLinear); % is
% taken normalized, since we then need to apply percentages for: 99-99.9%
lowerIndex99  = find(cumulativePowerNormalized >= 0.005,  1);
upperIndex99  = find(cumulativePowerNormalized >= 0.995,  1);
lowerIndex999 = find(cumulativePowerNormalized >= 0.0005, 1);
upperIndex999 = find(cumulativePowerNormalized >= 0.9995, 1);
occupiedBandwidth99_Hz = frequencyAxis_Hz(upperIndex99) ...
    - frequencyAxis_Hz(lowerIndex99);
occupiedBandwidth999_Hz = frequencyAxis_Hz(upperIndex999) ...
    - frequencyAxis_Hz(lowerIndex999);

%% 4. Spectrum width at -40, -50 and -60 dBc
psdSmoothed = movmean(psdLinear,25); % moving average to PSD so to smooth 
% the curve so to avoid small oscillations due to PSD estimate (25 window 
% length)
% calculate again the dBc scale, but on the smoothed PSD version 
psdSmoothed_dBc = 10*log10(psdSmoothed/max(psdSmoothed));
% chosen reference levels
referenceLevels_dBc = [-40 -50 -60];
% vector to store results
spectralWidths_Hz = zeros(size(referenceLevels_dBc));
% minimum of the absolute values of the frequencies
[~,centerIndex] = min(abs(frequencyAxis_Hz)); % ~ ignores the first value, 
% instead take just the position
for k = 1:length(referenceLevels_dBc)
    level = referenceLevels_dBc(k);
    % find the last sample on the left side of the spectrum that is still
    % under the threshold
    leftIndex = find( ... % returns the indices of only True values
        psdSmoothed_dBc(1:centerIndex) <= level, ... % array of True/False
        1, ... % I need just one value
        'last'); % in the specific return only the last index that satisfies the condition
    % likewise, but this time returns only the FIRST one that is under the 
    % threshold
    rightIndex = centerIndex - 1 + find( ... % -1, for scaling on the original vector based on the center found
        psdSmoothed_dBc(centerIndex:end) <= level, ...
        1, ...
        'first');
    % extrcat the width
    spectralWidths_Hz(k) = frequencyAxis_Hz(rightIndex) ...
                - frequencyAxis_Hz(leftIndex);
end

%% 5. Peak-to-Average Power Ratio
instantaneousPower = abs(txWaveformIdeal).^2;
% peak-average ratio in dB
papr_dB = 10*log10(max(instantaneousPower) / mean(instantaneousPower));

%% 6. Different Throughputs Calculations
% there are two distinct throughput measures (with or without Duty Cycle), 
% but always as Gross measures, since the link-level (PHY) does not 
% differentiate bits.
% throughput calculated: symbolRate x bitsPerSymbol (without Duty Cycle)
groos_throughput_without_DC = symbolRate*log2(M); 
% gross thoughput as: numTotalBits / TotalTime (with duty cycle)
gross_throughput_with_DC = symbolsPerBurstWithCP*bitsPerSymbol / ...
            ((symbolsPerBurstWithCP/symbolRate)/dutyCycle);

%% Results
fprintf('\n--- SPECTRAL PRE-CHARACTERIZATION ---\n'); % title
fprintf('Roll-off: %.0f%%, Filter Span: %.0f symbols\n', rolloff*100, ... 
                        filterSpan); % fundamental parameters used
fprintf('Symbol rate: %.3f Msymbol/s\n', symbolRate/1e6); % symbolRate used
fprintf('OBW 99%%: %.3f MHz\n', occupiedBandwidth99_Hz/1e6); % OBW 99%
fprintf('OBW 99.9%%: %.3f MHz\n', occupiedBandwidth999_Hz/1e6); % OBW 99.9%
for k = 1:length(referenceLevels_dBc) % spectrum width at desired levels
    fprintf( ...
        'Width at %d dBc: %.3f MHz\n', ...
        referenceLevels_dBc(k), ...
        spectralWidths_Hz(k)/1e6);
end
fprintf('PAPR: %.3f dB\n', papr_dB); % PAPR
fprintf('---------------------'); % div
fprintf('\n--- List of Different Throughputs ---\n'); % sub-title
fprintf(['Gross throughput without Duty Cycle ' ...
    '(symbolRate X bitsPerSymbol): %.3f Mb/s\n'], ...
            groos_throughput_without_DC/1e6);
fprintf('Gross throughput with Duty Cycle: %.3f Mb/s\n', ...
                                             gross_throughput_with_DC/1e6);
fprintf('------------------------------------------\n'); % report end

%% Single summary plot
figure; % open the graphical window
plot(frequencyAxis_Hz/1e6, ...
     psdSmoothed_dBc, ...
    'LineWidth',1.2);
grid on; % activate the grid on the graph
hold on; % allow to write over the already present plot without erasing it
% Nominal MHz channel limits
xline(-channelBandwidth/2/1e6, '--', ...
    sprintf('-%.3f MHz', channelBandwidth/2/1e6));
xline(channelBandwidth/2/1e6, '--', ...
    sprintf('%.3f MHz', channelBandwidth/2/1e6));
% Occupied bandwidth limits 
xline(frequencyAxis_Hz(lowerIndex99)/1e6,'-.','OBW 99%');
xline(frequencyAxis_Hz(upperIndex99)/1e6,'-.');
xline(frequencyAxis_Hz(lowerIndex999)/1e6,':','OBW 99.9%');
xline(frequencyAxis_Hz(upperIndex999)/1e6,':');
% Reference levels
yline(-40,'--','-40 dBc');
yline(-50,'--','-50 dBc');
yline(-60,'--','-60 dBc');
xlabel('Frequency relative to carrier [MHz]');
ylabel('PSD relative to peak [dB]');
title(sprintf(['Waveform Spectral Pre-characterization\n' ...
               'Nominal Channel Bandwidth: %.3f MHz\n' ...
               'Roll-off: %.0f%%, Filter Span: %.0f symbols'], ...
               channelBandwidth/1e6, ...
               rolloff*100, ...
               filterSpan));
xlim([-25 25]); % view limits
ylim([-100 5]); % view limits
hold off; % release the hold on



%% Time-domain whole waveform envelope
timeAxis_us = ((0:length(txWaveformIdeal)-1).' / Fs) * 1e6; % microseconds
figure;
% absolute value of a complex number is the envelope or magnitude
plot(timeAxis_us, abs(txWaveformIdeal), 'LineWidth', 1.1); 
grid on;
xlabel('Time [\mus]'); % since on the first line we multiplied by 1e6
ylabel('Magnitude'); % magnitude, abs or envelope
title('Whole Baseband Waveform Envelope');
xlim([timeAxis_us(1), timeAxis_us(end)]);

%% Select one burst
burstIndex = 1;
burstWaveform = idealShapedBursts{burstIndex};
burstTime_us = ((0:length(burstWaveform)-1).' / Fs) * 1e6;

%% Single-burst envelope
figure;
plot(burstTime_us, abs(burstWaveform), 'LineWidth', 1.1);
grid on;
xlabel('Time [\mus]');
ylabel('Magnitude');
title(sprintf('Baseband Waveform Envelope — Burst %d', burstIndex));
xlim([burstTime_us(1), burstTime_us(end)]);

%% Time-domain I/Q components of one burst
figure;
plot(burstTime_us, real(burstWaveform), 'LineWidth', 1.0);
hold on;
plot(burstTime_us, imag(burstWaveform), 'LineWidth', 1.0);
grid on;
xlabel('Time [\mus]');
ylabel('Amplitude');
title(sprintf('Baseband Waveform — Burst %d', burstIndex));
legend('In-phase I', 'Quadrature Q', 'Location', 'best');
xlim([burstTime_us(1), burstTime_us(end)]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  END SPECTRAL PRE-CHARACTERIZATION                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%