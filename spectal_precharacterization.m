
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
nominalChannelBandwidth_Hz = 13e6;
symbolRate = nominalChannelBandwidth_Hz / (1 + rolloff);
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
[~,centerIndex] = min(abs(frequencyAxis_Hz)); % ~ ignores the value, 
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

%% 6. Net throughput
% TODO! pay attention that these values remain calculated on the basis of
% the Symbol Rate set on the workspace by the LinkLevel.m simulator!
% fix these calculations accordingly!
% the same is for the PAPR! and maybe even other parameters!!! check
% everyone of them!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Net_Throughput_bps = packetLength / totalPacketTime; 
% Gross_Throughput_bps = numEncodedBits / totalPacketTime; 
% Gross_actual_ perBurst_Throughput_bps = bitsPerFramedBurst / burstTxTime;
% Net_perBurst_Throughput_bps =  (packetLength/numBursts) / burstTxTime;

%% Results
fprintf('\n--- SPECTRAL PRE-CHARACTERIZATION ---\n'); % title
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
fprintf('Net throughput: %.3f Mb/s\n', Net_Throughput_bps/1e6);

%% Single summary plot
figure; % open the graphical window
plot(frequencyAxis_Hz/1e6, ...
     psdSmoothed_dBc, ...
    'LineWidth',1.2);
grid on; % activate the grid on the graph
hold on; % allow to write over the already present plot without erasing it
% Nominal 13 MHz channel
%xline(-6.5,'--','-6.5 MHz');
%xline( 6.5,'--','+6.5 MHz');
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
ylabel('Normalized PSD [dBc/Hz]');
title('Waveform Spectral Pre-characterization nominalChannelBandwidth = 20.75 MHz');
xlim([-25 25]); % view limits
ylim([-100 5]); % view limits
hold off; % release the hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  END SPECTRAL PRE-CHARACTERIZATION                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%