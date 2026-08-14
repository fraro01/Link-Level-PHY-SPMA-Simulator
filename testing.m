clc;
clear;
close all;

%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
* 
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           TRANSMITTER                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Originate the Packet
% Set the length of the packet in bits
packetLength = 1331;
% Generate random bits (0 or 1)
packet = randi([0 1], packetLength, 1); %column vector not row for applying 
% turbo coding


%% Turbo Coding
% TODO! add trellis and numeiterations
% NOTE! study the turbo coding functionings
% Interleaver turbo
interleaverIndices = randperm(packetLength);
% Turbo encoder rate 1/3, (it is the default)
% Use as "InterleaverIndices" the one defined above
turboEnc = comm.TurboEncoder('InterleaverIndices', interleaverIndices);
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
                        length(phaseReference) + numIndexBits, numBursts);
for k = 1:numBursts
    % decimal to binary, with numIndexBits as #bits of result, most
    % significant bit is put left (used in TLC), then we transpose it: .'
    burstIndexBits = de2bi(k-1, numIndexBits, 'left-msb').';
    framedBursts(:,k) = [ % take all rows of cloumn k
    % vertical concatenation of three vectors
        syncReference;        % 48 bit
        phaseReference;       % 48 bit
        burstIndexBits;       % 5 bit: burst identifier
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          END TRANSMITTER                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           AVIONIC CHANNEL                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Channel pipeline:
% txSignal
%    -> TX power normalization
%    -> FSPL and Link Budget
%    -> Doppler
%    -> Rician multipath fading
%    -> oscillator CFO and Phase Offset
%    -> thermal AWGN
%    -> RX front-end gain

% The complex-baseband samples are interpreted as voltage amplitudes
% across a 1-ohm reference load.
% For a generic reference load R:
%     instantaneousPower_W = abs(x[n])^2 / R
% Since R = 1 ohm:
%     instantaneousPower_W = abs(x[n])^2
% This convention must be used consistently for both signal scaling (3.)
% and thermal-noise generation (9.)

%% 1. Avionic scenario parameters
% ASSUMPTION!
% Distance and velocities remain approximately constant during one packet.
% we could implement a sort of time variance for future implementations.
linkDistance_m = 50e3;          % TX-RX distance
aircraftSpeed_mps = 686;        % aircraft relative speed magnitude [m/s]
% Convention used here:
% radialVelocity_mps > 0  --> approaching
% radialVelocity_mps < 0  --> moving away

%% 2. Link-budget parameters
% Transmitter power
txPower_dBm = 48; % [dBm]                formula: x dBm = 10^(x/10) mW; 
% dBm (Decibels RELATIVE to 1 milliwatt)
% Antenna gains
txAntennaGain_dBi = 6; % [dBi] (Decibels RELATIVE to Isotropic antenna)
rxAntennaGain_dBi = 6; % RX antenna gain [dBi]
% Feeder, cable and connector losses
txFeederLoss_dB = 0.5;  % TX cable/connector loss [dB]
rxFeederLoss_dB = 0.5;  % RX cable/connector loss [dB]
% Additional propagation and implementation losses
otherPropagationLosses_dB = 1; % polarization, implementation, etc. [dB]
% Receiver noise and front-end parameters:
receiverNoiseFigure_dB = 2; % added by the receiver device, degrading SNR.
% Analog gain before the normalized ADC model. LNA/IF gain before ADC [dB]
rxFrontEndGain_dB = 90; % The received signal at the antenna is 
% approximately between -100 dBm and -80 dBm. Apply 80 dB of analog 
% front-end gain to bring it approximately between -20 dBm and 0 dBm before 
% the ADC. to make it compliant to our RX ADC.
% check it with: 10*log10(mean(abs(rxSignalAtAntenna).^2)) for properly set

%% 3. Normalize the transmitted waveform to the selected TX power
% txSignal currently has arbitrary normalized MATLAB units.
% Select only active transmitted samples (ommit empty tx intervals).
txActiveSamples = abs(txSignal) > 1e-12;
if ~any(txActiveSamples)
    error('The transmitted waveform contains no active samples.');
end
% Average power of the current normalized MATLAB waveform.
txSignalPowerArbitrary = mean(abs(txSignal(txActiveSamples)).^2);
% Convert transmitter power from dBm to watts. 
txPower_W = 10^((txPower_dBm - 30)/10); % standard formula
% FUNDAMENTAL SCALING: |txSignalPhysical|^2 is interpreted as instantaneous
% power in watts, and its active-sample average is txPower_W.
txSignalPhysical = txSignal .* sqrt(txPower_W / txSignalPowerArbitrary);
% square root is necessary since we are modifying the amplitude

%% 4. Free Space Path Loss and link budget
speedOfLight_mps = physconst('LightSpeed'); % for the wavelength calculus
% Sample-wise propagation attenuation.
largeScaleAmplitude = zeros(size(txSignalPhysical));
% one value for each burst.
fsplPerBurst_dB = zeros(numBursts,1);
% different bursts are on different carriers leading to different FSPL
for k = 1:numBursts %iterate over all the N bursts
    % Actual RF carrier used by burst k.
    carrierHz = carrierForBurst(k);
    % Wavelength corresponding to the current carrier.
    wavelength_m = speedOfLight_mps / carrierHz;
    % Free Space Path Loss in dB.
    fsplPerBurst_dB(k) = fspl(linkDistance_m, wavelength_m);
    % Propagation gain excluding txPower_dBm.
    % txPower_dBm is excluded because txSignalPhysical has already
    % been normalized to the actual transmitter power.
    % link budget calculation
    propagationGain_dB = ...
            txAntennaGain_dBi ...
            + rxAntennaGain_dBi ...
            - txFeederLoss_dB ...
            - rxFeederLoss_dB ...
            - otherPropagationLosses_dB ...
            - fsplPerBurst_dB(k); 
    % Convert power gain in dB into an amplitude multiplier.
    propagationAmplitudeGain = 10^(propagationGain_dB/20);
    % First sample occupied by the pulse-shaped burst.
    burstStartSample = round(startTimes(k)*Fs) + 1;
    % Include useful waveform and RRC filter tail.
    burstLengthSamples = ...
        (symbolsPerBurstWithCP + extraSymbols) * samplesPerSymbol;
    burstEndSample = min( ...
        burstStartSample + burstLengthSamples - 1, ...
        length(txSignalPhysical));
    % Assign the attenuation to samples belonging to this burst.
    largeScaleAmplitude(burstStartSample:burstEndSample) = ...
        propagationAmplitudeGain;
end
% Apply large-scale gains and losses. (expressed as amplitude, not power)
txSignalAfterLinkBudget = txSignalPhysical .* largeScaleAmplitude;

%% 5. Multipath parameters
% for a Rice model:
%   - If K = 0, the dominant path disappears, and the Rician model reduces 
%     perfectly to a Rayleigh model.
%   - As K → ∞, the signal resembles an uncorrupted, stable signal with 
%     almost zero fading.
% PathDelays values are expressed in seconds, using symbolRate.
% Here:
%   path 1 -> LOS, zero delay ASSUMPTION!
% The absolute propagation delay is neglected in the proposed channel 
% model. As the simulator targets link-level physical-layer performance, 
% the receiver is assumed to be synchronized to the first arriving path. 
% Consequently, only the excess delays between the multipath components are 
% modeled, since they determine the channel impulse response and directly 
% affect equalization and detection performance (relative delays).
%   path 2 -> reflection delayed by N symbol;
%   path 3 -> reflection delayed by N+X symbols.
% delays are taken relative from the symbol rate and samples per symbol and
% more importantly: are strictly related to the "Lest" parameter.
% TODO! in case of tweakings, set "Lest" parameter accordingly
pathDelays_s = [0, ...
                1/symbolRate, ...
                2/symbolRate];
% Average relative power of the paths relative to the strongest path, 
% LOS/echo paths [dB]
averagePathGains_dB = [0, ... % 0 gain with respect to reference (itself is taken as reference indeed it makes sense zero)
                       -12, ... % with -12 dB --> the path has 6.31% of the power of the first path, or equivalently a 16th
                       -18]; % with -18 dB --> the path has 1.58% of the power of the first path
% First path is Rician with K>0 (component specular/LOS is dominant), other
% paths are Rayleigh (because KFactor is supplied as a vector), with a 
% diffuse component without LOS. 
% here we've set "KFactor" as a linear ratio, not in dB, e.g.: linear_K=X
% --> the deterministic component is X times stronger (power) than the 
% diffuse casual component.
% the K-factor Rician determines the ratio between: power of the
% deterministic component (LOS) and the diffuse casual component. Since a
% Rician channel contains a deterministic component and a diffuse casual
% component
ricianKFactor = [30, ...
                 0, ... % the deterministic compoment LOS power is 0 --> Rayleigh channel
                 0]; % the deterministic compoment LOS power is 0 --> Rayleigh channel
% because a Rician channel with K --> 0 approaches a Rayleigh channel
% NOTE! actually the others two minor taps are reflected by the sea and
% probably by the wing of the aircraft --> we could model them with a
% predominant deterministic component with respect to the diffuse component
% as well! --> Rice instead of Rayleigh (K>0).
% Notice that the doppler shift is usually applied only for the Rice,
% because the Rayleigh is constituted only by many diffuse components by
% which attributing a single neat frequency shift is too difficult, hence
% in this case the doppler is only accounted as doppler spread spectrum.

%% 6. Doppler parameters
% Corresponding wavelength.
minimumWavelength_m = speedOfLight_mps / fMax;
centralWavelength_m = speedOfLight_mps / fCenter;
% Conservative maximum Doppler magnitude used to parametrize the
% diffuse fading process of all channel paths, for doppler spread spectrum
maximumDopplerShift_Hz = aircraftSpeed_mps / minimumWavelength_m;
% Signed Doppler shift of the direct LOS component
directPathDopplerShift_Hz = aircraftSpeed_mps / centralWavelength_m;

%% 7. Rician multipath fading channel
ricianChannel = comm.RicianChannel( ...
    'SampleRate', Fs, ... % how many samples per second the channel contains
    'PathDelays', pathDelays_s, ... % seconds
    'AveragePathGains', averagePathGains_dB, ... % average power of each path
    'NormalizePathGains', true, ... % normalize the total average power of all the paths of the multipath to 1 (otherwise is greater), see below
    'KFactor', ricianKFactor, ... % ratio between power of LOS component / power casual of diffuse component
    'DirectPathDopplerShift', [directPathDopplerShift_Hz*cosd(10) directPathDopplerShift_Hz*cosd(45) directPathDopplerShift_Hz*cosd(5)], ... % doppler of the LOS component of each Rician path; it is applied only for Rician paths
    'DirectPathInitialPhase', [0 pi/3 -pi/2], ... % initial phase of LOS component of each Rician path; they depend on the distance TX-RX. It represents the phase with which the specific path arrives at the RX when the simulation starts.
    'MaximumDopplerShift', maximumDopplerShift_Hz, ... % maximum Doppler magnitude of the diffuse fading components paths
    'PathGainsOutputPort', true, ... % returns also the instaneous gains of all the paths
    'RandomStream', 'Global stream'); % keeps generating casual fading components
% Apply time-varying multipath channel.
[rxSignalFaded, pathGains] = ricianChannel(txSignalAfterLinkBudget);
% rxSignalFaded includes:
% - LOS component;
% - echoes;
% - delays of multipath;
% - temporal variation of fading;
% - Doppler LOS;
% - Doppler spectrum of the diffused components.
% "pathGains" includes the complex instaneous coefficients of the paths

% regarding: "'NormalizePathGains', true":
% normalize the total average power of all the paths of the multipath to 1 
% (otherwise is greater), hence 0 dB; [1+0,0631+0,01585]≈1,07895 --> 
% P_1≈0,927;P_2≈0,0585;P_3≈0,0147; it avoids that the multipath introduces 
% more power with respect to the one already determined by link-budget, 
% since the channel can temporarily increase or decrease power due to 
% interferences with paths. We do this because the link-budget already 
% contains all attenuations, the multipath must not take into account 
% others arbitrary powers. Hence distribute the total power throught the 
% paths, but keep the average power equal to the one already determined 
% by the link-budget. But notice that the istantaneous value is not fixed!

%% 8. Oscillator CFO and phase offset
% Hardware-related oscillator Carrier Frequency Offset (CFO) and Phase 
% Impairments.
% TODO! set the following Offset parameters according to the simulation
% These are distinct from physical propagation Doppler.
trueCarrierPhaseOffset_rad = deg2rad(15); % the Phase offset introduced
trueCarrierFrequencyOffset_Hz = 1e3; % the CFO introduced in Hz
sampleIndex = (0:length(rxSignalFaded)-1).'; % we start from 0 --> -1
timeAxis = sampleIndex / Fs; % in seconds
% formula used below: ϕ(t) = ϕ0 + 2πΔft | CFO makes rotating continuously 
% the phase, Δf rotations per second; a complete rotation values 2π radians
% hence: dϕ/dt ​= 2π⋅1000
carrierImpairmentPhase = ...
    trueCarrierPhaseOffset_rad ...
    + 2*pi * trueCarrierFrequencyOffset_Hz * timeAxis;
% Impair the signal with the parameters defined above
rxSignalImpaired = rxSignalFaded .* exp(1j*carrierImpairmentPhase);

%% 9. Colored Noise from Interferers (external auxiliary script)
% colored noise of the intererfers, comes from the auxiliar external script
numInterferers = 5; % set the number of interferer nodes, NOTE! it
% resembles an ALOHA network!
interferer_waveform = interferers(numInterferers, txSCFDEBursts); % calls
% the external script
maxLength = max(length(interferer_waveform), length(rxSignalImpaired));
rxSignalImpaired(end+1:maxLength, 1) = 0;
interferer_waveform(end+1:maxLength, 1) = 0; 
rxSignalWithInterference = rxSignalImpaired + interferer_waveform;

%% 10. Thermal AWGN and receiver noise figure
% the SNR as a consequence of: transmitted power, distance, antennas gains,
% losses, bandwidth, temperature, noise figure.
% Convert receiver noise figure from dB to linear noise factor.
% Receiver-input noise power, standard formula:
% Pn = k * T * B * F    [W]
% where:
% - kT : noise spectral density;
% - kTB : total noise power ideal;
% - F is the linear noise factor, it adds to the ideal total noise power  
% the noise added in reception by RX device.
% notice that B is the width of the frequencies window through which the RX
% let pass the noise, (the greater --> the more noise we let get into RX).
% here it is calculated by comm.ThermalNoise() taking advantage of its
% relation with the sampling frequency Fs.

% Physical convention used in this simulation:
%   ReferenceLoad = 1 ohm
% Therefore, for every complex-baseband sample x[n]:
%   instantaneousPower_W = abs(x[n])^2 
% since P = V^2/R
% This is coherent with the previous TX normalization, where the active
% waveform was scaled so that its mean abs(x)^2 equals txPower_W.
% NoiseMethod = 'Noise figure':
% the receiver noise figure characterizes the noise introduced by the
% receiver electronics, assuming the standard reference temperature 290 K.
% Add290KAntennaNoise = true, adds also: thermal noise received from an
% equivalent 290 K antenna, hence it takes into account also the EM noise 
% of the environment; other than the equivalent receiver circuit noise 
% associated with the noise figure.
% Consequently, the total input-referred noise power is equivalent to:
%   Pnoise = k * 290 K * Fs * F 
% for ReferenceLoad = 1 ohm.
referenceLoad_Ohm = 1;
thermalNoiseChannel = comm.ThermalNoise( ...
    'NoiseMethod', 'Noise figure', ...
    'NoiseFigure', receiverNoiseFigure_dB, ... % equivalent noise of the RX electronics
    'SampleRate', Fs, ...
    'ReferenceLoad', referenceLoad_Ohm, ... % how to interpret the numerical amplitude of the signal
    'Add290KAntennaNoise', true); % whether we want to introdue the thermal 
% noise of the antenna at 290 K other than the one of the RX device.
% Add antenna and receiver thermal noise at the receiver input.
rxSignalAtAntenna = thermalNoiseChannel(rxSignalWithInterference);

%% 11. Receiver Analog Front-End FIXED Gain (Opposed to AGC, but 
% Complementary)
% Convert gain in dB into an amplitude multiplier.
% so this block performs: LNA, mixer, IF amplifier, baseband amplifier
rxFrontEndAmplitudeGain = 10^(rxFrontEndGain_dB/20);
% Apply the same gain to signal AND noise. So to bring the signal in a
% level compatible for my ADC. It is a fixed front-end analog gain.
rxSignal = rxSignalAtAntenna * rxFrontEndAmplitudeGain;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        END AVIONIC CHANNEL                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             RECEIVER                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Conversion in reception of EQUIVALENT-BASEBAND signal with RRC filter 
%% Realistic Burst Time Synchronization Using Preamble Correlation
%% Coarse AGC --> ADC --> Fine AGC
%% Carrier Synchronization (Frequency Synch. and Phase Synch.)
% The Rx knows:
%   - "Slot" duration
%   - Frequency hopping sequence
%   - Sync symbols
% The Rx DOES NOT know:
%   - startTimes(k) hence the start time of transmission of each burst
%   - intraSlotOffset(k) the offset of the bursts inside their own slot
%   - the exact sample of burst start

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
% where to store the syncohronized bursts
rxSCFDEBursts = zeros(symbolsPerBurstWithCP, numBursts);
% Set the threshold in between 0 e 1. The higher the amount of sync symbols 
% the higher could be set, and viceversa. e.g. with only 5 symbols --> 0.55
% for diagnostic purposes
syncDetectionThreshold = 0.3; % TODO! make it variable on the basis of the 
% Blanking Compensation!

% COARSE and FINE AGC configuration parameters (2.1. and 8.1.)
% Approximate desired power at the output of the Coarse AGC.
% A power equal to 0.25 correspondes to an RMS amplitude of sqrt(0.25)=0.5
% We do not use power=1, so to leave a certain margin to an eventual
% saturation of RX.
% coaorseGain adds a variable analog gain added to the fixed 80 dB frontend
coarseTargetPower = 0.25; % adimensional as a MATLAB unit (nomore Watt/dB)
% NOTE! it must be compliant to the variable: "adcFullScale" of the ADC,
% for the clipping, in order to avoid clipping complex samples.
% Coarse gain limits.
% the Coarse AGC can attenuate up to a factor of: 0.05
% or amplify up to a factor of: 30.
coarseMinGain = 0.05;
coarseMaxGain = 30;
% Fine gain limits.
% fineGain is a digital gain added after ADC conversion
% The Fine AGC has to perform only a small final correction, because the
% main part has been already corrected by the Coarse AGC.
fineMinGain = 0.5;
fineMaxGain = 2;

% Carrier Synchronization (Frequency and Phase  Synchronization) Parameters
% Configiguration (used also for Fine AGC)
% We use both:
%   - synchronization-reference symbols;
%   - phase-reference symbols.
% Since they are both known by the receiver.
% The symbols known by the receiver are the ones generated by:
% - 48 synchronization reference bits;
% - 48 phase reference bits.
% Hence in total we have:
% 48 + 48 = 96 known bits.
% With 8-PSK each symbol contains 3 bits:
% 96 / 3 = 32 known symbols.
numCarrierPilotSymbols = floor( ...
        (length(syncReference) + length(phaseReference)) ...
        / bitsPerSymbol);   
% Known transmitted symbols.
% They are the first symbols after the CP and are identical in every burst.
knownCarrierPilotSymbols = modulatedBursts(1:numCarrierPilotSymbols, 1);
% Positions of the known carrier-pilot symbols inside rxBurstSymbols.
% rxBurstSymbols has the structure:
% CP | sync | phase reference | index | payload
carrierPilotIndicesInBurst = cpLength + (1:numCarrierPilotSymbols);
% Storage for the Carrier Frequency Synchronization Estimates
estimatedCarrierPhaseOffset_rad = nan(numBursts, 1);
estimatedCarrierFrequencyOffset_Hz = nan(numBursts, 1);

% For burst erasure
burstDetected = false(numBursts,1);

% NOTE! this is the main external cycle that iterates over all the bursts, 
% and it's divided into 8 distinct sections, with some sub-sections as well
% it is very big!
for k = 1:numBursts
    reset(rxRRCsync); % clean the filter

    % 1. EXTRACTION OF THE ENTIRE TIME SLOT
    % We do not use startTimes(k).
    % the RX looks at the whole slot k.
    slotStartTime = (k-1) * slotDuration;
    slotEndTime   = k * slotDuration;
    slotStartSample = floor(slotStartTime * Fs) + 1;
    slotEndSample   = min(ceil(slotEndTime * Fs), length(rxSignal));
    rxSlotUpconverted = rxSignal(slotStartSample : slotEndSample);

    % 2. DOWNCONVERSION
    % RX knows the FH sequence. But doesn't know the exact time of start 
    % burst
    globalSampleIndices = (slotStartSample-1 : slotEndSample-1).';
    tSlot = globalSampleIndices / Fs; % Fs = [samples/second], do the math
    fRel = relativeCarrierFrequencies(frequencyChannels(k)); % FH scheme
    % Downconversion bring the entire time slot back to baseband
    % raw because we still have to apply Coarse AGC
    rxSlotBasebandRaw = rxSlotUpconverted .* exp(-1j * 2*pi * fRel *tSlot);

       % 2.1. SIMPLE COARSE AGC
        % NOTE! it is not causal since the code reads the whole slot, hence it
        % takes advantage of samples that will actually arrive in the future,
        % but a real recevier cannot do this in reality in real time!
        % instead it would update the gain sample by sample.
        % AGC Coarse is required in order to bring the received signal in an
        % amplitude interval compatible with the front-end of the RX, before
        % synchronization. It must guarantee that the analog signal does not
        % overcome the full scale of the ADC.
        % Put the coarse AGC after downconversion since it's easier in the
        % equivalent-baseband, because AGC depends on the modulus of the
        % signal, the frequency traslation does not modify amplitude.
        % Calculate the istantaneous power of all the samples of the slot.
        slotSamplePower = abs(rxSlotBasebandRaw).^2;
        % The slot contains:
        % - a part with the burst;
        % - a part mainly with noise.
        % In order to avoid that the long noise part lowers too much the 
        % estimate, we only consider samples with relatively high power.
        % the 75th percentile divides the samples as follow:
        % - the 75% has power lower than the threshold;
        % - the 25% has power higher than the threshold.
        % NOTE! TODO! adjust the percentile based on the Duty Cycle set.
        powerThreshold = prctile(slotSamplePower,75);
        % Select the samples probably belonging to the burst.
        activeSamples = slotSamplePower >= powerThreshold; % array of: 1/0
        % Safety check.
        % Normally activeSamples contains at least 25% of samples.
        % If for some unforseen reasons no sample is selected, we use the whole 
        % slot to calculate the average.
        if any(activeSamples) % if at least there is a: 1 in the array --> True
            measuredCoarsePower = mean( ...
                slotSamplePower(activeSamples));
        else
            measuredCoarsePower = mean(slotSamplePower);
        end
        % Calculate the required gain to obtain the objective power.
        % Remember: power is the square of the amplitude --> sqrt()
        % If:
        % P_out = gain^2 * P_in
        % Then:
        % gain = sqrt(P_out / P_in)
        coarseGain = sqrt( ...
            coarseTargetPower / ...
            max(measuredCoarsePower,1e-12));
        % Limit the gain to the allowed interval.
        coarseGain = min( ...
            max(coarseGain,coarseMinGain), ...
            coarseMaxGain);
        % Apply the same gain to the whole slot.
        rxSlotBaseband = coarseGain * rxSlotBasebandRaw; % so now we have 
        % tried to set the average power of the "rxSlotBaseband" to the 
        % "coarseTargetPower" that is adimensional, is just a MATLAB scale,
        % not anymore in Watt or dB as we did in the channel modeling.
    
       % 2.2. SIMPLE ADC MODEL
        adcNumBits = 10; % maximum number of available bits to represent values
        adcFullScale = 1; % full scale, hence range: [-1;+1], for clipping
        % ADC maximum positive integer code.
        adcMaxCode = 2^(adcNumBits-1) - 1; % maximum representable number of 
        % ADC, in order to be symmetric to 0 --> -511,..,0,...,+511 (total 1023
        % not 1024).
        % Quantization step. distance between two consecutive quantized levels
        adcStep = adcFullScale / adcMaxCode;    
        % Separate I and Q components.
        rxADCReal = real(rxSlotBaseband);
        rxADCImag = imag(rxSlotBaseband);    
        % Input clipping. [-1;+1], very clever with only one row of code
        rxADCReal = min(max(rxADCReal,-adcFullScale),adcFullScale);
        rxADCImag = min(max(rxADCImag,-adcFullScale),adcFullScale);    
        % Uniform quantization. Uniform since the quantization step is equal.
        % how many quantization steps corresponds the input value, round them
        % to the nearest integer, then go back to the normalized value
        rxADCReal = adcStep * round(rxADCReal/adcStep);
        rxADCImag = adcStep * round(rxADCImag/adcStep);    
        % Complex ADC output. after clipping and quantization
        rxSlotBasebandADC = rxADCReal + 1j*rxADCImag;   

    % 3. RRC FILTER WITHOUT DECIMATION
    % The exit still keeps the #samplesPerSymbol samples per symbol
    rxSlotMatched = rxRRCsync(rxSlotBasebandADC);

    % 4. CORRELATION METRIC CALCULUS
    % For each possible sample q, we take:
    % q,
    % q + samplesPerSymbol,
    % q + 2*samplesPerSymbol,
    % ...
    % So to obtain a possible sequence of received symbols
    % We confront it with: knownSyncSymbols.
    numMatchedSamples = length(rxSlotMatched);
    % the first sync symbol must have as header before itslef the whole CP,
    % hence we go after the CP for position correlation.
    firstCandidate = cpLength * samplesPerSymbol + 1;
    % After the first sync symbol there must be available:
    %   - All sync symbols;
    %   - The whole rest part of the burst.
    % In other words we are extracting the entire burst from its "slot".
    lastCandidate = numMatchedSamples ...
        - (symbolsPerBurstWithCP - cpLength - 1) * samplesPerSymbol; % here 
    % we are taking the POSSIBLE LAST POSITION of the burst inside its slot 
    % in terms of samples, starting directly AFTER the CP.
    numCandidates = lastCandidate - firstCandidate + 1; % number of 
    % candidated POSITIONS for starting the preamble (in samples)
    syncMetric = zeros(numCandidates, 1);
    for candidateIdx = 1:numCandidates
        q = firstCandidate + candidateIdx - 1;
        % possible samples of sync symbols 
        syncSampleIndices = q + ...
            (0:numSyncSymbols-1).' * samplesPerSymbol;
        % Possible samples candidated
        receivedSyncCandidate = ...
            rxSlotMatched(syncSampleIndices);
        % Local energy of the received sequence
        receivedEnergy = sum(abs(receivedSyncCandidate).^2);
        % if the energy is above this threshold we proceed, otherwise no
        if receivedEnergy > 1e-12
            % normalize it, to compare with the normalized known sequence
            receivedSyncNormalized = ...
                receivedSyncCandidate / sqrt(receivedEnergy);
            % Complex correlation normalized.
            % the absolute value, renders the metric insensible to a phase
            % rotation in common.
            correlationValue = ...
                knownSyncSymbolsNormalized' * receivedSyncNormalized; % the 
            % knownSyncSymbolsNormalized' makes the conjugated transpose
            syncMetric(candidateIdx) = abs(correlationValue)^2; % in order 
            % to cope with phase rotation of the channel, the modulus
            % reimans unchanged while the angle could change
        end
    end

    % 5. LOOK FOR THE MAXIMUM CORRELATION
    % take the maximum and its position
    [bestMetric, bestCandidateIdx] = max(syncMetric);
    % position of the first sync symbol in the output of the adapted filter
    firstSyncSample = firstCandidate + bestCandidateIdx - 1;

    % 6. DETECTION DECISION
    if bestMetric < syncDetectionThreshold
        warning( ...
            ['Burst %d not detected: synchronization metric ' ...
             '= %.3f.'], ...
            k, ...
            bestMetric);
        % we mark the burst as erased
        burstDetected(k) = false;
        continue;
    end
    % we mark the burst as found
    burstDetected(k) = true;

    % 7. FIND WHERE THE CP START
    % The preamble is the first element after the CP
    firstCPSample = firstSyncSample - cpLength * samplesPerSymbol;

    % 8. EXTRACTION OF THE SYMBOLS OF THE BURST
    % now we know both:
    %   - temporal position;
    %   - correct sampling phase.
    % We sample a value each #samplePerSymbols samples
    burstSampleIndices = firstCPSample + ...
        (0:symbolsPerBurstWithCP-1).' * samplesPerSymbol;

    if burstSampleIndices(1) < 1 || ... % || (or)
        burstSampleIndices(end) > numMatchedSamples
        warning( ...
            'Burst %d detected, but not entirely contained in the slot.', ...
            k);
        burstDetected(k) = false; % burst marked as erased
        continue;
    end

    rxBurstSymbols = rxSlotMatched(burstSampleIndices);

       % 8.1. SIMPLE FINE AGC 
        % AGC Fine is required in order to finalize the coarse AGC, over the
        % precise found synchronization symbols.
        % Fine AGC must be implemented once we have understood which samples
        % belong to the burst, the timing must be known!
        % Extract the pilots received.
        receivedFinePilotSymbols = rxBurstSymbols(carrierPilotIndicesInBurst);
        % Average power of the received pilots.
        measuredFinePilotPower = mean(abs(receivedFinePilotSymbols).^2);
        % Average power of the transmitted pilots.
        % for a normal 8-PSK modulation, should be around 1.
        % In any case it is better to directly calculate it.
        targetFinePilotPower = mean(abs(knownCarrierPilotSymbols).^2);
        % Calculate the fine gain.
        fineGain = sqrt( ...
            targetFinePilotPower / ...
            max(measuredFinePilotPower,1e-12));
        % Limit the correction
        fineGain = min( ...
            max(fineGain,fineMinGain), ...
            fineMaxGain);
        % apply the same gain to the whole burst.
        rxBurstSymbols = fineGain * rxBurstSymbols;
    
           % 8.2. FEED-FORWARD CARRIER SYNCHRONIZATION (ACQUISITION)
            % 8.2.1. EXTRACT THE KNOWN CARRIER PILOTS
            receivedCarrierPilotSymbols = ... % take only the known pilots!
                rxBurstSymbols(carrierPilotIndicesInBurst);
            
            % 8.2.2. REMOVE THE KNOWN PILOT MODULATION
            % keep in mind the error we must correct: ϕ(t) = ϕ0 + 2πΔft
            % If:
            % receivedPilot(m) ≈ knownPilot(m) * exp(j*(2πΔft+ϕ0​))
            % then:
            % receivedPilot(m) * conj(knownPilot(m)) ≈ exp(j*(2πΔft+ϕ0​))
            % because 8-PSK symbols have unit magnitude!
            carrierPilotProducts = ...
                receivedCarrierPilotSymbols .* conj(knownCarrierPilotSymbols); % so
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
            
            % 8.2.3. MEASURE AND UNWRAP THE PILOT PHASE
            wrappedPilotPhase = angle(carrierPilotProducts);
            % angle() returns values between -pi and +pi.
            % unwrap() removes artificial jumps of approximately 2*pi, converting
            % the phase evolution into an approximately continuous line.
            unwrappedPilotPhase = unwrap(wrappedPilotPhase);
            
            % 8.2.4. ESTIMATE PHASE INTERCEPT AND PHASE SLOPE
            % The first known pilot has relative symbol index zero.
            carrierPilotRelativeIndex = (0:numCarrierPilotSymbols-1).';
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
        
            % 8.2.5. CONVERT PHASE SLOPE TO CFO
            % phase_variation_per_symbol = [2π * ​Δf] / Symbol_Rate
            % thus: Δf = [phase_variation_per_symbol * Symbol_Rate] / 2*pi​
            estimatedCarrierFrequencyOffset_Hz(k) = ...
                estimatedPhaseSlope_radPerSymbol * symbolRate / (2*pi);
            % Wrap the initial phase estimate to the interval [-pi, +pi].
            estimatedCarrierPhaseOffset_rad(k) = ...
                angle(exp(1j*estimatedPhaseIntercept_rad));
            
            % 8.2.6. CORRECT CP and useful burst symbols
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
                rxBurstSymbols .* ...
                exp(-1j*estimatedCarrierPhaseEvolution);

    % finally at the very end now save the corrected received burst
    rxSCFDEBursts(:,k) = rxBurstSymbolsCarrierCorrected;

end

%testing
%{%}
errMax = max(abs(rxSCFDEBursts(:) - txSCFDEBursts(:)));
errRMS = sqrt(mean(abs(rxSCFDEBursts(:) - txSCFDEBursts(:)).^2));
fprintf('Max error: %.3e\n', errMax);
fprintf('RMS error: %.3e\n', errRMS);



%% PILOT-BASED Channel Estimation using Sync + Phase reference bits and 
%% SC-FDE Demodulation with MMSE Equalization
N = symbolsPerBurst;
Ncp = cpLength;
% Known pilot bits available in each burst
knownPilotBits = [syncReference; phaseReference];
numPilotSymbols = floor(length(knownPilotBits)/bitsPerSymbol);
knownPilotSymbols = modulatedBursts(1:numPilotSymbols, 1);
numPilotSymbols = length(knownPilotSymbols);

% NOTE! the pilot tones (symbols) must be >> number of taps (Lest). 
% there are two possible ways to operate:
% 1. Use dedicated Pilot Bursts;
% 2. Augment the length of the sync and reference bits, at least covering a
% total of 16-32 symbols.
% by doing so you will notice a strong decrease in numErrors8PSK even with
% higher AWGN.
% Length ESTimate
Lest = 3; % length of the discrete impulsive response we want to estimate. 
% it represents #Lest taps at symbols level, hence we estimate the channel 
% memory length in symbols, hence it is totally different from what we did 
% in the Channel implementation before the RRC RX. Indeed we are doing this
% after the RRC RX

rxEqualizedBursts = zeros(symbolsPerBurst, numBursts); % equalized bursts
H_est_bursts = zeros(N, numBursts); % channel estimation in frequency
h_est_bursts = zeros(Lest, numBursts); % channel estimation in time
% estimated noise variance after MMSE equalization for 8-PSK LLR demod
noiseVarAfterEq = zeros(numBursts,1);
for k = 1:numBursts % cycle over all bursts
    % if the current burst was not found (erased), we directly skip its
    % channel estimation/MMSE
    if ~burstDetected(k)
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
    % indeed! for the MMSE denominator: (σ_n)^2/(σ_x)^2
    pilotResidual = y - A*h_hat;
    % Estimate effective noise + unmodelled interference variance
    % If: (σ_n)^2 --> 0, MMSE converges to ZF
    noiseVarHat = mean(abs(pilotResidual).^2);
    noiseVarHat = max(noiseVarHat, 1e-12); % floor for safety
    % Symbol power (again for MMSE)
    Px = mean(abs(knownPilotSymbols).^2);
    h_est_bursts(:,k) = h_hat; % saves the times response estimates by 
    % burst k
    % Convert time-domain channel estimate to frequency domain
    h_hat_padded = [h_hat; zeros(N-Lest,1)];
    H_hat = fft(h_hat_padded, N); % channel frequency response
    H_est_bursts(:,k) = H_hat; % saves the frequency response estimates by 
    % burst k
    % SC-FDE Equalization
    Y = fft(rxBurst, N);
    % MMSE equalizaer
    W = conj(H_hat) ./ (abs(H_hat).^2 + noiseVarHat/Px);
    Xhat = W .* Y; % equalization
    % Back to time/symbol domain
    rxBurstEqualized = ifft(Xhat, N);    

    % POST-MMSE residual gain estimation using pilots (for potential 
    % ampitude gain and phase rotations on the burst)
    % Extract equalized pilots
    rxPilotEqualized = rxBurstEqualized(1:numPilotSymbols);    
    % Estimate residual complex gain:
    % rxPilotEqualized ~= gHat * knownPilotSymbols + residual
    % Least Sqaures estimate of gHat, (notice we cannot simplify since are 
    % vectors and matrices operations)
    gHat = (knownPilotSymbols' * rxPilotEqualized) / ...
           (knownPilotSymbols' * knownPilotSymbols);    
    % Numerical safety if too low gain
    if abs(gHat) < 1e-8
        gHat = 1;
    end
    % Normalize the whole burst    
    rxBurstEqualized = rxBurstEqualized / gHat;
    % Store normalized equalized burst
    rxEqualizedBursts(:,k) = rxBurstEqualized;   
    % Effective post-MMSE variance estimation    
    % Equalized and normalized pilots
    rxPilotEqualized = rxBurstEqualized(1:numPilotSymbols);    
    % Effective residual actually seen by the soft demapper (8-PSK demod)
    % variance from ground-truth!
    postEqResidual = rxPilotEqualized - knownPilotSymbols;    
    % Effective post-equalizatoin error variance (is a power since it is 
    % the square modulus of a complex number representing symbols)
    noiseVarAfterEq(k) = ...
    sum(abs(postEqResidual).^2) / max(numPilotSymbols - 1, 1);  % -1 since 
    % we have already consumed one degree of freedom in Least Square estim
    % Numerical floor, to avoid a too low variance that could explose LLR
    noiseVarAfterEq(k) = max(noiseVarAfterEq(k), 1e-8);
end


%% 8-PSK Demodulation
rxModulatedBursts = rxEqualizedBursts; 
% due to channel effects
rxFramedBurstsPaddedLLR = zeros(symbolsPerBurst*bitsPerSymbol, numBursts);
for k = 1:numBursts
    % if the current burst was not found (erased), we directly skip its
    % demodulation
    if ~burstDetected(k)
        % Erased burst --> leave all LLRs equal to zero, initializated in:
        % 'rxFramedBurstsPaddedLLR'; sync failure⇒burst erasure⇒LLR=0⇒
        % ⇒Turbo decoder receives no information from that burst.
        continue;
    end
    rxSymbols_k = rxModulatedBursts(:,k);
    % estimated noise variance at the MMSE equalizer output
    noiseVariance_k = noiseVarAfterEq(k); 
    % MATLAB pskdemod LLR convention:
    %   positive LLR -> bit 0 more likely
    %   negative LLR -> bit 1 more likely
    burstLLR = pskdemod( ...
        rxSymbols_k, ...
        M, ...
        pi/M, ...
        'gray', ...
        'OutputType','approxllr', ... % instead of 'llr' is simpler and more robust, it uses max-log approximation instead of the sum of the lilelihoods
        'NoiseVariance', noiseVariance_k);
    rxFramedBurstsPaddedLLR(:,k) = burstLLR; % a soft value per each bit
end
% Remove 8-PSK padding
rxFramedBurstsLLR = rxFramedBurstsPaddedLLR(1:end-numPadBits8PSK, :);


%% Turbo decoder
% Recover encoded packet from bursts
headerLength = length(syncReference) + ...
               length(phaseReference) + ...
               numIndexBits;
%remove headers
rxBurstsPayloadLLR = rxFramedBurstsLLR(headerLength+1:end,:);
% stack into column
rxEncodedPacketPaddedLLR = reshape(rxBurstsPayloadLLR,[],1); 
% leave padding
rxEncodedPacketLLR = rxEncodedPacketPaddedLLR(1:numEncodedBits); 

% Diagnostic check-point: convert LLR sign back to hard decisions
%{%}
rxEncodedPacketHardFromLLR = rxEncodedPacketLLR < 0;
numErrorsLLRHard = sum(rxEncodedPacketHardFromLLR ~= encodedPacket);
fprintf('Hard errors inferred from LLR sign: %d\n', numErrorsLLRHard);

% Random burst loss simulation
numLostBursts = 0; % TODO! choose manually (remove it once added colored N)
lostBurstIdx = randperm(numBursts, numLostBursts);
payloadBitsPerBurst = size(rxBurstsPayloadLLR,1);
for b = lostBurstIdx % iterate over erased bursts
    firstBit = (b-1)*payloadBitsPerBurst + 1; % take also the first burst, 
    % hence 0, furthermore MATLAB starts from 1 not 0 
    lastBit  = min(b*payloadBitsPerBurst, numEncodedBits); % handle also 
    % the case of the last burst in which there was padding (minimum len)
    % LLR=0 means that neither bit 0 nor bit 1 is preferred
    rxEncodedPacketLLR(firstBit:lastBit) = 0; % set these bits to zero
end
% numIterations = 6; % the default is 6
turboDec = comm.TurboDecoder( ...
    'InterleaverIndices', interleaverIndices);
    %'NumIterations', numIterations);
% convention used by turboDecoder MATLAB is with the reversed values (-)
decodedPacket = turboDec(-rxEncodedPacketLLR); 

% Diagnostic chek-point: just to compare the whole packet from Tx to Rx
%{%}
numFinalErrors = sum(packet ~= decodedPacket)