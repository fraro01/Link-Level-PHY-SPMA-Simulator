% This is a secondary auxiliary file for the simulator, hence there must 
% not be any code before the 'function' declaration.

%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
* 
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            RX LOCAL TRANSMIT ACTIVITY / BLANKING MASK                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function blankingMask = blanking_mask(rxSignal, burstTxTime, ...
                        numBursts, slotDuration, Fs)
% 1 --> receiver available
% 0 --> local transceiver is transmitting, hence: receiver blanked
blankingMask = ones(size(rxSignal)); % the length must be the same of TX!
% Duration during which RX is considered blanked.
% We include the complete shaped transmitted burst.
rxTxBlankingDuration = burstTxTime; % inlcude also the tails of RRC
% Independent TX timing of the local transceiver
rxTxStartTimes = zeros(numBursts,1);
rxTxEndTimes   = zeros(numBursts,1);
for k = 1:numBursts
    slotStartTime = (k-1)*slotDuration;
    % Maximum start offset that keeps the TX pulse inside the slot
    maxRxTxOffset = slotDuration - rxTxBlankingDuration;
    rxTxOffset = rand() * maxRxTxOffset;
    % in the time domain
    rxTxStartTimes(k) = slotStartTime + rxTxOffset;
    rxTxEndTimes(k)   = rxTxStartTimes(k) + rxTxBlankingDuration;
    % sample conversion
    startSampleBlanking = floor(rxTxStartTimes(k)*Fs) + 1;
    endSampleBlanking   = ceil(rxTxEndTimes(k)*Fs);
    % constraints to keep them inside the boundaries
    startSampleBlanking = max(startSampleBlanking,1);
    endSampleBlanking   = min(endSampleBlanking,length(blankingMask));
    % the TX time is set to 0
    blankingMask(startSampleBlanking:endSampleBlanking) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            END RX LOCAL TRANSMIT ACTIVITY / BLANKING MASK               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%