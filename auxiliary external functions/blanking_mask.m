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

% ASSUMPTION! we blank only the actual time in which the RX is actively
% transmitting. Random receiver erasures occupying X% of observation time

function blankingMask = blanking_mask(rxSignal, ...
                                      blankingPercentage, ...
                                      dutyCycle, ...
                                      burstDurationSamples)

    % blankingPercentage is expressed from 0 to 1.
    % Example: blankingPercentage = 0.25 means 25% blanking.

    % Check that the percentage is valid
    if blankingPercentage < 0 || blankingPercentage > 1
        error('blankingPercentage must be between 0 and 1.');
    end

    % Length of the received signal
    numSamples = length(rxSignal);
    % Number of samples that must be blanked
    numBlankedSamples = round(numSamples * blankingPercentage);    
    
    % 1 --> receiver available
    % 0 --> local transceiver is transmitting, hence: receiver blanked
    % Initialize the mask with all ones
    blankingMask = ones(size(rxSignal)); % the length must be the same of TX!

    % If there are samples to blank
    if numBlankedSamples > 0    
        if dutyCycle == 1    
            % With dutyCycle = 1, keep the original behavior:
            % all blanked samples belong to one consecutive block.    
            % Last possible starting index that allows the entire zero block
            % to remain inside the mask
            maxStartIndex = numSamples - numBlankedSamples + 1;    
            % Randomly choose the starting position
            startIndex = randi(maxStartIndex);    
            % Compute the ending position
            endIndex = startIndex + numBlankedSamples - 1;    
            % Insert the consecutive block of zeros
            blankingMask(startIndex:endIndex) = 0;    
        else
            % Minimum blanking duration, 0.5 --> half of a burst
            minBlankingBlockSamples = round(1.5 * burstDurationSamples);            
            % Number of blocks
            numBlocks = floor(numBlankedSamples / minBlankingBlockSamples);
            numBlocks = max(1, numBlocks);    % at least one block        
            % Length of each block
            blockLength = floor(numBlankedSamples / numBlocks);            
            % Divide the entire signal into numBlocks regions
            regionLength = floor(numSamples / numBlocks);            
            for b = 1:numBlocks            
                % Limits of the region assigned to this block
                regionStart = (b-1) * regionLength + 1;            
                if b < numBlocks
                    regionEnd = b * regionLength;
                else
                    regionEnd = numSamples;
                end            
                % Random start inside the assigned region
                maxStartIndex = regionEnd - blockLength + 1;            
                startIndex = randi([regionStart, maxStartIndex]);
                endIndex = startIndex + blockLength - 1;            
                % Apply blanking
                blankingMask(startIndex:endIndex) = 0;            
            end        
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            END RX LOCAL TRANSMIT ACTIVITY / BLANKING MASK               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%