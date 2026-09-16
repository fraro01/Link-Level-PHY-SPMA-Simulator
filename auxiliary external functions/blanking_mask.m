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
% transmitting.

function blankingMask = blanking_mask(rxSignal, blankingPercentage, dutyCycle)
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
            % With dutyCycle < 1, distribute the blanking randomly
            % throughout the complete received signal.    
            % randperm returns numBlankedSamples different random indices
            % between 1 and numSamples.
            blankedIndices = randperm(numSamples, numBlankedSamples);    
            % Blank exactly the selected samples.
            blankingMask(blankedIndices) = 0;    
        end    
    end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            END RX LOCAL TRANSMIT ACTIVITY / BLANKING MASK               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%