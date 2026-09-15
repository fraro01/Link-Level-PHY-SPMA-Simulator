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

function blankingMask = blanking_mask(rxSignal, blankingPercentage)
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

        % Last possible starting index that allows the entire zero block
        % to remain inside the mask
        maxStartIndex = numSamples - numBlankedSamples + 1;

        % Randomly choose the starting position
        startIndex = randi(maxStartIndex); % random integer between 1 and 
        % 'maxStartIndex' 

        % Compute the ending position
        endIndex = startIndex + numBlankedSamples - 1;

        % Insert the consecutive block of zeros
        blankingMask(startIndex:endIndex) = 0;
    end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            END RX LOCAL TRANSMIT ACTIVITY / BLANKING MASK               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%