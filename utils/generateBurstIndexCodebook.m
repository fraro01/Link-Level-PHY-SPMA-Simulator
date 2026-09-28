% This is a secondary auxiliary file for the simulator, hence there must 
% not be any code before the 'function' declaration.

%{
* ctrl + f --> TODO! for coding;
* ctrl + f --> NOTE! something over which we must pay attention;
* ctrl + f --> ASSUMPTION! when assumptions are introduced;
* 
%}



% function for generating the codeword for the indexing bits, in order to
% avoid a turbo coding scheme for just 5 indexing bits (2^5>#28)!
% this is based on the Hamming distance, because for few indexing bits does
% not make any sense the computational overload of turbo coders...

% with a Monte Carlo search of the codebook we maximize the minimum Hamming
% distance

function [codebook, ...
          friendIndexCodeword] = generateBurstIndexCodebook(numBursts, ...
                                                            codewordLength)

    % Fixed seed: TX and RX must generate exactly the same codebook.
    % Without altering the current state of the random generator, at the
    % end of this script we'll assign bac again the 'previousRng' state to
    % the simulation.
    previousRng = rng;
    rng(42);

    % Number of random candidate codebooks to test.
    numTrials = 1; % NOTE! with 'rng(42)' I have already tested that in the 
    % very first iteration I already get: 
    % "Burst-index codebook: minimum Hamming distance = 6" , hence it does
    % not make any sense to increase 'numTrials'

    bestCodebook = [];
    bestMinDistance = -inf;

    % Generate the complete binary space only once.
    allWords = de2bi(0:(2^codewordLength - 1), ...
                     codewordLength, ...
                     'left-msb');

    for trial = 1:numTrials
        % extract a random codebook to compare at first
        firstIndex = randi(size(allWords,1)); % it will remain as the first
        candidateCodebook = allWords(firstIndex,:); % extract it
        % where to extract for the comparisons
        availableWords = allWords;
        % set the comparison to the current extraction to null
        availableWords(firstIndex,:) = []; 
    
        % since at each iteration it enlarges --> keep iterating as long as
        % we are below #numBursts
        while size(candidateCodebook,1) < numBursts
            % array for storing the distances
            minDistances = zeros(size(availableWords,1),1);
            % min-max logic
            for c = 1:size(availableWords,1)
                distances = sum(candidateCodebook ~= availableWords(c,:), 2);
                minDistances(c) = min(distances);
            end
            % extract the max from the mins
            bestDistance = max(minDistances);
            % in case we have multiple distances we the same values, we
            % choose one randomly
            bestCandidates = find(minDistances == bestDistance);
            selectedIndex = bestCandidates(randi(length(bestCandidates)));
            % store the found codebook
            candidateCodebook(end+1,:) = availableWords(selectedIndex,:); %#ok<AGROW>
            availableWords(selectedIndex,:) = []; % remove it from comparison
        end

        % Generate distinct random binary codewords.
        % it generates an array of length 28, with its 28 cells populated
        % by random decimal numbers ranging between 0 and (2^15 - 1)
       % candidateDecimals = randperm(2^codewordLength, numBursts) - 1;
        
        % from the above array, develop the binary conversion column by
        % column, through its cells. Then trasponse in rows.
       % candidateCodebook = de2bi(candidateDecimals, ...
       %                           codewordLength, ...
       %                           'left-msb');

        % Compute minimum pairwise Hamming distance, initialize it to inf
        minDistance = inf;
        % if within the current loop we hvae already found a bad candidate
        discardCandidate = false; % Early Stopping
        
        % NOTE! this is part of the old structure that I used to generate
        % the indexing, with a brute force Monte Carlo basically...
        for i = 1:numBursts-1
            for j = i+1:numBursts % we skip all the candidates before --> avoid counting them twice
                distance = sum(candidateCodebook(i,:) ~= candidateCodebook(j,:));
                minDistance = min(minDistance, distance);
                % in case we have already found a bad candidate, mark it!
                % Early Stopping
                if minDistance <= bestMinDistance
                    discardCandidate = true;
                    break;
                end
            end
            % exit the loop in case of alreayd found bad candidate
            % Early Stopping
            if discardCandidate
                break;
            end
        end

        % Keep the best codebook found.
        if minDistance > bestMinDistance
            bestMinDistance = minDistance;
            bestCodebook = candidateCodebook;
        end
    end
    
    % fprintf('Burst-index codebook: minimum Hamming distance = %d\n', bestMinDistance);
    codebook = bestCodebook;

    %% Reserved FRIEND / INTERFERER codeword
    % Generate all possible binary words having the same length
    % as the protected burst-index codewords.
    allCandidates = allWords;
    % For every possible candidate, calculate its minimum Hamming
    % distance from the valid desired-transmitter codewords.
    minDistanceFromDesiredCodebook = zeros(size(allCandidates,1),1);    
    for c = 1:size(allCandidates,1)    
        % Hamming distance between candidate c and EVERY valid
        % desired burst-index codeword.
        distances = sum(codebook ~= allCandidates(c,:), 2); % 2 --> sum over the columns for each row!
        % We care about the closest desired codeword.
        minDistanceFromDesiredCodebook(c) = min(distances);
    end    
    % Select the candidate having the largest possible minimum
    % Hamming distance from the desired burst-index codebook.
    % take the maximum of the minimum distances!
    [~, bestFriendCandidate] = max(minDistanceFromDesiredCodebook);
    friendIndexCodeword = allCandidates(bestFriendCandidate,:);

    % Restore the previous random state.
    rng(previousRng);

end

