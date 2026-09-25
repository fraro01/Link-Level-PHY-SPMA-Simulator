% function for generating the codeword for the 5 indexing bits, in order to
% avoid a turbo coding scheme for just 5 indexing bits!
% this is based on the Hamming distance, because for few indexing bits does
% not make any sense the computational overload of turbo coders...

% with a Monte Carlo search of the codebook we maximize the minimum Hamming
% distance


function [codebook, friendIndexCodeword] = generateBurstIndexCodebook(numBursts, codewordLength)

    % Fixed seed: TX and RX must generate exactly the same codebook.
    previousRng = rng;
    rng(42);

    % Number of random candidate codebooks to test.
    numTrials = 5000;

    bestCodebook = [];
    bestMinDistance = -inf;

    for trial = 1:numTrials

        % Generate distinct random binary codewords.
        % it generates an array of length 15, with its 15 cells populated
        % by random decimal numbers ranging between 0 and 2^15
        candidateDecimals = randperm(2^codewordLength, numBursts) - 1;

        candidateCodebook = de2bi(candidateDecimals, ...
                                  codewordLength, ...
                                 'left-msb');

        % Compute minimum pairwise Hamming distance.
        minDistance = inf;

        for i = 1:numBursts-1
            for j = i+1:numBursts
                distance = sum( ...
                    candidateCodebook(i,:) ~= candidateCodebook(j,:));
                minDistance = min(minDistance, distance);
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
    allCandidates = de2bi(0:(2^codewordLength - 1), ...
                          codewordLength, ...
                          'left-msb');
    % For every possible candidate, calculate its minimum Hamming
    % distance from the valid desired-transmitter codewords.
    minDistanceFromDesiredCodebook = zeros(size(allCandidates,1),1);    
    for c = 1:size(allCandidates,1)    
        % Hamming distance between candidate c and every valid
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