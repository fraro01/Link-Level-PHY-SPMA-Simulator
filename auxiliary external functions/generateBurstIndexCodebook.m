% function for generating the codeword for the 5 indexing bits, in order to
% avoid a turbo coding scheme for just 5 indexing bits!
% this is based on the Hamming distance, because for few indexing bits does
% not ìmake any sense the computational overload of turbo coders...



function codebook = generateBurstIndexCodebook(numBursts, codewordLength)

    % Fixed seed: TX and RX must generate exactly the same codebook.
    previousRng = rng;
    rng(12345);

    % Number of random candidate codebooks to test.
    numTrials = 5000;

    bestCodebook = [];
    bestMinDistance = -inf;

    for trial = 1:numTrials

        % Generate distinct random binary codewords.
        candidateDecimals = randperm(2^codewordLength, numBursts) - 1;

        candidateCodebook = de2bi( ...
            candidateDecimals, ...
            codewordLength, ...
            'left-msb');

        % Compute minimum pairwise Hamming distance.
        minDistance = inf;

        for i = 1:numBursts-1
            for j = i+1:numBursts

                distance = sum( ...
                    candidateCodebook(i,:) ~= ...
                    candidateCodebook(j,:));

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

    % Restore the previous random state.
    rng(previousRng);

end