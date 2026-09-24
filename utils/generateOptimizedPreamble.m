function [syncReference, phaseReference, info] = ...
    generateOptimizedPreamble(M, LestDesign, numTrials)

% GENERATEOPTIMIZEDPREAMBLE
%
% Generates the 96-bit preamble used by the Link-Level simulator:
%
%       48 syncReference bits
%       48 phaseReference bits
%
% The two 48-bit blocks are optimized JOINTLY as a single 96-bit
% preamble.
%
% In the current simulator:
%
%       96 bits --> 32 Gray-coded 8-PSK symbols
%
% The optimization is performed directly on the MODULATED symbols,
% because these are the quantities actually used by:
%
%   1. burst synchronization;
%   2. carrier synchronization;
%   3. LS channel estimation.
%
%
% PRIMARY OPTIMIZATION CRITERION
% --------------------------------
%
% Minimize:
%
%       J_LS = trace(inv(A' * A))
%
% where A is constructed exactly as in the LS channel estimator of
% receiver.m.
%
% For white noise:
%
%       Cov(h_hat - h) = sigma_n^2 * inv(A' * A)
%
% therefore J_LS is proportional to the total variance of the
% estimated channel taps.
%
%
% SECONDARY OPTIMIZATION CRITERIA
% --------------------------------
%
% Among candidates having essentially the same J_LS, prefer:
%
%   1. lower condition number of A'*A;
%   2. lower maximum normalized autocorrelation sidelobe.
%
% The autocorrelation criterion is useful for burst synchronization.
%
%
% INPUTS
% --------------------------------
%
% M
%       Modulation order.
%       Current simulator: M = 8.
%
% LestDesign
%       Maximum symbol-spaced channel length for which the preamble
%       is designed.
%
%       IMPORTANT:
%       this is a DESIGN parameter and must NOT depend on the
%       instantaneous realization of the channel.
%
% numTrials
%       Number of random 96-bit candidate preambles tested.
%       Suggested:
%
%           1e5       quick design
%           5e5       better search
%           1e6       offline final design
%
%
% OUTPUTS
% --------------------------------
%
% syncReference
%       First 48 optimized bits.
%
% phaseReference
%       Last 48 optimized bits.
%
% info
%       Structure containing optimization diagnostics.
%
%
% NOTE
% --------------------------------
%
% This function restores the MATLAB random-number-generator state
% before returning. Therefore calling this function does not alter
% the randomness of the rest of the Link-Level simulation.


    %% ===============================================================
    %  1. PARAMETERS
    % ================================================================

    numPreambleBits = 96;

    bitsPerSymbol = log2(M);


    %% ===============================================================
    %  2. INPUT CHECKS
    % ================================================================

    if mod(bitsPerSymbol, 1) ~= 0
        error('M must be a power of two.');
    end

    if mod(numPreambleBits, bitsPerSymbol) ~= 0
        error(['The 96-bit preamble length must be divisible by ' ...
               'log2(M).']);
    end

    numPilotSymbols = numPreambleBits / bitsPerSymbol;

    if LestDesign < 1 || LestDesign > numPilotSymbols
        error(['LestDesign must satisfy 1 <= LestDesign <= ' ...
               'numPilotSymbols.']);
    end

    if numTrials < 1 || mod(numTrials,1) ~= 0
        error('numTrials must be a positive integer.');
    end


    %% ===============================================================
    %  3. SAVE CURRENT RNG STATE
    %
    %  The optimization is deterministic and does not modify the
    %  random evolution of the main simulation.
    % ================================================================

    previousRngState = rng;

    rng(12345, 'twister');


    %% ===============================================================
    %  4. INITIALIZE SEARCH
    % ================================================================

    bestBits = [];
    bestSymbols = [];
    bestA = [];

    bestJls = inf;
    bestConditionNumber = inf;
    bestMaxSidelobe = inf;

    % Numerical tolerance used when comparing J_LS values.
    relativeTolerance = 1e-10;


    %% ===============================================================
    %  5. MONTE CARLO SEARCH
    % ================================================================

    for trial = 1:numTrials


        %% -----------------------------------------------------------
        %  5.1 Generate random 96-bit candidate
        % ------------------------------------------------------------

        candidateBits = randi([0 1], numPreambleBits, 1);


        %% -----------------------------------------------------------
        %  5.2 Map candidate EXACTLY as transmitter.m
        %
        %  96 bits
        %       -->
        %  groups of log2(M) bits
        %       -->
        %  Gray-coded M-PSK symbols with pi/M rotation
        % ------------------------------------------------------------

        bitGroups = reshape( ...
            candidateBits, ...
            bitsPerSymbol, ...
            []).';

        symbolIndices = bi2de( ...
            bitGroups, ...
            'left-msb');

        candidateSymbols = pskmod( ...
            symbolIndices, ...
            M, ...
            pi/M, ...
            'gray');


        %% -----------------------------------------------------------
        %  5.3 Build LS matrix A
        %
        %  EXACTLY the same structure used in receiver.m:
        %
        %  y(n) =
        %      h(1)x(n)
        %    + h(2)x(n-1)
        %    + ...
        %    + h(L)x(n-L+1)
        %
        % ------------------------------------------------------------

        numRows = numPilotSymbols - LestDesign + 1;

        A = zeros(numRows, LestDesign);

        row = 1;

        for n = LestDesign:numPilotSymbols

            for ell = 0:LestDesign-1

                A(row, ell+1) = ...
                    candidateSymbols(n-ell);

            end

            row = row + 1;

        end


        %% -----------------------------------------------------------
        %  5.4 Gram matrix
        %
        %               G = A^H A
        %
        %  Ideally:
        %
        %               G ~= c I
        %
        % ------------------------------------------------------------

        G = A' * A;


        %% -----------------------------------------------------------
        %  5.5 Reject badly conditioned / singular candidates
        % ------------------------------------------------------------

        reciprocalCondition = rcond(G);

        if ~isfinite(reciprocalCondition) || ...
                reciprocalCondition < 1e-12

            continue;

        end


        %% -----------------------------------------------------------
        %  5.6 LS QUALITY METRIC
        %
        %  Do NOT explicitly calculate inv(G).
        %
        %       trace(G \ I)
        %
        %  is numerically preferable to:
        %
        %       trace(inv(G))
        %
        % ------------------------------------------------------------

        Jls = real( ...
            trace( ...
                G \ eye(LestDesign) ...
            ) ...
        );


        %% -----------------------------------------------------------
        %  5.7 CONDITION NUMBER
        % ------------------------------------------------------------

        conditionNumber = cond(G);


        %% -----------------------------------------------------------
        %  5.8 APERIODIC AUTOCORRELATION
        %
        %  This is calculated on the actual M-PSK symbols, NOT on
        %  the binary sequence.
        % ------------------------------------------------------------

        r = xcorr(candidateSymbols);

        center = numPilotSymbols;

        mainPeak = abs(r(center));

        % All non-zero-delay sidelobes.
        sidelobes = [ ...
            r(1:center-1);
            r(center+1:end) ...
        ];

        % Normalize with respect to zero-delay correlation.
        normalizedSidelobes = ...
            abs(sidelobes) / mainPeak;

        maxSidelobe = max(normalizedSidelobes);


        %% -----------------------------------------------------------
        %  5.9 DECIDE WHETHER THIS CANDIDATE IS BETTER
        %
        %  Lexicographic optimization:
        %
        %  1. minimum J_LS
        %  2. minimum cond(A'*A)
        %  3. minimum autocorrelation sidelobe
        %
        % ------------------------------------------------------------

        candidateIsBetter = false;


        % ---- Primary criterion: LS quality -------------------------

        if Jls < bestJls * (1 - relativeTolerance)

            candidateIsBetter = true;


        % ---- J_LS approximately equal ------------------------------

        elseif abs(Jls - bestJls) <= ...
                relativeTolerance * max(1, abs(bestJls))


            % Secondary criterion: condition number

            if conditionNumber < ...
                    bestConditionNumber * (1 - relativeTolerance)

                candidateIsBetter = true;


            % Condition numbers approximately equal

            elseif abs(conditionNumber - bestConditionNumber) <= ...
                    relativeTolerance * ...
                    max(1, abs(bestConditionNumber))


                % Tertiary criterion:
                % synchronization autocorrelation

                if maxSidelobe < bestMaxSidelobe

                    candidateIsBetter = true;

                end

            end

        end


        %% -----------------------------------------------------------
        %  5.10 SAVE BEST CANDIDATE
        % ------------------------------------------------------------

        if candidateIsBetter

            bestBits = candidateBits;

            bestSymbols = candidateSymbols;

            bestA = A;

            bestJls = Jls;

            bestConditionNumber = conditionNumber;

            bestMaxSidelobe = maxSidelobe;

        end

    end


    %% ===============================================================
    %  6. CHECK THAT SEARCH SUCCEEDED
    % ================================================================

    if isempty(bestBits)

        rng(previousRngState);

        error('No numerically valid preamble candidate was found.');

    end


    %% ===============================================================
    %  7. SPLIT OPTIMIZED 96-BIT PREAMBLE INTO 48 + 48
    %
    %  This maintains complete compatibility with the current
    %  transmitter.m and receiver.m interfaces.
    % ================================================================

    syncReference = bestBits(1:48);

    phaseReference = bestBits(49:96);


    %% ===============================================================
    %  8. FINAL DIAGNOSTICS
    % ================================================================

    Gbest = bestA' * bestA;


    % ---------------------------------------------------------------
    % Full normalized autocorrelation
    % ---------------------------------------------------------------

    [rBest, lagsBest] = xcorr(bestSymbols);

    center = numPilotSymbols;

    rBestNormalized = ...
        rBest / abs(rBest(center));


    % ---------------------------------------------------------------
    % Autocorrelation relevant specifically to LS channel memory
    %
    % delays:
    %
    %       1,...,LestDesign-1
    %
    % ---------------------------------------------------------------

    if LestDesign > 1

        lsSidelobes = ...
            rBest(center+1 : center+LestDesign-1);

        maxLSSidelobe = ...
            max(abs(lsSidelobes)) / abs(rBest(center));

    else

        maxLSSidelobe = 0;

    end


    %% ===============================================================
    %  9. OUTPUT INFORMATION STRUCTURE
    % ================================================================

    info.bits = bestBits;

    info.symbols = bestSymbols;

    info.A = bestA;

    info.gramMatrix = Gbest;

    info.Jls = bestJls;

    info.conditionNumber = bestConditionNumber;

    info.maxSyncSidelobe = bestMaxSidelobe;

    info.maxLSSidelobe = maxLSSidelobe;

    info.autocorrelation = rBestNormalized;

    info.autocorrelationLags = lagsBest;

    info.numPreambleBits = numPreambleBits;

    info.numPilotSymbols = numPilotSymbols;

    info.LestDesign = LestDesign;

    info.numTrials = numTrials;


    %% ===============================================================
    %  10. RESTORE RNG
    % ================================================================

    rng(previousRngState);


    %% ===============================================================
    %  11. DISPLAY RESULTS
    % ================================================================

    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('       OPTIMIZED LINK-LEVEL PREAMBLE\n');
    fprintf('====================================================\n');

    fprintf('Preamble length             : %d bits\n', ...
        numPreambleBits);

    fprintf('Modulation                  : %d-PSK\n', ...
        M);

    fprintf('Number of pilot symbols     : %d\n', ...
        numPilotSymbols);

    fprintf('Design channel length       : %d taps\n', ...
        LestDesign);

    fprintf('Monte Carlo candidates      : %d\n', ...
        numTrials);

    fprintf('\n');

    fprintf('LS metric trace((A''A)^-1)  : %.8f\n', ...
        bestJls);

    fprintf('Condition number cond(A''A) : %.8f\n', ...
        bestConditionNumber);

    fprintf('Max sync sidelobe           : %.6f\n', ...
        bestMaxSidelobe);

    fprintf('Max LS-range sidelobe       : %.6f\n', ...
        maxLSSidelobe);

    fprintf('\n');

    fprintf('syncReference:\n');
    disp(syncReference.');

    fprintf('phaseReference:\n');
    disp(phaseReference.');

    fprintf('====================================================\n\n');

end