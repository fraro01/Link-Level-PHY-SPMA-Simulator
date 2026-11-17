%% BASELINE A: Validation of uncoded Gray 8-PSK over AWGN (fixed-sample Monte Carlo)
% Reference: Proakis & Salehi, Digital Communications, 5th ed.
% Mapping convention: left-msb, phase offset pi/M, Gray PSK.
% Communications Toolbox required: pskmod, pskdemod, berawgn.
% Run: validate_8PSK_AWGN_v2
clear; 
clc; 
close all;

rng(42,'twister');

M = 8;
k = log2(M);
phaseOffset = pi/M;
EbN0_dB = 0:2:14;
% Sample sizes FIXED BEFORE observing any errors, separately for each SNR.
% Higher Eb/N0 requires more bits to obtain useful confidence intervals.
bitsPerPoint = [600000 600000 600000 600000 600000 600000 6000000 60000000];
bitsPerBlock = 300000; % divisible by k; bounded RAM consumption
assert(all(mod(bitsPerPoint,k)==0) && mod(bitsPerBlock,k)==0);

berTheory = berawgn(EbN0_dB,'psk',M,'nondiff');
bitErrors = zeros(size(EbN0_dB));
bitsTested = zeros(size(EbN0_dB));

for ii = 1:numel(EbN0_dB)
    gammaB = 10^(EbN0_dB(ii)/10);
    % Es=1; Eb=Es/log2(M); complex AWGN variance E|n|^2=N0.
    N0 = 1/(k*gammaB);
    nBlocks = bitsPerPoint(ii)/bitsPerBlock;
    for ib = 1:nBlocks
        txBits = randi([0 1],bitsPerBlock,1);
        txGroups = reshape(txBits,k,[]).';
        txIndices = bi2de(txGroups,'left-msb');
        txSymbols = pskmod(txIndices,M,phaseOffset,'gray');
        noise = sqrt(N0/2)*(randn(size(txSymbols))+1j*randn(size(txSymbols)));
        rxIndices = pskdemod(txSymbols+noise,M,phaseOffset,'gray');
        rxGroups = de2bi(rxIndices,k,'left-msb');
        rxBits = reshape(rxGroups.',[],1);
        bitErrors(ii) = bitErrors(ii)+nnz(txBits~=rxBits);
    end
    bitsTested(ii) = bitsPerPoint(ii);
    fprintf('Eb/N0 = %2g dB | errors = %6d | bits = %9d\n', ...
        EbN0_dB(ii),bitErrors(ii),bitsTested(ii));
end
berMeasured = bitErrors./bitsTested;

% Two-sided 95% Wilson confidence intervals for fixed-n binomial sampling.
z = 1.95996398454005;
p = berMeasured;
n = bitsTested;
den = 1+z^2./n;
center = (p+z^2./(2*n))./den;
half = z*sqrt(p.*(1-p)./n+z^2./(4*n.^2))./den;
ciLow = max(0,center-half);
ciHigh = min(1,center+half);
withinCI = (berTheory>=ciLow) & (berTheory<=ciHigh);

T = table(EbN0_dB(:),berTheory(:),berMeasured(:),bitErrors(:), ...
    bitsTested(:),ciLow(:),ciHigh(:),withinCI(:), ...
    'VariableNames',{'EbN0_dB','BER_theory','BER_measured', ...
    'BitErrors','BitsTested','CI95_Lower','CI95_Upper','TheoryInsideCI95'});
disp(T);
writetable(T,'baseline_A_validation_8PSK_AWGN_results.csv');

fig = figure('Color','w','Position',[100 100 980 610]);
ax = axes(fig);
hTheory = semilogy(ax,EbN0_dB,berTheory,'k-','LineWidth',1.8); hold(ax,'on');
valid = bitErrors>0;
hMonte = semilogy(ax,EbN0_dB(valid),berMeasured(valid),'bo', ...
    'LineWidth',1.2,'MarkerFaceColor','b','MarkerSize',6);
% Confidence intervals drawn as vertical bars on the log-BER axis.
for ii = find(valid)
    x = EbN0_dB(ii);
    line(ax,[x x],[max(ciLow(ii),realmin) ciHigh(ii)], ...
        'Color',[0.12 0.35 0.82],'LineWidth',1.1,'HandleVisibility','off');
    cap = 0.12;
    line(ax,[x-cap x+cap],[max(ciLow(ii),realmin) max(ciLow(ii),realmin)], ...
        'Color',[0.12 0.35 0.82],'LineWidth',1.1,'HandleVisibility','off');
    line(ax,[x-cap x+cap],[ciHigh(ii) ciHigh(ii)], ...
        'Color',[0.12 0.35 0.82],'LineWidth',1.1,'HandleVisibility','off');
end
zero = ~valid;
if any(zero)
    hZero = semilogy(ax,EbN0_dB(zero),ciHigh(zero),'rv', ...
        'MarkerSize',7,'LineWidth',1.3);
end
xlabel(ax,'E_b/N_0 [dB]'); ylabel(ax,'Bit Error Rate (BER)');
title(ax,'Uncoded Gray 8-PSK in AWGN: theoretical vs simulated (95% CI)');
grid(ax,'on');
set(ax,'YScale','log');
if any(zero)
    legend(ax,[hTheory hMonte hZero], ...
        {'Theoretical BER (berawgn)','Monte Carlo BER (95% Wilson CI)', ...
        '95% upper bound (zero errors)'},'Location','southwest');
else
    legend(ax,[hTheory hMonte], ...
        {'Theoretical BER (berawgn)','Monte Carlo BER (95% Wilson CI)'}, ...
        'Location','southwest');
end
% Suppress toolbar icon in the exported image where supported.
if isprop(ax,'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
exportgraphics(fig,'baseline_A_validation_8PSK_AWGN.png','Resolution',300);
fprintf('Theory inside 95%% Wilson CI at %d/%d SNR points.\n', ...
    nnz(withinCI),numel(withinCI));
fprintf('Note: 95%% pointwise CI coverage is probabilistic, not guaranteed.\n');
