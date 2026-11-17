%% BASELINE B: Turbo-coded Gray 8-PSK, ideal AWGN receiver
% Uses the same MATLAB TurboEncoder/Decoder defaults and max-log LLR convention
% as the thesis code, but bypasses synchronization, preamble, CP and hopping.
% Eb/N0 refers to information bits, counting coded DATA symbols only.
% MATLAB Communications Toolbox required.
clear; 
clc; 

rng(42,'twister');

M=8; 
k=3; 
phase=pi/M; 
K=1331;
EbN0dB=0:0.25:2;
packetsPerPoint=1000;  % Fixed sample size at every Eb/N0
nS=numel(EbN0dB);
BER=zeros(nS,1); PER=zeros(nS,1); bits=zeros(nS,1); errors=zeros(nS,1);
packetErrors=zeros(nS,1);
for s=1:nS
    gamma=10^(EbN0dB(s)/10);
    for trial=1:packetsPerPoint
        u=randi([0 1],K,1);
        piIdx=randperm(K);
        enc=comm.TurboEncoder('InterleaverIndices',piIdx);
        dec=comm.TurboDecoder('InterleaverIndices',piIdx);
        c=enc(u); Nc=numel(c);
        pad=mod(-Nc,k);
        cp=[c;zeros(pad,1)];
        idx=bi2de(reshape(cp,k,[]).','left-msb');
        x=pskmod(idx,M,phase,'gray');
        % Energy per information bit counts the actual padded DATA symbols.
        % Es=1, Eb_info = numel(x)/K; N0=Eb_info/gamma.
        N0=numel(x)/(K*gamma);
        y=x+sqrt(N0/2)*(randn(size(x))+1j*randn(size(x)));
        llr=pskdemod(y,M,phase,'gray','OutputType','approxllr', ...
            'NoiseVariance',N0);
        % Thesis RX negates pskdemod LLR before MATLAB TurboDecoder.
        uHat=dec(-llr(1:Nc));
        e=nnz(uHat~=u);
        errors(s)=errors(s)+e; bits(s)=bits(s)+K;
        packetErrors(s)=packetErrors(s)+(e>0);
    end
    BER(s)=errors(s)/bits(s);
    PER(s)=packetErrors(s)/packetsPerPoint;
    fprintf('B: Eb/N0=%g dB, BER=%.5g, PER=%.5g (%d packets)\n', ...
        EbN0dB(s),BER(s),PER(s),packetsPerPoint);
end
% Wilson 95% CI for PER (packets are independent experimental units).
z=1.959963984540054;
n=packetsPerPoint; p=PER;
den=1+z^2/n;
center=(p+z^2/(2*n))/den;
half=z*sqrt(p.*(1-p)/n+z^2/(4*n^2))/den;
PER_Lower=max(0,center-half); PER_Upper=min(1,center+half);
% BER counts are correlated within Turbo-coded packets: do NOT attach
% naive bitwise-binomial confidence intervals to BER.
T=table(EbN0dB(:),BER,PER,errors,bits,packetErrors, ...
 PER_Lower,PER_Upper, ...
 'VariableNames',{'EbN0_dB','BER','PER','BitErrors','BitsTested', ...
 'PacketErrors','PER_CI95_Lower','PER_CI95_Upper'});
disp(T); writetable(T,'baseline_B_coded_AWGN.csv');
% Two-panel figure: BER positive measurements and PER Wilson intervals.
fig=figure('Color','w','Position',[100 100 1000 450]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile;
positive=BER>0;
semilogy(EbN0dB(positive),BER(positive),'o-','LineWidth',1.5); hold on;
zero=~positive;
if any(zero)
    % 95% upper bound from zero observed packet errors implies a
    % conservative upper bound for BER: BER <= PER (since BER <= PER).
    ub=PER_Upper(zero);
    semilogy(EbN0dB(zero),ub,'v','MarkerSize',7,'LineWidth',1.5);
    legend('Measured BER','Upper bound (BER <= PER_{95%,upper})', ...
        'Location','southwest');
else
    legend('Measured BER','Location','southwest');
end
grid on; xlabel('E_b/N_0 per information bit [dB]'); ylabel('BER');
title('Turbo-coded 8-PSK: BER');
nexttile;
errorbar(EbN0dB,PER,PER-PER_Lower,PER_Upper-PER, ...
    'o-','LineWidth',1.2,'CapSize',4);
ylim([0 1]); grid on;
xlabel('E_b/N_0 per information bit [dB]'); ylabel('PER');
title('PER and pointwise Wilson 95% CI');
try
    ax=findall(fig,'Type','axes');
    for a=1:numel(ax)
        if isprop(ax(a),'Toolbar'), ax(a).Toolbar.Visible='off'; end
    end
catch
end
exportgraphics(fig,'baseline_B_waterfall_0to2dB.png','Resolution',250);
fprintf('Saved baseline_B_coded_AWGN.csv and baseline_B_waterfall_0to2dB.png\n');
