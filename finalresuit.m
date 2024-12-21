bit = 4800000; %%送るデータ量
vbit = bit*2; %%畳み込み後のデータ量
y = vbit/48; 
Slength = y/2;%%サブキャリア１つあたりのデータ量
Subcarrier = 48; %%サブキャリア数
FFTsize = 64; %FFTサイズ
SI = 312500; %%帯域幅
snr_dB = 100; % SNRを指定（db)
attenuation_dB = 68; % 減衰量を指定(dB)
tbl=32; %ビタビ符号化率
Qnoise = 10^(snr_dB/10); %ここで真値にSNR_dBを真値に直しているよ 復調するときに発生するノイズ

%%送信側
%%PN符号出力(組み込み関数）
pnSequence = comm.PNSequence("Polynomial",[21 19 0],'InitialConditions',1,'SamplesPerFrame',bit);
PN = pnSequence();

%%畳み込み符号化
x = PN; %%PN符号を代入
n = x(1:bit); %%出力される畳み込み符号の数だけ配列を準備
%%畳み込み符号の組み込み関数
trellis = poly2trellis(7,[171 133]);
TT = convenc(n,trellis);

%%S/P変換
k = 1;
SPdata = zeros(y,Subcarrier); %%中身が0のS/P後のデータを入れる配列を作る
%%S/Pを行う２重ループ
for SPloop = 1:y
    for j = 1:Subcarrier
        SPdata(SPloop,j) = TT(k);
        k = k + 1;
    end
end


%%QPSK変調
z = SPdata; %%S/Pデータを代入
%%QPSKの組み込み関数
QPSKdata = zeros(Slength,Subcarrier);%用意したデータ長でサブキャリア数の配列を用意する
qpskmod = comm.QPSKModulator('PhaseOffset',pi/4,'BitInput',1);

%上で準備した配列に入れていく
for Qloop = 1:Subcarrier
    QPSKdata(:,Qloop) = qpskmod(z(:,Qloop));
end

%%OFDM変調
hx = QPSKdata.'; %%行列変換

%%OFDM変調の組み込み関数　(FFT→GI挿入→P/S変換)
ofdmmod = comm.OFDMModulator("FFTLength",FFTsize,"InsertDCNull",1,"CyclicPrefixLength",0,"PilotInputPort",1);

PilotS = [1 1 1 1]; %%パイロットサブキャリア
q = PilotS.'; %%パイロットサブキャリアを行列変換

OFDMdata = zeros(Subcarrier,1);
for Oloop = 1:Slength
    %%64=FFTポイント数+GI 僕の場合はGIなしなので64 Wi-Fiの規格でGI挿入するのならば80
    OFDMdata(1 + 64*(Oloop - 1):64*Oloop,1) = ofdmmod(hx(1:Subcarrier,Oloop),q);
end


%衛星Bの受信機側
% AWGNを生成
signal = OFDMdata; % AWGN用の信号を生成
awgn_signal = generate_awgn(signal, snr_dB);


%OFDM復調
OFDMdata2 = zeros(Subcarrier,Slength); %%OFDM復調のデータを入れる配列を用意
deMod = comm.OFDMDemodulator("FFTLength",FFTsize,"RemoveDCCarrier",1,"CyclicPrefixLength",0,"PilotOutputPort",1);
for Zloop = 1:Slength
    %%64=FFTポイント数+GI 僕の場合はGIなしなので64 Wi-Fiの規格でGI挿入するのならば80
    OFDMdata2(:,Zloop) = deMod(awgn_signal(1 + 64*(Zloop - 1):64*Zloop,1));
end
OFDMdata3 = OFDMdata2.'; %配列が横向きになっているので縦に直してるよ

%QPSK復調
QPSKdata2 = zeros(y,Subcarrier); %QPSK復調のデータを入れる配列を用意
hDemod = comm.QPSKDemodulator('BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio','Variance',1);
for Sloop = 1:Subcarrier
    QPSKdata2(:,Sloop) = hDemod(OFDMdata3(:,Sloop));
end


%%P/S変換
l = 1;
PSdataA = zeros(1,vbit); %%中身が0のP/S後のデータを入れる配列を作る
%S/Pの逆のことをしているよ
for PSloopA = 1:y
    for  PSloopB = 1:Subcarrier
        PSdataA(l) = QPSKdata2(PSloopA,PSloopB);
        l = l+1;
    end
end

PSdataB = PSdataA.';%配列の縦と横が逆になっているから直した
partitionPoints = (-2.1:0.7:2.1)/0.96;
quantizedValue1 = quantiz(-PSdataB,partitionPoints);

%ビタビ復号をしているよ
decSoft = comm.ViterbiDecoder(trellis,'InputFormat','Soft', ...
    'SoftInputWordLength',3,'TracebackDepth',tbl);
PNdata2 = decSoft(double(quantizedValue1));

%BER特性の出力
errSoft = comm.ErrorRate('ReceiveDelay',tbl);
berSoft = errSoft(PN,PNdata2);

%衛星Bの送信側(再変調)
%%畳み込み符号化
PNdataB = zeros(bit,1);
for PNsatelliteB = 33:bit
    PNdataB(PNsatelliteB-32,1) = PNdata2(PNsatelliteB,1);
end
x_satelliteB =  PNdataB; %%PN符号を代入
n_satelliteB = x_satelliteB(1:bit); %%出力される畳み込み符号の数だけ配列を準備
%%畳み込み符号の組み込み関数
trellis = poly2trellis(7,[171 133]);
TT_satelliteB = convenc(n_satelliteB,trellis);

%S/P変換
k_satelliteB = 1;
SPsatelliteB = zeros(y,Subcarrier); %%中身が0のS/P後のデータを入れる配列を作る
%%S/Pを行う２重ループ
for SPloop_satelliteB = 1:y
    for j_satelliteB = 1:Subcarrier
        SPsatelliteB(SPloop_satelliteB,j_satelliteB) = TT_satelliteB(k_satelliteB);
        k_satelliteB = k_satelliteB + 1;
    end
end

%%QPSK変調
z_satelliteB = SPsatelliteB; %%S/Pデータを代入
%%QPSKの組み込み関数
QPSK_satelliteB = zeros(Slength,Subcarrier);%用意したデータ長でサブキャリア数の配列を用意する
qpskmod = comm.QPSKModulator('PhaseOffset',pi/4,'BitInput',1);
%上で準備した配列に入れていく
for Qloop_satelliteB = 1:Subcarrier
    QPSK_satelliteB(:,Qloop_satelliteB) = qpskmod(z_satelliteB(:,Qloop_satelliteB));
end

%%OFDM変調
hx_satelliteB = QPSK_satelliteB.'; %%行列変換
%%OFDM変調の組み込み関数　(FFT→GI挿入→P/S変換)
ofdmmod = comm.OFDMModulator("FFTLength",FFTsize,"InsertDCNull",1,"CyclicPrefixLength",0,"PilotInputPort",1);
OFDM_satelliteB = zeros(Subcarrier,1);
for Oloop_satelliteB = 1:Slength
    %%64=FFTポイント数+GI 僕の場合はGIなしなので64 Wi-Fiの規格でGI挿入するのならば80
    OFDM_satelliteB(1 + 64*(Oloop_satelliteB - 1):64*Oloop_satelliteB,1) = ofdmmod(hx_satelliteB(1:Subcarrier,Oloop_satelliteB),q);
end

%角度オフセット
%どれだけ信号のレベルが下がるかを計算しているよ
angleOffset = 10^(0/20);
%信号のレベルを下げているよ
angleOffset_data = OFDM_satelliteB * angleOffset;
%%スペクトルアナライザーに出力
%scope = spectrumAnalyzer("SampleRate",(FFTsize*SI),"SpectrogramChannel",2);
%scope(OFDMdata,angleOffset_data);
%衛星Cの受信機側(再復調）
% AWGNを生成
signal_satelliteB = OFDM_satelliteB; % AWGN用の信号を生成
awgn_satelliteB = generate_awgnB(signal_satelliteB,angleOffset_data, snr_dB);

%OFDM復調
OFDM_satelliteC_2 = zeros(Subcarrier,Slength); %%OFDM復調のデータを入れる配列を用意
deMod = comm.OFDMDemodulator("FFTLength",FFTsize,"RemoveDCCarrier",1,"CyclicPrefixLength",0,"PilotOutputPort",1);
for Zloop_satelliteC = 1:Slength
    %%64=FFTポイント数+GI 僕の場合はGIなしなので64 Wi-Fiの規格でGI挿入するのならば80
    OFDM_satelliteC_2(:,Zloop_satelliteC) = deMod(awgn_satelliteB(1 + 64*(Zloop_satelliteC - 1):64*Zloop_satelliteC,1));
end
OFDM_satelliteC_3 = OFDM_satelliteC_2.'; %配列が横向きになっているので縦に直してるよ

%QPSK復調
QPSK_satelliteC_2 = zeros(y,Subcarrier); %QPSK復調のデータを入れる配列を用意
hDemod = comm.QPSKDemodulator('BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio','Variance',1);
for Sloop_satelliteC = 1:Subcarrier
    QPSK_satelliteC_2(:,Sloop_satelliteC) = hDemod(OFDM_satelliteC_3(:,Sloop_satelliteC));
end

%%P/S変換
l_satelliteC = 1;
PS_satelliteC = zeros(1,vbit); %%中身が0のP/S後のデータを入れる配列を作る
%S/Pの逆のことをしているよ
for PSloop_satelliteC = 1:y
    for  PSloop_satelliteC_A = 1:Subcarrier
        PS_satelliteC(l_satelliteC) = QPSK_satelliteC_2(PSloop_satelliteC,PSloop_satelliteC_A);
        l_satelliteC = l_satelliteC+1;
    end
end

PS_satelliteC_A = PS_satelliteC.';%配列の縦と横が逆になっているから直した

quantizedValue2 = quantiz(-PS_satelliteC_A,partitionPoints);

%
PS_satelliteC_B = zeros(vbit-64,1);
for PNsatelliteC = 1:vbit-64
    PS_satelliteC_B(PNsatelliteC,1) = quantizedValue2(PNsatelliteC,1);
end

%ビタビ復号をしているよ
decSoft = comm.ViterbiDecoder(trellis,'InputFormat','Soft', ...
    'SoftInputWordLength',3,'TracebackDepth',tbl);
PNdata_satelliteC = decSoft(double(PS_satelliteC_B));

PN_cut = zeros(bit-32,1);
for PNLoop = 1:bit-32
    PN_cut(PNLoop,1) = PN(PNLoop,1);
end

%BER特性の出力
errSoft2 = comm.ErrorRate('ReceiveDelay',tbl);
berSoft2 = errSoft(PN_cut,PNdata_satelliteC);

%オーバーリーチ信号
% 減衰された信号を生成
input_signal = OFDMdata; % オーバーリーチ信号用の信号を生成
attenuated_signal = attenuate_signal(input_signal, attenuation_dB);

%%スペクトルアナライザーに出力
%scope = spectrumAnalyzer("SampleRate",(FFTsize*SI),"SpectrogramChannel",2);
%scope(OFDMdata,attenuated_signal);

% オーバーリーチ信号のAWGNを生成
OLsignal = attenuated_signal; % AWGN用の信号を生成
OLawgn_signal = OLgenerate_awgn(OLsignal, snr_dB);

%OFDM復調（オーバーリーチ信号）
OFDMdata_OL = zeros(Subcarrier,Slength); %%OFDM復調のデータを入れる配列を用意
deMod = comm.OFDMDemodulator("FFTLength",FFTsize,"RemoveDCCarrier",1,"CyclicPrefixLength",0,"PilotOutputPort",1);
for Zloop_OL = 1:Slength
    %%64=FFTポイント数+GI 僕の場合はGIなしなので64 Wi-Fiの規格でGI挿入するのならば80
    OFDMdata_OL(:,Zloop_OL) = deMod(OLawgn_signal(1 + 64*(Zloop_OL - 1):64*Zloop_OL,1));
end
OL_OFDMdata = OFDMdata_OL.'; %配列が横向きになっているので縦に直してるよ

%QPSK復調
QPSKdata_OL = zeros(y,Subcarrier); %QPSK復調のデータを入れる配列を用意
hDemod = comm.QPSKDemodulator('BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio','Variance',1);
for Sloop_OL = 1:Subcarrier
    QPSKdata_OL(:,Sloop_OL) = hDemod(OL_OFDMdata(:,Sloop_OL));
end

%%P/S変換
l_OL = 1;
PSdata_OL = zeros(1,vbit); %%中身が0のP/S後のデータを入れる配列を作る
%S/Pの逆のことをしているよ
for PSloop_OL = 1:y
    for  OL_PSloop = 1:Subcarrier
        PSdata_OL(l_OL) = QPSKdata_OL(PSloop_OL,OL_PSloop);
        l_OL = l_OL+1;
    end
end
OL_PSdata = PSdata_OL.';%配列の縦と横が逆になっているから直した
partitionPoints2 = (-2.1:0.7:2.1)/4180.103;%減衰量ごとに変える
quantizedValue3 = quantiz(-OL_PSdata,partitionPoints2);

%
PS_OL = zeros(vbit-64,1);
for PN_OLloop = 1:vbit-64
    PS_OL(PN_OLloop,1) = quantizedValue3(PN_OLloop,1);
end

%ビタビ復号をしているよ
decSoft = comm.ViterbiDecoder(trellis,'InputFormat','Soft', ...
    'SoftInputWordLength',3,'TracebackDepth',tbl);
PNdata_OL= decSoft(double(PS_OL));

%BER特性の出力
errSoft3 = comm.ErrorRate('ReceiveDelay',tbl);
berSoft3 = errSoft(PN_cut,PNdata_OL);


%ダイバーシチ
DS_PSdata = OL_PSdata + PS_satelliteC_A;
partitionPoints3 = (-2.1:0.7:2.1)/1.05;%減衰量ごとに変える
quantizedValue4 = quantiz(-DS_PSdata,partitionPoints3);

%
PS_DS = zeros(vbit-64,1);
for PN_DSloop = 1:vbit-64
    PS_DS(PN_DSloop,1) = quantizedValue4(PN_DSloop,1);
end

%ビタビ復号をしているよ
decSoft = comm.ViterbiDecoder(trellis,'InputFormat','Soft', ...
    'SoftInputWordLength',3,'TracebackDepth',tbl);
PNdata_DS= decSoft(double(PS_DS));

%BER特性の出力
errSoft4 = comm.ErrorRate('ReceiveDelay',tbl);
berSoft4 = errSoft(PN_cut,PNdata_DS);



%使用した関数
function awgn_signal = generate_awgn(signal, snr_dB)
    % signal: 入力信号
    % snr_dB: 信号対雑音比（dB単位）
    %Offsetsignal_data:ノイズをかける信号
    % 信号のパワーを計算
    signal_power = sum(abs(signal).^2) / length(signal);
    % 信号対雑音比（SNR）を計算
    snr_linear = 10^(snr_dB / 10);
    % 雑音のパワーを計算
    noise_power = (signal_power / snr_linear)/2;
    % 雑音を生成（ガウス分布に従う）
    noise = sqrt(noise_power) * (randn(size(signal)) + 1i*(randn(size(signal))));
    % 入力信号と雑音を加算してAWGN信号を生成
    awgn_signal = signal + noise;
end

function awgn_satelliteB = generate_awgnB(signal_satelliteB, angleOffset_data, snr_dB)
    % signal: 入力信号
    % snr_dB: 信号対雑音比（dB単位）
    % 信号のパワーを計算
    signal_power = sum(abs(signal_satelliteB).^2) / length(signal_satelliteB);
    % 信号対雑音比（SNR）を計算
    snr_linear = 10^(snr_dB / 10);
    % 雑音のパワーを計算
    noise_power = (signal_power / snr_linear)/2;
    % 雑音を生成（ガウス分布に従う）
    noise = sqrt(noise_power) * (randn(size(signal_satelliteB)) + 1i*(randn(size(signal_satelliteB))));
    % 入力信号と雑音を加算してAWGN信号を生成
    awgn_satelliteB = angleOffset_data + noise;
    
end

function OLawgn_signal = OLgenerate_awgn(OLsignal, snr_dB)
    % signal: 入力信号
    % snr_dB: 信号対雑音比（dB単位）
    % 信号のパワーを計算
    signal_power2 = sum(abs(OLsignal).^2) / length(OLsignal);
    % 信号対雑音比（SNR）を計算
    snr_linear2 = 10^(snr_dB / 10);
    % 雑音のパワーを計算
    noise_power2 = (signal_power2 / snr_linear2)/2;
    % 雑音を生成（ガウス分布に従う）
    noise2 = sqrt(noise_power2) * (randn(size(OLsignal)) + 1i*(randn(size(OLsignal))));
    % 入力信号と雑音を加算してAWGN信号を生成
    OLawgn_signal =  OLsignal+ noise2;
end

%減衰器
function attenuated_signal = attenuate_signal(input_signal, attenuation_dB)
    % input_signal: 入力信号
    % attenuation_dB: 減衰量（dB単位）
    % 減衰量を線形スケールに変換
    attenuation_linear = 10^(-attenuation_dB / 20);
    % 入力信号を減衰
    attenuated_signal = attenuation_linear * input_signal;
end