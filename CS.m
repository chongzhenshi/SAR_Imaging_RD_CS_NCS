function [s_Imaging_By_CS]=CS(Signal,v,squint_angle,lambda,PRF,Doppler_centroid,Fs,Kr,R0_SceneCenter)
%% 参数设定
[Nrn,Nan]=size(Signal);
% Nrn_mid = fix(Nrn/2);               %--场景中心对应的距离点数 
% Nan_mid = fix(Nan/2);               %--方位向中心点数
c=3e8;                                                            
Range_Res=c/2/Fs;                   %--距离分辨率（m）
PRT=1/PRF;                          %--脉冲重复周期（s） 
fc=c/lambda;                        %--载频（Hz）
% B = abs(Kr*Tp);                     %--发射信号带宽
% deltaR = c/2/Fs;                    %--距离分辨率
Rs = R0_SceneCenter;
% Rs = R0 + deltaR* Nrn_mid;          %--整个幅宽场景中心斜距（前提是R0是场景最近斜距）
% fr_step = Fs/Nrn;                   %--距离频率分辨率
% fa_step = PRF/Nan;                  %--方位频率分辨率

tm=([0:Nan-1]-ceil(Nan/2))*PRT;     %--慢时间变量
fa=([0:Nan-1]-ceil(Nan/2))/Nan*PRF; %--多普勒频率变量
fam=2.*v/lambda.*cos(squint_angle);   %--载机正前方目标回波多普勒（最大多普勒）
% sintheta = fa/fam;
% costheta = sqrt(1-sintheta.^2);
fr=([0:Nrn-1]'-ceil(Nrn/2))/Nrn*Fs; %--距离频率变量   
t=2* Rs/c + [-Nrn/2:Nrn/2-1]'/Fs;   %--快时间变量
R0=R0_SceneCenter+([0:Nrn-1]'-ceil(Nrn/2))*Range_Res;   %--每个距离单元对应的斜视维距离（m）
%% 频谱搬移，对准多普勒中心
Signal = Signal.*exp(-2j*pi*Doppler_centroid*tm);
%% 利用波束中心点处的速度进行距离走动矫正
H_RMC=exp(-1j*4*pi/c.*(fr+fc).*v.*sin(squint_angle).*tm);%--距离走动校正函数
S_fr_tm_RMC=fftshift(fft(Signal,[],1),1).*H_RMC;        %--距离频域进行校正
s_t_tm_RMC=ifft(ifftshift(S_fr_tm_RMC,1),[],1);         %--再进行距离IFFT
% figure;imagesc(abs(s_t_tm_RMC));colorbar;
% xlabel('方位向');ylabel('距离向');
% title('距离脉压前距离走动校正后的数据');
%% 方位FFT==>CS相位函数==>距离FFT
tao_d_fa_R0=2/c.*R0_SceneCenter./sqrt(1-(fa./fam).^2);           %--参考距离处的时延
Kr_e=1./(1/Kr-R0_SceneCenter*2*lambda/c^2 ...
    *((fa./fam).^2)./(sqrt(1-(fa./fam).^2)).^3);                %--参考距离处新调频率
a_fa=1./sqrt(1-(fa./fam).^2)-1;                                 %--CS因子
S_t_fa_Signal=fftshift(fft(s_t_tm_RMC,[],2),2);                 %--方位FFT
clear s_t_tm_RMC;
CS=exp(1j*pi*(Kr_e.*a_fa).* ...
    (t*ones(1,Nan)-tao_d_fa_R0).^2);             %--改变调频率尺度的CS二次相位函数
S_fr_fa_CS=S_t_fa_Signal.*CS;
S_fr_fa_CS=fftshift(fft(S_fr_fa_CS,[],1),1);                  %--距离FFT
clear S_t_fa_Signal;clear CS;
% figure;imagesc(abs(S_fr_fa_CS));
% xlabel('方位频域向');ylabel('距离频域向');
% title('二维频域数据');
%% 距离压缩
H_RangeComp=exp(1j*pi.*fr.^2./((1+a_fa).*Kr_e));              %--距离压缩函数
S_fr_fa_RangeComp=S_fr_fa_CS.*H_RangeComp; 
clear H_RangeComp;
% S_t_fa_RangeComp=ifft(ifftshift(S_fr_fa_RangeComp,1),[],1);
% figure;imagesc(abs(S_t_fa_RangeComp));
% xlabel('方位向');ylabel('距离向');
% title('距离脉压后的数据');
%% 距离徙动动校正==>距离IFFT
H_match_RCMC=exp(1j*4*pi*R0_SceneCenter.*a_fa.*fr/c);                     %--距离徙动校正函数
S_fr_fa_RCMC=S_fr_fa_RangeComp.*H_match_RCMC;
clear H_match_RCMC;
S_t_fa_RCMC=ifft(ifftshift(S_fr_fa_RCMC,1),[],1);             %--距离IFFT
clear S_fr_fa_RCMC;
% figure;imagesc(abs(S_t_fa_RCMC));
% xlabel('方位向');ylabel('距离向');
% title('距离徙动校正后的数据');
%% 方位压缩
H_AziComp=exp(1j*2*pi*R0.*sqrt(fam.^2-fa.^2)./(v.*cos(squint_angle)));
S_t_fa_match=S_t_fa_RCMC.*H_AziComp;   % 这里可以对不进行残余相位补偿的数据进行输出
clear S_t_fa_RCMC
%% 残余相位补偿
H_RPC=exp(1j*4*pi/c^2*Kr_e.*a_fa.*(1+a_fa).*(R0-R0_SceneCenter).^2);
S_t_fa_RPC=S_t_fa_match.*H_RPC;
s_t_tm_RPC=ifft(ifftshift(S_t_fa_RPC,2),[],2);
s_Imaging_By_CS=s_t_tm_RPC;
% figure;imagesc(abs(s_t_tm_RPC));
% xlabel('方位向');ylabel('距离向');
% title('CS成像结果');
%% 计算峰值旁瓣比和积分旁瓣比和冲击响应宽度
tanew=linspace(t(1),t(end),100*Nrn);
y0 = interp1(tm,s_Imaging_By_CS(Nr_x,:),tanew,'spline');
figure;plot(20*log10(abs(y0./max(y0))));
%title('方位维切片');
grid on;set_type;ylabel('Normalized Amplitude(dB)');xlabel('Azimuth Gates')
tnew=linspace(t(1),t(end),100*Nrn);
y1 = interp1(t,Image_Result(:,Na_x),tnew,'spline');
figure;plot(20*log10(abs(y1./max(y1))));%title('距离维切片');
grid on;set_type;xlabel('Range Gates');ylabel('Normalized Amplitude(dB)')

