function [s4_t_tm,Sac_RD]=RD(srnm,V,lambda,PRF,Fs,Kr,R0,squint_angle)
%% 函数功能：时域校正距离走动，频域校正距离弯曲的距离多普勒算法
%%
c=3e8;           %光速（m/s）
fc=c/lambda;     %载频（Hz）
[Nfast,Nslow]=size(srnm);
Ta=Nslow/PRF;
fa=([0:Nslow-1]-ceil(Nslow/2))/Nslow*PRF;    %多普勒频率
fam=2*V*cos(squint_angle)/lambda;          
fr=([0:Nfast-1]'-ceil(Nfast/2))/Nfast*Fs;   % 距离向频率
DeltaR=c/(2*Fs);                            % SAR图像的距离分辨率（m）
rs=([0:Nfast-1]'-ceil(Nfast/2))*DeltaR;     % SAR图像对应的距离门矢量（m）
fd_center=2*V*sin(squint_angle)/lambda;     % 多普勒中心频率（Hz）
tm=([0:Nslow-1]-ceil(Nslow/2))/PRF;
%% 距离FFT后矫正距离走动
H1=exp(-1j*4*pi*V*sin(squint_angle)*(fr+fc)*tm/c);
S1_fr_tm=fftshift(fft(srnm,[],1),1).*H1;
%% 方位FFT变换到二维频域
S1_fr_fa=fftshift(fft(S1_fr_tm,[],2),2);
sintheta=fa/fam;              % 慢时间对应的斜视角正弦矢量
costheta=sqrt(1-sintheta.^2); 
%% 距离压缩和距离徙动校正==>距离IFFT==>距离时域-多普勒域
H21=exp(1j*pi*fr.^2*(1/Kr-2*lambda*R0.*sintheta.^2./(c.^2.*(costheta).^3))); % 距离向脉冲压缩频率函数
H22=(exp(1j*2*pi*R0.*fr*sintheta.^2/c));                                     % 距离徙动校正频率函数(0，2π)
S2_fr_fa=S1_fr_fa.*H21.*H22;          
S2_t_fa=(ifft(ifftshift(S2_fr_fa,1),[],1));
%% RD方位压缩
H3=exp(1j*2*pi*(R0+rs)*sqrt(fam.^2-...
     fa.^2)/V/cos(squint_angle)); % 方位向脉冲压缩频率函数 
S3_t_fa=S2_t_fa.*H3;clear H3;
s4_t_tm=ifft(ifftshift(S3_t_fa,2),[],2);
