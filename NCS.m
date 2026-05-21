function [s_t_tm_RCMC_By_NCS,s_Imaging_By_NCS]=NCS(Echo,v,theta_0,lambda,PRF,Fs,Kr,R0_SceneCenter,t)
%% 参数准备
[Nrn,Nan]=size(Echo);
theta_0=theta_0/pi*180;
c=3e8;
PRT=1/PRF;                         %--脉冲重复周期（s）                                                             
tm=([0:Nan-1]-ceil(Nan/2))*PRT;    %--慢时间变量
Range_Res=c/2/Fs;                  %--距离分辨率（m）
fa=([0:Nan-1]-ceil(Nan/2))/Nan*PRF;                   %--多普勒频率变量
fr=([0:Nrn-1]'-ceil(Nrn/2))/Nrn*Fs;%--距离频率变量                                       
fc=c/lambda;                       %--载频（Hz）
%% 利用波束中心点处的速度进行距离走动矫正
echo_fr_tm_match1=fftshift(fft(Echo,[],1),1)...
    .*exp(-1j*4*pi/c*(fc+fr)*(v*sind(theta_0)*tm));
echo_RWC=ifft(ifftshift(echo_fr_tm_match1,1),[],1);
figure;imagesc(abs(echo_RWC));
title('中间过程')
%%
fdc=2*v*cosd(theta_0)/lambda*sind(theta_0);                         %--波束中心多普勒频率
R0=R0_SceneCenter+([0:Nrn-1]'-ceil(Nrn/2))*Range_Res; %--每个距离单元对应的斜视维距离（m）
faM=2*v*cosd(theta_0)/lambda;                         %--最大多普勒频率（Hz）
fa=([0:Nan-1]-ceil(Nan/2))/Nan*PRF;                   %--多普勒频率变量
sin_theta=fa./faM;
cos_theta=sqrt(1-sin_theta.^2);
Gama_e=1./(1/Kr-R0*2*lambda*(sin_theta.^2/c^2./cos_theta.^3));%--距离维等效调频率
Phi3=-2*pi*R0*lambda^2*(sin_theta.^2./(c^3.*cos_theta.^5));  
Ks=Gama_e.^2./fc.*repmat(sin_theta.^2./cos_theta.^2,Nrn,1);           
Alpha=sqrt(1-(fdc/faM)^2)./cos_theta;               
Ym=-Ks./Gama_e.^3.*repmat((Alpha-0.5)./(Alpha-1),Nrn,1);           
Y=-Ym-3/2/pi*Phi3;      
q2=Gama_e.*(repmat(Alpha,Nrn,1)-1); 
q3=-Ks.*(repmat(Alpha,Nrn,1)-1)/2;   
Taod_fa_Rs=2/c*R0_SceneCenter./cos_theta;                      %--计算taod_fa_Rs：场景中心延时         
Taod_fdc_Rs=2/c*R0_SceneCenter./sqrt(1-fdc./faM.^2);           %--计算taod_fdc_Rs：场景中心延时     
deltaTao_fa_R0=2/c*R0*(1./cos_theta)-repmat(Taod_fa_Rs,Nrn,1); %--相对于场景中心的每个距离单元对应的延时差 
%% H1三次相位滤波
H1_fr=exp(1i*2*pi/3*(fr*ones(1,Nan)).^3.*Y);
Echo_fr_fa_match1=fftshift(fft2(echo_RWC)).*H1_fr;
Echo_t_fa_match1=ifft(ifftshift(Echo_fr_fa_match1,1),[],1);
figure;imagesc(abs(Echo_fr_fa_match1));
title('中间过程')
%% H2 NCS相位函数
H2_t_fa=exp(1i*pi*q2.*(t*ones(1,Nan)-ones(Nrn,1)*Taod_fa_Rs).^2-1i*2*pi/3*q3.*(t*ones(1,Nan)-ones(Nrn,1)*Taod_fa_Rs).^3);
Echo_t_fa_match2=Echo_t_fa_match1.*H2_t_fa;
Echo_fr_fa_match2=fftshift(fft(Echo_t_fa_match2,[],1),1);
figure;imagesc(abs(Echo_fr_fa_match2));
title('NCS处理后二维频域数据');
%% 距离压缩、H3二次距离压缩以及距离徙动矫正
H3_fr_fa_R0=exp(1i*2*pi.*fr*(Taod_fa_Rs-Taod_fdc_Rs)).*exp(1i*2*pi*(q3+Ym.*Gama_e.^3).*...
    (fr*ones(1,Nan)).^3./(3*(ones(Nrn,1)*Alpha).^3.*Gama_e.^3)).*exp(1i*pi.*(fr*ones(1,Nan)).^2./((ones(Nrn,1)*Alpha).*Gama_e));
Echo_fr_fa_match3=Echo_fr_fa_match2.*H3_fr_fa_R0;
figure;imagesc(abs(Echo_fr_fa_match3));
title('中间过程')
Echo_t_fa_match3=ifft(ifftshift(Echo_fr_fa_match3,1),[],1); 
s_t_tm_RCMC_By_NCS=ifft(ifftshift(Echo_t_fa_match3,2),[],2);%--NCS进行距离脉压后的结果
%% 方位压缩处理
s_Imaging_By_NCS=ifft(ifftshift(Echo_t_fa_match3.*exp(1i*2*pi/(v*cosd(theta_0)).*R0*sqrt(faM.^2-fa.^2))...
    .*exp(-1i.*(pi*Gama_e.*(1-1./(ones(Nrn,1)*Alpha)).*deltaTao_fa_R0.^2+...
    pi*Ks/3.*(1-1./(ones(Nrn,1)*Alpha)).*deltaTao_fa_R0.^3)),2),[],2);




















