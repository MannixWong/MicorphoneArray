%波束形成实际测试程序
%16阵元与36阵元切换只需替换数据文件即可(E16/E36)
%阵列中心高度1.2845m
%声源高度1.27
%此程序为远场
clear;close all;clc
load('D:\Programming\Matlab\MicrophoneArray\E36.mat')
%全局变量声明
global C N M FMAIN Array_X Array_Y
C=340;%sound speed
% $c=\sqrt(1.4*287.06*T)$T为开尔文热力学温度，单位为K
L=0.82159;%阵列孔径，之后再实际测量一下
r=5.6;%声源距麦克风阵列的距离
Alogrithm=1;%1:Beamforming 2:MUSIC
%%

Fs=44100;%采样频率
N=length(X(1,:));%快拍数
M=length(X(:,1));%阵元数
disp(['阵元数为: ',num2str(M)])
Array_X=[0.000 -2.468 -1.510 1.535 2.458 -0.015 -12.985 -1.276 12.196 8.814 ...%the first column of excel,unit:meter!!!
    -6.749 -22.675 -8.026 17.714 18.974 -5.987 -21.711 -25.279 6.088 ...
    29.042 11.861 -12.662 -34.434 -8.619 29.107 26.609 -0.842 -37.427 ...
    -22.327 23.628 36.929 10.459 -36.729 -33.159 16.236 43.193]'*0.01;
Array_Y=[0.000 -0.786 2.104 2.086 -0.815 -2.590 2.877 13.239 5.305 -9.960 ...
    -11.461 -1.072 21.234 14.133 -12.461 -21.896 -19.526 14.615 28.558 ...
    3.035 -26.682 -32.093 2.125 33.406 18.521 -21.959 -39.090 -11.352 ...
    32.099 31.153 -12.845 -42.018 -22.932 27.845 40.141 -3.037]'*0.01;

figure
for i=1:M
    scatter(Array_X(i),Array_Y(i),'r')%2d array figure
    hold on
    text(Array_X(i),Array_Y(i),num2str(i))
end
grid on
%% SECTION TITLE
% DESCRIPTIVE TEXT

figure
plot(X(1,:));
hold on
plot(X(3,:));
% legend('通道1','通道16')
grid on
ax=gca;
ax.LineWidth=0.5;%设置网格线宽为0.5磅
ax.XAxis.LineWidth=1;%设置x坐标轴线宽为1磅
ax.YAxis.LineWidth=1;%设置y坐标轴线宽为1磅
for i=1:M
    delay(i)=finddelay(X(1,:),X(i,:))/Fs;%这里算的延时单位是s
end
[a,b]=xcorr(X(1,:),X(7,:),'coeff');
figure
plot(b,a)

X1=fft(X(1,:));
figure
ff=(0:N/2-1)*Fs/N;
F=abs(X1(1:N/2));
plot(ff,F,'Color',[0 0 0]);
xlabel('频率/Hz','FontName','宋体');ylabel('幅值','FontName','宋体');
xlim([0 10000]);
title('信号频率','FontName','黑体')
[maxFvalue,maxFindex]=max(F(:));
FMAIN=ff(maxFindex);%信号主频
lambda=C/FMAIN;%波长
disp(['信号的主频为: ',num2str(FMAIN),'Hz'])
if r<=2*L^2/lambda
    disp('近场')
else
    disp('远场')
end
%% 理论延时计算
% DESCRIPTIVE TEXT
ArrayPosition=[0 0 0];
SpeakerPosition=[2.5 -0.0145 3];%y为阵列中心高度减扬声器中心高度
r2central=sqrt(sum((SpeakerPosition-ArrayPosition).^2));%声源距阵列中心距离
ele_d=asind(SpeakerPosition(3)/r2central);%理论俯仰角，与XOY平面夹角
azi_d=atand(SpeakerPosition(2)/SpeakerPosition(1));%理论方位角，与X轴夹角
for m=1:M
    tau(m)=1/C*(Array_X(m)*cosd(azi_d)*cosd(ele_d)+Array_Y(m)*sind(azi_d)*cosd(ele_d));%求理论延时
end
figure
error=delay-tau;
m=linspace(1,M,M);
plot(m,error,'-o')
y_val=get(gca,'YTick');   %为了获得y轴句柄
y_str=num2str(y_val');    %为了将数字转换为字符数组
set(gca,'YTickLabel',y_str);    %显示
xlim([1 M])
xlabel('阵元序号','FontName','宋体');ylabel('延时误差','FontName','宋体');
title('理论延时与实际延时误差','FontName','黑体')
%% 声源定位算法
% DESCRIPTIVE TEXT

P=DOAalogrithm(X,Alogrithm);%定位算法


[azi_max,ele_max]=find(P==max(P(:)));
disp(['预测方位角: ',num2str(azi_max),'°','预测俯仰角: ',num2str(ele_max),'°'])
% P_CB=20*log10(P_CB/max(P_CB(:)));
figure
h=pcolor(P);%等高线图
xlabel('俯仰角/(\circ)','FontName','宋体');ylabel('方位角/(\circ)','FontName','宋体')
set(h,'edgecolor','none','facecolor','interp');%去掉网格，平滑网络
colorbar %添加色标
map=jet(256);%扩展到256色，颜色条无分层
colormap(map);
hold on
plot(ele_max,azi_max,'g*')
text(ele_max,azi_max,['  ele=',num2str(ele_max),newline,'  azi=',num2str(azi_max),newline],'Color','g');
% title('螺旋阵波束形成单声源4000Hz','FontName','黑体')
titlestring=strcat('螺旋阵单声源',num2str(FMAIN),'Hz');
title(titlestring,'fontname','黑体','FontSize',12);
ax=gca;
ax.LineWidth=0.5;%设置网格线宽为0.5磅
ax.XAxis.LineWidth=1;%设置x坐标轴线宽为1磅
ax.YAxis.LineWidth=1;%设置y坐标轴线宽为1磅

figure%三维图
mesh(P,'FaceColor','interp');%3D figure
xlabel('俯仰角/(\circ)','FontName','宋体');ylabel('方位角/(\circ)','FontName','宋体');zlabel('幅值/dB','FontName','宋体');
colorbar %添加色标
map=jet(256);%扩展到256色，颜色条无分层
colormap(map);
titlestring=strcat('螺旋阵单声源',num2str(FMAIN),'Hz');
title(titlestring,'fontname','黑体','FontSize',12);
ax=gca;
ax.LineWidth=0.5;%设置网格线宽为0.5磅
ax.XAxis.LineWidth=1;%设置x坐标轴线宽为1磅
ax.YAxis.LineWidth=1;%设置y坐标轴线宽为1磅