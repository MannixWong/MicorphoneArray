%�����γ�ʵ�ʲ��Գ���
%16��Ԫ��36��Ԫ�л�ֻ���滻�����ļ�����(E16/E36)
%�������ĸ߶�1.2845m
%��Դ�߶�1.27
%�˳���ΪԶ��
clear;close all;clc
load('D:\Programming\Matlab\MicrophoneArray\E36.mat')
%ȫ�ֱ�������
global C N M FMAIN Array_X Array_Y
C=340;%sound speed
% $c=\sqrt(1.4*287.06*T)$TΪ����������ѧ�¶ȣ���λΪK
L=0.82159;%���п׾���֮����ʵ�ʲ���һ��
r=5.6;%��Դ����˷����еľ���
Alogrithm=1;%1:Beamforming 2:MUSIC
%%

Fs=44100;%����Ƶ��
N=length(X(1,:));%������
M=length(X(:,1));%��Ԫ��
disp(['��Ԫ��Ϊ: ',num2str(M)])
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
% legend('ͨ��1','ͨ��16')
grid on
ax=gca;
ax.LineWidth=0.5;%���������߿�Ϊ0.5��
ax.XAxis.LineWidth=1;%����x�������߿�Ϊ1��
ax.YAxis.LineWidth=1;%����y�������߿�Ϊ1��
for i=1:M
    delay(i)=finddelay(X(1,:),X(i,:))/Fs;%���������ʱ��λ��s
end
[a,b]=xcorr(X(1,:),X(7,:),'coeff');
figure
plot(b,a)

X1=fft(X(1,:));
figure
ff=(0:N/2-1)*Fs/N;
F=abs(X1(1:N/2));
plot(ff,F,'Color',[0 0 0]);
xlabel('Ƶ��/Hz','FontName','����');ylabel('��ֵ','FontName','����');
xlim([0 10000]);
title('�ź�Ƶ��','FontName','����')
[maxFvalue,maxFindex]=max(F(:));
FMAIN=ff(maxFindex);%�ź���Ƶ
lambda=C/FMAIN;%����
disp(['�źŵ���ƵΪ: ',num2str(FMAIN),'Hz'])
if r<=2*L^2/lambda
    disp('����')
else
    disp('Զ��')
end
%% ������ʱ����
% DESCRIPTIVE TEXT
ArrayPosition=[0 0 0];
SpeakerPosition=[2.5 -0.0145 3];%yΪ�������ĸ߶ȼ����������ĸ߶�
r2central=sqrt(sum((SpeakerPosition-ArrayPosition).^2));%��Դ���������ľ���
ele_d=asind(SpeakerPosition(3)/r2central);%���۸����ǣ���XOYƽ��н�
azi_d=atand(SpeakerPosition(2)/SpeakerPosition(1));%���۷�λ�ǣ���X��н�
for m=1:M
    tau(m)=1/C*(Array_X(m)*cosd(azi_d)*cosd(ele_d)+Array_Y(m)*sind(azi_d)*cosd(ele_d));%��������ʱ
end
figure
error=delay-tau;
m=linspace(1,M,M);
plot(m,error,'-o')
y_val=get(gca,'YTick');   %Ϊ�˻��y����
y_str=num2str(y_val');    %Ϊ�˽�����ת��Ϊ�ַ�����
set(gca,'YTickLabel',y_str);    %��ʾ
xlim([1 M])
xlabel('��Ԫ���','FontName','����');ylabel('��ʱ���','FontName','����');
title('������ʱ��ʵ����ʱ���','FontName','����')
%% ��Դ��λ�㷨
% DESCRIPTIVE TEXT

P=DOAalogrithm(X,Alogrithm);%��λ�㷨


[azi_max,ele_max]=find(P==max(P(:)));
disp(['Ԥ�ⷽλ��: ',num2str(azi_max),'��','Ԥ�⸩����: ',num2str(ele_max),'��'])
% P_CB=20*log10(P_CB/max(P_CB(:)));
figure
h=pcolor(P);%�ȸ���ͼ
xlabel('������/(\circ)','FontName','����');ylabel('��λ��/(\circ)','FontName','����')
set(h,'edgecolor','none','facecolor','interp');%ȥ������ƽ������
colorbar %���ɫ��
map=jet(256);%��չ��256ɫ����ɫ���޷ֲ�
colormap(map);
hold on
plot(ele_max,azi_max,'g*')
text(ele_max,azi_max,['  ele=',num2str(ele_max),newline,'  azi=',num2str(azi_max),newline],'Color','g');
% title('���������γɵ���Դ4000Hz','FontName','����')
titlestring=strcat('��������Դ',num2str(FMAIN),'Hz');
title(titlestring,'fontname','����','FontSize',12);
ax=gca;
ax.LineWidth=0.5;%���������߿�Ϊ0.5��
ax.XAxis.LineWidth=1;%����x�������߿�Ϊ1��
ax.YAxis.LineWidth=1;%����y�������߿�Ϊ1��

figure%��άͼ
mesh(P,'FaceColor','interp');%3D figure
xlabel('������/(\circ)','FontName','����');ylabel('��λ��/(\circ)','FontName','����');zlabel('��ֵ/dB','FontName','����');
colorbar %���ɫ��
map=jet(256);%��չ��256ɫ����ɫ���޷ֲ�
colormap(map);
titlestring=strcat('��������Դ',num2str(FMAIN),'Hz');
title(titlestring,'fontname','����','FontSize',12);
ax=gca;
ax.LineWidth=0.5;%���������߿�Ϊ0.5��
ax.XAxis.LineWidth=1;%����x�������߿�Ϊ1��
ax.YAxis.LineWidth=1;%����y�������߿�Ϊ1��