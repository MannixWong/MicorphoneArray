%�����γ�ʵ�ʲ��Գ���
%�������ĸ߶�1.287m
%��Դ�߶�1
%�˳���ΪԶ��
clear;close all;clc
load('D:\Programming\Matlab\MicrophoneArray\E36.mat')
%ȫ�ֱ�������
global C N M FMAIN Array_X Array_Y Array_Z Central_X Central_Y Central_Z Fs lambda
global gridX gridY gridZ
C=340;%sound speed
% $c=\sqrt(1.4*287.06*T)$TΪ����������ѧ�¶ȣ���λΪK
L=0.82159;%���п׾���֮����ʵ�ʲ���һ��
Alogrithm=3;%1:Beamforming 2:MUSIC 3:Functional Beamforming
if Alogrithm==1
    alostring='CB';
elseif Alogrithm==2
    alostring='MUSIC';
elseif Alogrithm==3
    alostring='FB';
end
%%

Fs=44100;%����Ƶ��
N=length(X(1,:));%������
M=length(X(:,1));%��Ԫ��
disp(['��Ԫ��Ϊ: ',num2str(M)])
CentralPosition=[0 0 0];%�������
SpeakerPosition=[0.6 0 5];%xΪ�������0;yΪ���������ĸ߶ȼ��������ĸ߶�;zΪ�������
Array_X=[0.000 -2.468 -1.510 1.535 2.458 -0.015 -12.985 -1.276 12.196 8.814 ...%the first column of excel,unit:meter!!!
    -6.749 -22.675 -8.026 17.714 18.974 -5.987 -21.711 -25.279 6.088 ...
    29.042 11.861 -12.662 -34.434 -8.619 29.107 26.609 -0.842 -37.427 ...
    -22.327 23.628 36.929 10.459 -36.729 -33.159 16.236 43.193]'*0.01;
Array_Y=[0.000 -0.786 2.104 2.086 -0.815 -2.590 2.877 13.239 5.305 -9.960 ...
    -11.461 -1.072 21.234 14.133 -12.461 -21.896 -19.526 14.615 28.558 ...
    3.035 -26.682 -32.093 2.125 33.406 18.521 -21.959 -39.090 -11.352 ...
    32.099 31.153 -12.845 -42.018 -22.932 27.845 40.141 -3.037]'*0.01;
Array_Z=zeros(M,1);
Central_X = 0;
Central_Y = 0;
Central_Z = 0;

stepX = 0.01; % ���ƽ����������
stepY = 0.01;

gridX = (-1.5:stepX:1.5); %���ƽ�淶Χ
gridY = (-1.5:stepY:1.5);
gridZ = 5;%���ƽ��������֮��ľ���
% h5create('MicGeo36.h5','/front/myDataset',[36 1 1])
% h5write('MicGeo36.h5','/front/myDataset',Array_X)
% h5disp('MicGeo36.h5')
% h5read('MicGeo36.h5','/front/myDataset')
% figure
% for i=1:M
%     scatter(Array_X(i),Array_Y(i),'r')%2d array figure
%     hold on
%     text(Array_X(i),Array_Y(i),num2str(i))
% end
% grid on
%������Ԫ�Ų�ģʽ����Դ����ļ���ʾ��
figure
scatter3(Array_X,Array_Y,Array_Z,'r')%3d array figure
hold on
scatter3(SpeakerPosition(1),SpeakerPosition(2),SpeakerPosition(3),20,'gp')%sound source figure
hold on
plot3([SpeakerPosition(1),0],[SpeakerPosition(2),0],[SpeakerPosition(3),0])%connect the sound source point & central point
for i=1:M
    text(Array_X(i),Array_Y(i),num2str(i))
end

%% SECTION TITLE
% DESCRIPTIVE TEXT

% figure
% plot(X(1,:));
% hold on
% plot(X(3,:));
% % legend('ͨ��1','ͨ��16')
% grid on
% ax=gca;
% ax.LineWidth=0.5;%���������߿�Ϊ0.5��
% ax.XAxis.LineWidth=1;%����x�������߿�Ϊ1��
% ax.YAxis.LineWidth=1;%����y�������߿�Ϊ1��

for i=1:M
    delay(i)=finddelay(X(1,:),X(i,:))/Fs;%���������ʱ��λ��s
    [a,b]=xcorr(X(1,:),X(i,:),'coeff');
    [~,v]=max(abs(a));
    SampleDiff=b(v);
    %delayT=(b+(N-1))/Fs;
    %delayxcorr(i)=delayT(v);
    delayxcorr(i)=SampleDiff/Fs;
end



% figure
% plot(delayxcorr)
% hold on
% plot(delay)


XafFilter=m_BandPassFilter(X);%�˲�

%% ������ʱ����
% DESCRIPTIVE TEXT

r2central=sqrt(sum((SpeakerPosition-CentralPosition).^2));%��Դ���������ľ���
for i=1:M
    tau(i)=sqrt((SpeakerPosition(1)-Array_X(i))^2+(SpeakerPosition(2)-Array_Y(i))^2+(SpeakerPosition(3)-Array_Z(i))^2)/C-...
        sqrt((SpeakerPosition(1)-Central_X)^2+(SpeakerPosition(2)-Central_Y)^2+(SpeakerPosition(3)-Central_Z)^2)/C;
end

disp(['��Դ���������ľ���Ϊ: ',num2str(r2central),'m'])
if r2central<=2*L^2/lambda
    disp('����')
else
    disp('Զ��')
end
disp(['��������X: ',num2str(SpeakerPosition(1)),'��������Y: ',num2str(SpeakerPosition(2))])

% figure
% error=delay-tau;
% m=linspace(1,M,M);
% plot(m,error,'-o')
% y_val=get(gca,'YTick');   %Ϊ�˻��y����
% y_str=num2str(y_val');    %Ϊ�˽�����ת��Ϊ�ַ�����
% set(gca,'YTickLabel',y_str);    %��ʾ
% xlim([1 M])
% xlabel('��Ԫ���','FontName','����');ylabel('��ʱ���','FontName','����');
% title('������ʱ��ʵ����ʱ���','FontName','����')
% figure
% error=delayxcorr-tau;
% m=linspace(1,M,M);
% plot(m,error,'-o')
% y_val=get(gca,'YTick');   %Ϊ�˻��y����
% y_str=num2str(y_val');    %Ϊ�˽�����ת��Ϊ�ַ�����
% set(gca,'YTickLabel',y_str);    %��ʾ
% xlim([1 M])
% xlabel('��Ԫ���','FontName','����');ylabel('��ʱ���','FontName','����');
% title('finddelay��xcorr���','FontName','����')

%% ��Դ��λ�㷨
% DESCRIPTIVE TEXT

P=DOAalogrithm(XafFilter,Alogrithm);%��λ�㷨


[Y_max,X_max]=find(P==max(P(:)));
disp(['Ԥ������X: ',num2str(gridX(X_max)),'Ԥ������Y: ',num2str(gridY(Y_max))])
for m = 1:length(gridY)
    pp(m) = max(P(m,:)); % Pcbf �ĵ�k1�е����Ԫ�ص�ֵ
end

P = P/max(pp);  % ����Ԫ�س��������ֵ ��һ������
% P=20*log10(P/max(P(:)));
% ii=find(P<=max(P(:))-10);
% P(ii)=NaN;%������ʾ��Χ
figure
h=pcolor(gridX,gridY,P);%�ȸ���ͼ
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
xlabel('x/(m)','FontName','����');ylabel('y/(m)','FontName','����')
set(h,'edgecolor','none','facecolor','interp');%ȥ������ƽ������
colorbar %���ɫ��
map=jet(256);%��չ��256ɫ����ɫ���޷ֲ�
colormap(map);
hold on
plot([SpeakerPosition(1) SpeakerPosition(1)],[SpeakerPosition(2) SpeakerPosition(2)-1.5],'k','linewidth',2);%��˷�֧��
plot([-1.5 1.5],[-1.287 -1.287],'k','linewidth',2);%��˷�֧��
r = 0.1;%�뾶
a = SpeakerPosition(1);%Բ�ĺ�����
b = SpeakerPosition(2);%Բ��������
para = [a-r, b-r, 2*r, 2*r];
rectangle('Position', para, 'Curvature', [1 1],'EdgeColor', 'k','linewidth',2);
% rectangle('Position',[-1.5 -1.5 1 1],'edgecolor','k','linewidth',1.8) 
plot(gridX(X_max),gridY(Y_max),'g*')
text(gridX(X_max),gridY(Y_max),['  X=',num2str(gridX(X_max)),newline,'  Y=',num2str(gridY(Y_max)),newline],'Color','g');
scatter(SpeakerPosition(1),SpeakerPosition(2),20,'kp')%��Դλ��
titlestring=strcat('����Դ',num2str(FMAIN),'Hz ','��λ�㷨: ',alostring);
title(titlestring,'fontname','����','FontSize',12);
axis equal
ax=gca;
ax.LineWidth=0.5;%���������߿�Ϊ0.5��
ax.XAxis.LineWidth=1;%����x�������߿�Ϊ1��
ax.YAxis.LineWidth=1;%����y�������߿�Ϊ1��

figure
h=pcolor(gridX,gridY,P);%�ȸ���ͼ
% xlim([-2,2]);
% ylim([-2,2]);
% xlabel('x/(m)','FontName','����');ylabel('y/(m)','FontName','����')
set(h,'edgecolor','none','facecolor','interp');%ȥ������ƽ������
% colorbar %���ɫ��
map=jet(256);%��չ��256ɫ����ɫ���޷ֲ�
colormap(map);
% axis equal
set(gca,'position',[0 0 1 1])
set(gca,'XTick',[],'YTick',[]);
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca, 'Position', get(gca, 'OuterPosition') - ...
%     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
saveas(h,'D:\Programming\Matlab\MicrophoneArray\picture\soundfield.jpg');
audio=imread('D:\Programming\Matlab\MicrophoneArray\picture\soundfield.jpg');%��ȡ����ͼ��
photo=imread('D:\Programming\Matlab\MicrophoneArray\picture\r2_speaker.jpg');%��ȡ��Ƶͼ��
rect = [480 160 1600 1600];
photo=imcrop(photo,rect);   
% figure
% imshow(photo);    
audio=imresize(audio,[800,800]);
photo=imresize(photo,[800,800]);
% siz=size(audio);
% alpha=ones(siz(1),siz(2));
% alpha(B==0)=0;%������ɫ����Ϊ͸��
% imwrite(rgb,'test.png','Alpha',alpha)
addPhoto = imadd(0.4*audio,0.6*photo);
figure
imshow(addPhoto);%��ʾͼ��

% figure%��άͼ
% mesh(gridX,gridY,P,'FaceColor','interp');%3D figure
% xlabel('x/(m)','FontName','����');ylabel('y/(m)','FontName','����');zlabel('��ֵ/dB','FontName','����');
% colorbar %���ɫ��
% map=jet(256);%��չ��256ɫ����ɫ���޷ֲ�
% colormap(map);
% titlestring=strcat('��������Դ',num2str(FMAIN),'Hz');
% title(titlestring,'fontname','����','FontSize',12);
% ax=gca;
% ax.LineWidth=0.5;%���������߿�Ϊ0.5��
% ax.XAxis.LineWidth=1;%����x�������߿�Ϊ1��
% ax.YAxis.LineWidth=1;%����y�������߿�Ϊ1��