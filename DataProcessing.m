
clear;close all;clc
M=36
filen='C:\Users\MannixWong\Desktop\data\20210412_16_57_12.dat';
height=1000000;
fid=fopen(filen,'r');%����
img=fread(fid,[M,height],'uint16');%�����ݶ���1*40000�������У�ָ��Դ����Ϊuint16��
img1=mod(img,256);%ȡģ����(���������)������imgÿ��Ԫ�س���256������
img2=fix(img/256);%������������
Fs=44100;
img=img1*256+img2;
[m,n]=size(img);%1*40000
for i =1:m
    for j=1:n
        if img(i,j)>=32768
            img(i,j)=-(img(i,j)-32768)*5/32768;
        else
            img(i,j) = img(i,j)*5/32768;
        end
    end
end
Startpoint=3000;
Sample=1000;
x=1:1:Sample;
for i=1:M
y(i,:)=img(i ,Startpoint:Startpoint+Sample-1);
y(i,:)=y(i,:)-(max(y(i,:))+min(y(i,:)))/2;
[y(i,:),PS]=mapminmax(y(i,:));
end

%% ����Է���
% DESCRIPTIVE TEXT
[a,b]=xcorr(y(1,:),y(16,:),'coeff');
figure
plot(b,a)


figure
plot(x/Fs,y(1,:));
hold on
plot(x/Fs,y(16,:));
% legend('ͨ��1','ͨ��8')
grid on
ax=gca;
ax.LineWidth=0.5;%���������߿�Ϊ0.5��
ax.XAxis.LineWidth=1;%����x�������߿�Ϊ1��
ax.YAxis.LineWidth=1;%����y�������߿�Ϊ1��
X=y;
for i=1:M-1
delay(i)=finddelay(X(1,:),X(i+1,:))
end
string=strcat('E',num2str(M));
save(string,'X')

figure('Name','1-8ͨ��')
for i=1:8
    subplot(4,2,i)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
      plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %Ϊ�˻��y����
%     y_str=num2str(y_val');    %Ϊ�˽�����ת��Ϊ�ַ�����
%     set(gca,'YTickLabel',y_str);    %��ʾ
    titlestring=strcat('��',num2str(i),'ͨ��');
    xlabel('������','FontSize',12);
    ylabel('��ֵ','FontSize',12);
    title(titlestring,'fontname','����','FontSize',12);
    grid on
end
figure('Name','9-16ͨ��')
for i=9:16
    subplot(4,2,i-8)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
    plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %Ϊ�˻��y����
%     y_str=num2str(y_val');    %Ϊ�˽�����ת��Ϊ�ַ�����
%     set(gca,'YTickLabel',y_str);    %��ʾ
    titlestring=strcat('��',num2str(i),'ͨ��');
    xlabel('������','FontSize',12);
    ylabel('��ֵ','FontSize',12);
    title(titlestring,'fontname','����','FontSize',12);
    grid on
end
figure('Name','17-24ͨ��')
for i=17:24
    subplot(4,2,i-16)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %Ϊ�˻��y����
%     y_str=num2str(y_val');    %Ϊ�˽�����ת��Ϊ�ַ�����
%     set(gca,'YTickLabel',y_str);    %��ʾ
    titlestring=strcat('��',num2str(i),'ͨ��');
    xlabel('������','FontSize',12);
    ylabel('��ֵ','FontSize',12);
    title(titlestring,'fontname','����','FontSize',12);
    grid on
end
figure('Name','25-32ͨ��')
for i=25:32
    subplot(4,2,i-24)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
    plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %Ϊ�˻��y����
%     y_str=num2str(y_val');    %Ϊ�˽�����ת��Ϊ�ַ�����
%     set(gca,'YTickLabel',y_str);    %��ʾ
    titlestring=strcat('��',num2str(i),'ͨ��');
    xlabel('������','FontSize',12);
    ylabel('��ֵ','FontSize',12);
    title(titlestring,'fontname','����','FontSize',12);
    grid on
end
figure('Name','33-36ͨ��')
for i=33:36
    subplot(2,2,i-32)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %Ϊ�˻��y����
%     y_str=num2str(y_val');    %Ϊ�˽�����ת��Ϊ�ַ�����
%     set(gca,'YTickLabel',y_str);    %��ʾ
    titlestring=strcat('��',num2str(i),'ͨ��');
    xlabel('������','FontSize',12);
    ylabel('��ֵ','FontSize',12);
    title(titlestring,'fontname','����','FontSize',12);
    grid on
end