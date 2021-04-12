
clear;close all;clc
M=36
filen='C:\Users\MannixWong\Desktop\data\20210412_16_57_12.dat';
height=1000000;
fid=fopen(filen,'r');%读出
img=fread(fid,[M,height],'uint16');%将数据读到1*40000的数组中，指定源数据为uint16类
img1=mod(img,256);%取模运算(除后的余数)，计算img每个元素除以256的余数
img2=fix(img/256);%朝零四舍五入
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

%% 相关性分析
% DESCRIPTIVE TEXT
[a,b]=xcorr(y(1,:),y(16,:),'coeff');
figure
plot(b,a)


figure
plot(x/Fs,y(1,:));
hold on
plot(x/Fs,y(16,:));
% legend('通道1','通道8')
grid on
ax=gca;
ax.LineWidth=0.5;%设置网格线宽为0.5磅
ax.XAxis.LineWidth=1;%设置x坐标轴线宽为1磅
ax.YAxis.LineWidth=1;%设置y坐标轴线宽为1磅
X=y;
for i=1:M-1
delay(i)=finddelay(X(1,:),X(i+1,:))
end
string=strcat('E',num2str(M));
save(string,'X')

figure('Name','1-8通道')
for i=1:8
    subplot(4,2,i)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
      plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %为了获得y轴句柄
%     y_str=num2str(y_val');    %为了将数字转换为字符数组
%     set(gca,'YTickLabel',y_str);    %显示
    titlestring=strcat('第',num2str(i),'通道');
    xlabel('采样点','FontSize',12);
    ylabel('幅值','FontSize',12);
    title(titlestring,'fontname','黑体','FontSize',12);
    grid on
end
figure('Name','9-16通道')
for i=9:16
    subplot(4,2,i-8)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
    plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %为了获得y轴句柄
%     y_str=num2str(y_val');    %为了将数字转换为字符数组
%     set(gca,'YTickLabel',y_str);    %显示
    titlestring=strcat('第',num2str(i),'通道');
    xlabel('采样点','FontSize',12);
    ylabel('幅值','FontSize',12);
    title(titlestring,'fontname','黑体','FontSize',12);
    grid on
end
figure('Name','17-24通道')
for i=17:24
    subplot(4,2,i-16)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %为了获得y轴句柄
%     y_str=num2str(y_val');    %为了将数字转换为字符数组
%     set(gca,'YTickLabel',y_str);    %显示
    titlestring=strcat('第',num2str(i),'通道');
    xlabel('采样点','FontSize',12);
    ylabel('幅值','FontSize',12);
    title(titlestring,'fontname','黑体','FontSize',12);
    grid on
end
figure('Name','25-32通道')
for i=25:32
    subplot(4,2,i-24)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
    plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %为了获得y轴句柄
%     y_str=num2str(y_val');    %为了将数字转换为字符数组
%     set(gca,'YTickLabel',y_str);    %显示
    titlestring=strcat('第',num2str(i),'通道');
    xlabel('采样点','FontSize',12);
    ylabel('幅值','FontSize',12);
    title(titlestring,'fontname','黑体','FontSize',12);
    grid on
end
figure('Name','33-36通道')
for i=33:36
    subplot(2,2,i-32)
%     plot(x,img(i,Startpoint:Startpoint+Sample-1))
plot(x,y(i,:))
%     y_val=get(gca,'YTick');   %为了获得y轴句柄
%     y_str=num2str(y_val');    %为了将数字转换为字符数组
%     set(gca,'YTickLabel',y_str);    %显示
    titlestring=strcat('第',num2str(i),'通道');
    xlabel('采样点','FontSize',12);
    ylabel('幅值','FontSize',12);
    title(titlestring,'fontname','黑体','FontSize',12);
    grid on
end