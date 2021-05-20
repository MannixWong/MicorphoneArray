function XafFilter=m_BandPassFilter(X)
global N Fs C M FMAIN lambda
ff=(0:N/2-1)*Fs/N;
for i=1:M
    X(i,:)=X(i,:)-mean(X(i,:));    %去除直流分量
    XFbfFilter(i,:)=fft(X(i,:),N);
    FbfFilter(i,:)=abs(XFbfFilter(i,1:N/2));
end


FbfFilter=FbfFilter*2/N;
FbfFilter(1,:)=FbfFilter(1,:)/2;
for i=1:M
    [~,maxFindex]=max(FbfFilter(i,:));
    FMAIN(i)=ff(maxFindex);%信号主频
end
% FMAIN=sum(FMAIN)/M;
FMAIN=3886;
lambda=C/FMAIN;%波长
disp(['信号的主频为: ',num2str(FMAIN),'Hz'])
wp=[FMAIN-50 FMAIN+50].*2./Fs;%通带边界
ws=[FMAIN-100 FMAIN+100].*2./Fs;%阻带边界
rp=1;rs=40;
[n,wn]=ellipord(wp,ws,rp,rs);
[b,a]=ellip(n,rp,rs,wn);
for i=1:M
    XafFilter(i,:)=filter(b,a,X(i,:));%滤波
    XFafFilter(i,:)=fft(XafFilter(i,:),N);
    FafFilter(i,:)=abs(XFafFilter(i,1:N/2));
end


FafFilter=FafFilter*2/N;
FafFilter(1)=FafFilter(1)/2;

figure
subplot(2,1,1)
plot(ff,FbfFilter(2,:),'k');
xlabel('频率/Hz','FontName','宋体');ylabel('幅值','FontName','宋体');
xlim([0 10000]);
title('滤波前信号频率','FontName','黑体')

subplot(2,1,2)
plot(ff,FafFilter(2,:),'k');
xlabel('频率/Hz','FontName','宋体');ylabel('幅值','FontName','宋体');
xlim([0 10000]);
title('滤波后信号频率','FontName','黑体')
end