function XafFilter=m_BandPassFilter(X)
global N Fs C M FMAIN lambda
ff=(0:N/2-1)*Fs/N;
for i=1:M
    X(i,:)=X(i,:)-mean(X(i,:));    %ȥ��ֱ������
    XFbfFilter(i,:)=fft(X(i,:),N);
    FbfFilter(i,:)=abs(XFbfFilter(i,1:N/2));
end


FbfFilter=FbfFilter*2/N;
FbfFilter(1,:)=FbfFilter(1,:)/2;
for i=1:M
    [~,maxFindex]=max(FbfFilter(i,:));
    FMAIN(i)=ff(maxFindex);%�ź���Ƶ
end
% FMAIN=sum(FMAIN)/M;
FMAIN=3886;
lambda=C/FMAIN;%����
disp(['�źŵ���ƵΪ: ',num2str(FMAIN),'Hz'])
wp=[FMAIN-50 FMAIN+50].*2./Fs;%ͨ���߽�
ws=[FMAIN-100 FMAIN+100].*2./Fs;%����߽�
rp=1;rs=40;
[n,wn]=ellipord(wp,ws,rp,rs);
[b,a]=ellip(n,rp,rs,wn);
for i=1:M
    XafFilter(i,:)=filter(b,a,X(i,:));%�˲�
    XFafFilter(i,:)=fft(XafFilter(i,:),N);
    FafFilter(i,:)=abs(XFafFilter(i,1:N/2));
end


FafFilter=FafFilter*2/N;
FafFilter(1)=FafFilter(1)/2;

figure
subplot(2,1,1)
plot(ff,FbfFilter(2,:),'k');
xlabel('Ƶ��/Hz','FontName','����');ylabel('��ֵ','FontName','����');
xlim([0 10000]);
title('�˲�ǰ�ź�Ƶ��','FontName','����')

subplot(2,1,2)
plot(ff,FafFilter(2,:),'k');
xlabel('Ƶ��/Hz','FontName','����');ylabel('��ֵ','FontName','����');
xlim([0 10000]);
title('�˲����ź�Ƶ��','FontName','����')
end