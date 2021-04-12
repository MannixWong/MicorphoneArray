function P=DOAalogrithm(X,Alogrithm)
global N M FMAIN C Array_X Array_Y
switch(Alogrithm)
    case 1
        disp('��λ�㷨:beamforming')
        R=X*X'/N;
        for i=1:M%ȥ���Խ���Ӱ��
            R(i,i)=0;
        end
        A=zeros(M,1);
        thetac=zeros(1,180);
        phic=zeros(1,90);
        P_CB=zeros(180,90);
        for azi=1:1:180
            for ele=1:1:90
                thetac(azi)=azi-1;
                phic(ele)=ele-1;
                for m=1:M
                    A(m,1)=exp(-1i*2*pi*FMAIN/C*(Array_X(m)*cos(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)+Array_Y(m)*sin(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)));%�����󻡶Ȳ���
                    %A(m,1)=exp(-1i*2*pi*1/lambda*(Array_X(m)*cosd(thetac(azi))*cosd(phic(ele))+Array_Y(m)*sind(thetac(azi))*cosd(phic(ele))));%�����󻡶Ȳ���
                end
                P_CB(azi,ele)=A'*R*A;%the power of beamforming
            end
        end
        P_CB=abs(P_CB);
        P=P_CB;
    case 2
        disp('��λ�㷨:MUSIC')
        R=X*X'/N;
        [EV,D]=eig(R);% [V,D]=eig(A)�������A��ȫ������ֵ�����ɶԽ���D������A��������������V��������
        %EV ��������8*8=64��
        %D  ����ֵ���ɶԽǾ���  M��
        [EV,I]=sort(diag(D).');   %����ֵ����������
        EV=fliplr(EV(:,I));        %���ҷ�ת������ֵ����������
        Un=EV(:,2:end);          %�����ӿռ�
        A=zeros(M,1);
        thetac=zeros(1,180);
        phic=zeros(1,90);
        P_MUSIC=zeros(180,90);
        for azi=1:1:180
            for ele=1:1:90
                thetac(azi)=azi-1;
                phic(ele)=ele-1;
                for m=1:M
                    A(m,1)=exp(-1i*2*pi*FMAIN/C*(Array_X(m)*cos(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)+Array_Y(m)*sin(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)));%�����󻡶Ȳ���
                    %             A(m,1)=exp(-1i*2*pi*1/lambda*(Array_X(m)*cosd(thetac(azi))*cosd(phic(ele))+Array_Y(m)*sind(thetac(azi))*cosd(phic(ele))));%�����󻡶Ȳ���
                end
                P_MUSIC(azi,ele)=A'*R*A;%the power of beamforming
            end
        end
        P_MUSIC=abs(P_MUSIC);
        P=P_MUSIC;
    case 3
        disp('��λ�㷨:Functional beamforming')
        v=input('������v����ֵ: ');
        R=X*X'/N;%the covariance martix of received signals
        R=R^(1/v);
%         % for i=1:M%ȥ���Խ���Ӱ��
%         %     R(i,i)=0;
%         % end
        A=zeros(M,1);
        thetac=zeros(1,180);
        phic=zeros(1,90);
        P_FB=zeros(180,90);
        for azi=1:1:180
            for ele=1:1:90
                thetac(azi)=azi-1;
                phic(ele)=ele-1;
                for m=1:M
                    A(m,1)=exp(-1i*2*pi*FMAIN/C*(Array_X(m)*cos(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)+Array_Y(m)*sin(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)));%�����󻡶Ȳ���
                end
                P_FB(azi,ele)=(A'*R*A)^v;%the power of beamforming
                
            end
        end
        P_FB=abs(P_FB);
        P=P_FB;
        
end
end