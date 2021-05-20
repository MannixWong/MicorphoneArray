function P=DOAalogrithm(X,Alogrithm)
global N M FMAIN C Array_X Array_Y Array_Z Central_X Central_Y Central_Z
global gridX gridY gridZ
switch(Alogrithm)
    case 1
        disp('��λ�㷨:beamforming')
        R=X*X'/N;
        for i=1:M%ȥ���Խ���Ӱ��
            R(i,i)=0;
        end
        P_CB=zeros(length(gridY),length(gridX));
        for i=1:length(gridY)
            for j=1:length(gridX)
                Ri = sqrt((gridX(j)-Array_X).^2+(gridY(i)-Array_Y).^2+(gridZ-Array_Z).^2);% ��ɨ��㵽����Ԫ�ľ۽�����ʸ��
                Ri2 = sqrt((gridX(j)-Central_X).^2+(gridY(i)-Central_Y).^2+(gridZ-Central_Z).^2);% ɨ��㵽�ο���Ԫ�ĳ̲�ʸ��
                Rn = Ri-Ri2;
                A = exp(-1i*2*pi*FMAIN*Rn/C); % ��ѹ�۽�����ʸ��
                P_CB(i,j) = abs(A'*R*A); % CSM
            end
        end
        P=P_CB;
    case 2
        disp('��λ�㷨:MUSIC')
        R=X*X'/N;
        [EV,D]=eig(R);% [V,D]=eig(A)�������A��ȫ������ֵ�����ɶԽ���D������A��������������V��������
        %EV ��������8*8=64��
        %D  ����ֵ���ɶԽǾ���  M��
        [EVA,I]=sort(diag(D).');   %����ֵ����������
        EV=fliplr(EV(:,I));        %���ҷ�ת������ֵ����������
        Un=EV(:,2:end);          %�����ӿռ�
        P_MUSIC=zeros(length(gridY),length(gridX));
        for i=1:length(gridY)
            for j=1:length(gridX)
                Ri = sqrt((gridX(j)-Array_X).^2+(gridY(i)-Array_Y).^2+(gridZ-Array_Z).^2);% ��ɨ��㵽����Ԫ�ľ۽�����ʸ��
                Ri2 = sqrt((gridX(j)-Central_X).^2+(gridY(i)-Central_Y).^2+(gridZ-Central_Z).^2);% ɨ��㵽����Ԫ��ο���Ԫ�ĳ̲�ʸ��
                Rn = Ri-Ri2;
                A = exp(-1i*2*pi*FMAIN*Rn/C); % ��ѹ�۽�����ʸ��
                P_MUSIC(i,j) = abs((A'*A)/(A'*(Un*Un')*A)); % CSM
            end
        end
        P_MUSIC=abs(P_MUSIC);
        P=P_MUSIC;
    case 3
        disp('��λ�㷨:Functional beamforming')
        v=input('������v����ֵ: ');
        R=X*X'/N;%the covariance martix of received signals
        R=R^(1/v);
        % for i=1:M%ȥ���Խ���Ӱ��
        %     R(i,i)=0;
        % end
        P_FB=zeros(length(gridY),length(gridX));
        for i=1:length(gridY)
            for j=1:length(gridX)
                Ri = sqrt((gridX(j)-Array_X).^2+(gridY(i)-Array_Y).^2+(gridZ-Array_Z).^2);% ��ɨ��㵽����Ԫ�ľ۽�����ʸ��
                Ri2 = sqrt((gridX(j)-Central_X).^2+(gridY(i)-Central_Y).^2+(gridZ-Central_Z).^2);% ɨ��㵽�ο���Ԫ�ĳ̲�ʸ��
                Rn = Ri-Ri2;
                A = exp(-1i*2*pi*FMAIN*Rn/C); % ��ѹ�۽�����ʸ��
                P_FB(i,j) = abs(A'*R*A)^v; % CSM
            end
        end
        P=P_FB;
    case 4
        xmin=-3;xmax=3;ymin=-3;ymax=3;%���ƽ���С
        increment=0.05;
        gridZ=4;
        xsteps=(abs(xmax-xmin)+increment)/increment;
        ysteps=(abs(ymax-ymin)+increment)/increment;
        xsteps=int16(xsteps)
        ysteps=int16(ysteps)
        [U,V]=meshgrid(xmin:increment:xmax,ymin:increment:ymax);
        W=ones(xsteps,ysteps)*gridZ;%������Ե������ƽ����������ƽ��ľ���
        Scanrc=sqrt(U.^2+V.^2+W.^2);
        % mesh(U,V,Scanrc);%����Ǹ����������ĵ�ľ���
        eletest_r=asin(W./Scanrc);
        azitest_r=atan2(U,V);
        
        eletest_d=asind(W./Scanrc);
        % azitest_d=atand(U./V);
        azitest_d=atan2d(U,V);
        
        uu=find(azitest_d(:)<0);
        azitest_d(uu)=360+azitest_d(uu);
        
        R=X*X'/N;%the covariance martix of received signals
        for i=1:M%ȥ���Խ���Ӱ��
            R(i,i)=0;
        end
        
        for i=1:xsteps
            for j=1:ysteps
                for m=1:M
                    A(m,1)=exp(-1i*2*pi*FMAIN/C*(Array_X(m).*cos(azitest_r(i,j))*cos(eletest_r(i,j))+Array_Y(m).*sin(azitest_r(i,j))*cos(eletest_r(i,j))));
                end
                P_SP(i,j)=A'*R*A;
            end
        end
        P_SP=abs(P_SP);
        P=P_SP;
        
end
end