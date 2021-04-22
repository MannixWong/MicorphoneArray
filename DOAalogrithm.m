function P=DOAalogrithm(X,Alogrithm)
global N M FMAIN C Array_X Array_Y
switch(Alogrithm)
    case 1
        disp('定位算法:beamforming')
        R=X*X'/N;
        for i=1:M%去除对角线影响
            R(i,i)=0;
        end
        A=zeros(M,1);
        thetac=zeros(1,360);
        phic=zeros(1,90);
        P_CB=zeros(360,90);
        for azi=1:1:360
            for ele=1:1:90
                thetac(azi)=azi-1;
                phic(ele)=ele-1;
                for m=1:M
                    A(m,1)=exp(-1i*2*pi*FMAIN/C*(Array_X(m)*cos(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)+Array_Y(m)*sin(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)));%螺旋阵弧度测试
%                     A(m,1)=exp(-1i*2*pi*FMAIN/C*(Array_X(m)*cosd(thetac(azi))*cosd(phic(ele))+Array_Y(m)*sind(thetac(azi))*cosd(phic(ele))));%螺旋阵弧度测试
                end
                P_CB(azi,ele)=A'*R*A;%the power of beamforming
            end
        end
        P_CB=abs(P_CB);
        
        P=P_CB;
    case 2
        disp('定位算法:MUSIC')
        R=X*X'/N;
        [EV,D]=eig(R);% [V,D]=eig(A)：求矩阵A的全部特征值，构成对角阵D，并求A的特征向量构成V的列向量
        %EV 特征向量8*8=64个
        %D  特征值构成对角矩阵  M个
        [EV,I]=sort(diag(D).');   %特征值按升序排列
        EV=fliplr(EV(:,I));        %左右翻转，特征值按降序排列
        Un=EV(:,2:end);          %噪声子空间
        A=zeros(M,1);
        thetac=zeros(1,360);
        phic=zeros(1,90);
        P_MUSIC=zeros(360,90);
        for azi=1:1:360
            for ele=1:1:90
                thetac(azi)=azi-1;
                phic(ele)=ele-1;
                for m=1:M
                    A(m,1)=exp(-1i*2*pi*FMAIN/C*(Array_X(m)*cos(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)+Array_Y(m)*sin(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)));%螺旋阵弧度测试
                    %A(m,1)=exp(-1i*2*pi*1/lambda*(Array_X(m)*cosd(thetac(azi))*cosd(phic(ele))+Array_Y(m)*sind(thetac(azi))*cosd(phic(ele))));%螺旋阵弧度测试
                end
                P_MUSIC(azi,ele)=A'*R*A;%the power of beamforming
            end
        end
        P_MUSIC=abs(P_MUSIC);
        P=P_MUSIC;
    case 3
        disp('定位算法:Functional beamforming')
        v=input('请输入v的数值: ');
        R=X*X'/N;%the covariance martix of received signals
        R=R^(1/v);
        % for i=1:M%去除对角线影响
        %     R(i,i)=0;
        % end
        A=zeros(M,1);
        thetac=zeros(1,360);
        phic=zeros(1,90);
        P_FB=zeros(360,90);
        for azi=1:1:360
            for ele=1:1:90
                thetac(azi)=azi-1;
                phic(ele)=ele-1;
                for m=1:M
                    A(m,1)=exp(-1i*2*pi*FMAIN/C*(Array_X(m)*cos(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)+Array_Y(m)*sin(thetac(azi)*pi/180)*cos(phic(ele)*pi/180)));%螺旋阵弧度测试
                end
                P_FB(azi,ele)=(A'*R*A)^v;%the power of beamforming
                
            end
        end
        P_FB=abs(P_FB);
        P=P_FB;
    case 4
        xmin=-3;xmax=3;ymin=-3;ymax=3;%检测平面大小
        increment=0.05;
        z=4;
        xsteps=(abs(xmax-xmin)+increment)/increment;
        ysteps=(abs(ymax-ymin)+increment)/increment;
        xsteps=int16(xsteps)
        ysteps=int16(ysteps)
        [U,V]=meshgrid(xmin:increment:xmax,ymin:increment:ymax);
        W=ones(xsteps,ysteps)*z;%这个可以当作检测平面各点距阵列平面的距离
        Scanrc=sqrt(U.^2+V.^2+W.^2);
        % mesh(U,V,Scanrc);%这个是各坐标点距中心点的距离
        eletest_r=asin(W./Scanrc);
        azitest_r=atan2(U,V);
        
        eletest_d=asind(W./Scanrc);
        % azitest_d=atand(U./V);
        azitest_d=atan2d(U,V);
        
        uu=find(azitest_d(:)<0);
        azitest_d(uu)=360+azitest_d(uu);
        
        R=X*X'/N;%the covariance martix of received signals
        for i=1:M%去除对角线影响
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
        [x_max,y_max]=find(P==max(P(:)));
%         disp(['预测方位角: ',num2str(azitest_d(x_max,y_max)),'°','预测俯仰角: ',num2str(eletest_d(x_max,y_max)),'°'])
end
end