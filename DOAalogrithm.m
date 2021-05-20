function P=DOAalogrithm(X,Alogrithm)
global N M FMAIN C Array_X Array_Y Array_Z Central_X Central_Y Central_Z
global gridX gridY gridZ
switch(Alogrithm)
    case 1
        disp('定位算法:beamforming')
        R=X*X'/N;
        for i=1:M%去除对角线影响
            R(i,i)=0;
        end
        P_CB=zeros(length(gridY),length(gridX));
        for i=1:length(gridY)
            for j=1:length(gridX)
                Ri = sqrt((gridX(j)-Array_X).^2+(gridY(i)-Array_Y).^2+(gridZ-Array_Z).^2);% 该扫描点到各阵元的聚焦距离矢量
                Ri2 = sqrt((gridX(j)-Central_X).^2+(gridY(i)-Central_Y).^2+(gridZ-Central_Z).^2);% 扫描点到参考阵元的程差矢量
                Rn = Ri-Ri2;
                A = exp(-1i*2*pi*FMAIN*Rn/C); % 声压聚焦方向矢量
                P_CB(i,j) = abs(A'*R*A); % CSM
            end
        end
        P=P_CB;
    case 2
        disp('定位算法:MUSIC')
        R=X*X'/N;
        [EV,D]=eig(R);% [V,D]=eig(A)：求矩阵A的全部特征值，构成对角阵D，并求A的特征向量构成V的列向量
        %EV 特征向量8*8=64个
        %D  特征值构成对角矩阵  M个
        [EVA,I]=sort(diag(D).');   %特征值按升序排列
        EV=fliplr(EV(:,I));        %左右翻转，特征值按降序排列
        Un=EV(:,2:end);          %噪声子空间
        P_MUSIC=zeros(length(gridY),length(gridX));
        for i=1:length(gridY)
            for j=1:length(gridX)
                Ri = sqrt((gridX(j)-Array_X).^2+(gridY(i)-Array_Y).^2+(gridZ-Array_Z).^2);% 该扫描点到各阵元的聚焦距离矢量
                Ri2 = sqrt((gridX(j)-Central_X).^2+(gridY(i)-Central_Y).^2+(gridZ-Central_Z).^2);% 扫描点到各阵元与参考阵元的程差矢量
                Rn = Ri-Ri2;
                A = exp(-1i*2*pi*FMAIN*Rn/C); % 声压聚焦方向矢量
                P_MUSIC(i,j) = abs((A'*A)/(A'*(Un*Un')*A)); % CSM
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
        P_FB=zeros(length(gridY),length(gridX));
        for i=1:length(gridY)
            for j=1:length(gridX)
                Ri = sqrt((gridX(j)-Array_X).^2+(gridY(i)-Array_Y).^2+(gridZ-Array_Z).^2);% 该扫描点到各阵元的聚焦距离矢量
                Ri2 = sqrt((gridX(j)-Central_X).^2+(gridY(i)-Central_Y).^2+(gridZ-Central_Z).^2);% 扫描点到参考阵元的程差矢量
                Rn = Ri-Ri2;
                A = exp(-1i*2*pi*FMAIN*Rn/C); % 声压聚焦方向矢量
                P_FB(i,j) = abs(A'*R*A)^v; % CSM
            end
        end
        P=P_FB;
    case 4
        xmin=-3;xmax=3;ymin=-3;ymax=3;%检测平面大小
        increment=0.05;
        gridZ=4;
        xsteps=(abs(xmax-xmin)+increment)/increment;
        ysteps=(abs(ymax-ymin)+increment)/increment;
        xsteps=int16(xsteps)
        ysteps=int16(ysteps)
        [U,V]=meshgrid(xmin:increment:xmax,ymin:increment:ymax);
        W=ones(xsteps,ysteps)*gridZ;%这个可以当作检测平面各点距阵列平面的距离
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
        
end
end