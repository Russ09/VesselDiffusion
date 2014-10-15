function [ VnessOld, VnessNew ] = anIsoDiffImage( Image )

addpath('frangi_filter');

Image = squeeze(Image);
sigmas = 1:1:5
options.FrangiScaleRange=[sigmas(1),sigmas(end)];
options.FrangiScaleRatio = 1;
options.BlackWhite = 0;

Image = padarray(Image,[10,10,10]);

if(ndims(Image)==3)
    
    iOld = double(Image);
    clear Image;
    
    Dk = 0.1;
    omega = 50;
    epsilon = 0.0001;
    dt = 0.02;
    s = 10;
    sigs = sigmas(1:3);
    
%     options.FrangiScaleRange=[sigs sigs];
%     options.FrangiScaleRatio = 2;
%     imagesc(iOld);
   
    Vness1 = FrangiFilter3D(iOld,options);
    
    for iT = 1:50
       
        [Vness, whatScale,V1x,V1y,V1z] = FrangiFilter3DRuss(iOld,options);
        
        V2x = zeros(size(V1x));
        V2y = zeros(size(V1x));
        V2z = zeros(size(V1x));
        
        V3x = zeros(size(V1x));
        V3y = zeros(size(V1x));
        V3z = zeros(size(V1x));
        
        xyZero = (V1x==0).*(V1y==0);
        
        V2x(V1z~=0) = -V1y(V1z~=0);
        V2y(V1z~=0) = V1x(V1z~=0);
        V2z(V1z==0) = 1;
        
        V2x(xyZero==1) = 1;
        V2y(xyZero==1) = 0;
        V2z(xyZero==1) = 0;
        
        V3x = (V1y.*V2z) - (V1z.*V2y);
        V3y = (V1z.*V2x) - (V2z.*V1x);
        V3z = (V1x.*V2y) - (V1y.*V2x);
        
        norm1 = sqrt(V1x.^2 + V1y.^2 + V1z.^2);
        norm2 = sqrt(V2x.^2 + V2y.^2 + V2z.^2);
        norm3 = sqrt(V3x.^2 + V3y.^2 + V3z.^2);
        
        V1x = V1x./norm1;
        V1y = V1y./norm1;
        V1z = V1z./norm1;
        
        V2x = V2x./norm2;
        V2y = V2y./norm2;
        V2z = V2z./norm2;
        
        V3x = V3x./norm3;
        V3y = V3y./norm3;
        V3z = V3z./norm3;

        
        
        
        m1 = 1 + (omega - 1)*Vness.^(1/s).*(iOld./max(iOld(:)));
        m2 = 1 + (epsilon - 1)*Vness.^(1/s);
        m3 =  m2;
        
        
        K11 = (m1.*V1x.^2 + m2.*V2x.*V2x + m3.*V3x.*V3x);
        K21 = (m1.*V1y.*V1x + m2.*V2y.*V2x  + m3.*V3y.*V3x);
        K31 = (m1.*V1z.*V1x + m2.*V2z.*V2x  + m3.*V3z.*V3x);
        
        K22 = (m1.*V1y.*V1y + m2.*V2y.*V2y  + m3.*V3y.*V3y);
        K32 = (m1.*V1z.*V1y + m2.*V2z.*V2y  + m3.*V3z.*V3y);
        
        K33 = (m1.*V1z.*V1z + m2.*V2z.*V2z  + m3.*V3z.*V3z);
        
        K12 = K21;
        K23 = K32;
        K13 = K31;
        
        
        [gX,gY,gZ] = gradient(iOld);
        Step = dt.*divergence(K11.*gX + K12.*gY + K13.*gZ,K21.*gX + K22.*gY + K23.*gZ,K31.*gX + K32.*gY + K33.*gZ);
        iNew = iOld + Step;
        
        iOld = iNew;
        
        figure(1)
        subplot(221)
        imagesc(squeeze(Vness(:,:,50)));
        set(gca,'YDir','normal')
        subplot(222)
        imagesc(squeeze(Vness1(:,:,50)));
        set(gca,'YDir','normal')
        subplot(223);
        imagesc(squeeze((Vness(:,:,50)>0.1).*V1z(:,:,50)./norm1(:,:,50)));
        set(gca,'YDir','normal')
        subplot(224);
        imagesc(squeeze(Image(:,:,50)));
        set(gca,'YDir','normal')
        
        figure(2)
        
        subplot(121)
        cla
        h = patch(isosurface(Vness1));
        axis equal
%         axis([45,75,0,50,0,75]);
        lighting gouraud
        camlight
        view([50,15]);
        subplot(122)
        cla
        h2 = patch(isosurface(Vness));
        axis equal
%         axis([45,75,0,50,0,75]);
        lighting gouraud
        camlight
        view([50,15]);
        name = 'VesselDiffuse';
        
        set(h2,'EdgeColor','none','FaceColor',[1,0.1,0.2]);
        set(h,'EdgeColor','none','FaceColor',[1,0.1,0.2]);
        
%          if iT == 1
%                   writerObj = VideoWriter([name,'.avi']);
%                   open(writerObj);
%                   frame = getframe(2);
%                   writeVideo(writerObj,frame);
%         else
%                   frame = getframe(2);
%                   writeVideo(writerObj,frame);
%         end

    end
end
% 
% figure(1)
% subplot(121)
% isosurface(Vness1);
% axis manual
% axis([0,100,0,100,0,100]);
% 
% lighting gouraud
% camlight
% subplot(122)
% h = patch(isosurface(Vness));
% set(h,'axis',[0,100,0,100,0,100]);
% axis equal
% axis();
% lighting gouraud
% camlight

if(ndims(Image)==2)

figure(3)
Vness1 = FrangiFilter2DRuss(Image,options);

for iS = 1:length(sigmas)

Vness = FrangiFilter2DRuss(Image,options);

iOld = Image;
    s = 3;
    
    iOld = Image;
    Dk = 0.01;
    omega = 10;
    epsilon = 0.000001;
    sigs = sigmas(iS);
    options.FrangiScaleRange=[sigs sigs];
    options.FrangiScaleRatio = 2;
%     imagesc(iOld);
    T = 0.5
    dt = 0.001;
    Time = 0:dt:T
    
    for iT = 1:length(Time)
        
        [gX gY] = gradient(iOld);
        K = 1./(1+(gX.^2 + gY.^2)./Dk^2);
        
        
        Vness = FrangiFilter2DRuss(iOld,options);
        
        [Dxx,Dxy,Dyy] = Hessian2D(iOld,sigs);
        [Lambda2,Lambda1,Ix,Iy,I2x,I2y]=eig2imageRuss(Dxx,Dxy,Dyy);
        m1 = 1 + (omega - 1)*Vness.^(1/s).*(iOld./max(iOld(:)));
        m2 = 1 + (epsilon - 1)*Vness.^(1/s);
        
%         V = [permute(Ix,[1,3,2]),permute(I2x,[1,3,2]);permute(Iy,[1,3,2]),permute(I2y,[1,3,2])]
%         Vals = [permute(Ix,[1,3,2]),permute(I2x,[1,3,2]);permute(Iy,[1,3,2]),permute(I2y,[1,3,2])]
%         I2x = zeros(size(Ix));
%         I2y = zeros(size(Iy));
        
        K11 = m1.*Ix.^2 + m2.*I2x.*Iy;
        K12 = m1.*Ix.*Iy + m2.*Iy.*I2y;
        K21 = m1.*I2x.*Ix + m2.*I2x.*I2y;
        K22 = m1.*I2x.*Iy + m2.*I2y.^2;
        
%         Norm1 = sqrt(K11.^2 + K21.^2);
%         Norm2 = sqrt(K12.^2 + K22.^2);

%         FroNorm = sqrt(K11.^2 + K21.^2 + K12.^2 + K22.^2);
        
%         K11(FroNorm~=0) = K11(FroNorm~=0)./FroNorm(FroNorm~=0);
%         K21(FroNorm~=0) = K21(FroNorm~=0)./FroNorm(FroNorm~=0);
%         K12(FroNorm~=0) = K12(FroNorm~=0)./FroNorm(FroNorm~=0);
%         K22(FroNorm~=0) = K22(FroNorm~=0)./FroNorm(FroNorm~=0);
        
        K11 = abs(K11);
        K22 = abs(K22);
        K12 = abs(K12);
        K21 = abs(K21);
        
%         figure(2)
%         subplot(221);
%         imagesc(K11);
%         subplot(222);
%         imagesc(K12);
%         subplot(223);
%         imagesc(K21);
%         subplot(224);
%         imagesc(K22);
%         
        iNew = iOld + dt*divergence(K11.*gX + K12.*gY,K21.*gX + K22.*gY);
        iNew = iNew.*(iNew>0);
        iOld = iNew;
        
        sum(iOld(:))
        pause(0.01);
    end
end




end

VnessOld = Vness1>0.015;
VnessNew = Vness>0.015;

