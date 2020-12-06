%%%%    DebyeWolfDiffractionCalculations_For Test
%%%%    unit: um

clear all;
clc;

global lamda k n1 NA fo R

lamda=118.8;                                                               
k=2*pi/lamda;                                                               
n1 = 1;                                                                
NA = 0.93;                                                                  
fo=100*lamda;                                                                   
R=250*lamda;                                                               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  input plane  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nin=1024;   % 用DebyeWolf时Nin=1024,用二维矢量角谱时Nin=3000
Dx=2*R/Nin;                                                      
[xx,yy]=meshgrid(-(Nin-1)/2:(Nin-1)/2,-(Nin-1)/2:(Nin-1)/2);
Y=xx*Dx;
X=yy*Dx;
Aperture=sign(1-sign(xx.^2+yy.^2-((Nin-1)./2).^2));
g=-mod(2*pi/lamda*(sqrt(X.^2+Y.^2+fo^2)-fo),2*pi)+2*pi;% 双曲线相位分布——实现聚焦功能！
E=Aperture.*exp(1i.*g);         

clear  Aperture g         

th=asin(NA.*sqrt(xx.^2 + yy.^2)./((Nin-1)/2*n1));                             
thh=th;
th(thh>asin(NA./n1))=0;                                                     
E(thh>asin(NA./n1))=0;                                                     
phi = atan2 (yy,xx);                                                         
phi(phi<0) = phi(phi<0)+2.*pi;
z = 0e0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  45Line-polarization illuminating  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  using DebyeWolf                     
Ex = E./sqrt(cos(th)).*(cos(pi/4).*(1+(cos(th)-1).*cos(phi).^2)+sin(pi/4).*(cos(th)-1).*cos(phi).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th));
Ey = E./sqrt(cos(th)).*(cos(pi/4).*(cos(th)-1).*cos(phi).*sin(phi)+sin(pi/4).*(1+(cos(th)-1).*sin(phi).^2)).*exp (i.*k.*n1.*z.*cos (th));
Ez = E./sqrt(cos(th)).*(cos(pi/4).*sin(th).*cos(phi)+sin(pi/4).*sin(th).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th));
%%%%%%%%%%%%%%%%%%%  using 2D-VASD
% Ex=E.*cos(pi/4);
% Ey=E.*sin(pi/4);

clear E th  phi xx yy
Ex (real (thh)>asin (NA./n1))=0;                                          
Ey (real (thh)>asin (NA./n1))=0;
Ez (real (thh)>asin (NA./n1))=0;
where=isnan (Ex);Ex (where)=1;                                                             
where=isnan (Ey);Ey (where)=1;                                                              
where=isnan (Ez);Ez (where)=1;        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Diffraction field calculation  %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%  using 2D Vectorial Angular Spectrum method 
% Ex0=Ex;
% Ey0=Ey;
% [Ex,Ey,Ez]=Diffraction2DTransPolar(Ex0,Ey0,z+fo,lamda,n1,Dx,Nin);
%%%%%%%%%%%%%%%%%%%  using FFT method  
Ex0=padarray(Ex,[2*Nin 2*Nin],0);
Ey0=padarray(Ey,[2*Nin 2*Nin],0);
Ez0=padarray(Ez,[2*Nin 2*Nin],0);
Ex=fftshift(fft2(Ex0));
Ey=fftshift(fft2(Ey0));
Ez=fftshift(fft2(Ez0));
clear Ex0 Ey0 Ez0 thh where
[r,c]=size(Ex);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  xy plane Display  %%%%%%%%%%%%%%%%%%%%%%%%                                                 
Ix=abs(Ex).^2;
Iy=abs(Ey).^2;
Iz = abs(Ez).^2;                                                        
Itotal=Iz+Ix+Iy;
clear Ix Iy  Ex Ey Ez

%%%%%%%%%%%%%%%%%%%  for DebyeWolf method  
L = fo.*lamda./n1.*Nin./(2.*R);  
Dxx=L/r;
Ndisplay=round(10*lamda/Dxx);   %   截取所需展示的范围
x = linspace(-L/2,L/2,r);
y = x;
x=x(round((r-Ndisplay-1)/2):round((r+Ndisplay-1)/2));
y=x;
IzDisplay=Iz(round((r-Ndisplay-1)/2):round((r+Ndisplay-1)/2),round((r-Ndisplay-1)/2):round((r+Ndisplay-1)/2));
ItotalDisplay=Itotal(round((r-Ndisplay-1)/2):round((r+Ndisplay-1)/2),round((r-Ndisplay-1)/2):round((r+Ndisplay-1)/2));

figure
surf(x/lamda,y/lamda,IzDisplay); colormap(hot); shading interp  % DebyeWolf:x y IzDisplay; 2D-VASD:X Y Iz
colorbar
axis equal
set(0,'defaultfigurecolor','w') %设置图片的背景颜色为白色
set(gcf,'unit','normalized','position',[0.25,0.1,0.5,0.76]);%设置图片窗口的大小
set(gca,'FontName','Times New Roman','FontWeight','bold','FontSize',20);
view(2)    
xlabel('X(\lambda)');ylabel('Y(\lambda)');zlabel('Intensity (a.u.)');
title('Iz(45Line polarization-2D-DebyeWolf-FFT)');
set(gca,'XTick',-5:1:5,'XLim',[-5 5]);%X轴的数据显示范围%这句话对显示3D图很漂亮
set(gca,'YTick',-5:1:5,'YLim',[-5 5]);%X轴的数据显示范围


