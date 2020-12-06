function [Ex,Ey,Ez]=Diffraction2DTransPolar(Ex0,Ey0,Z,Wavelen0,n_refr,Dx,N)
     num=0:N-1;
     freq=1./(N*Dx)*(num-N/2+0.5);
     freq_x=freq'*ones(1,N);% ע�⣺���о����ˮƽ����ΪY����ֱ����ΪX
     freq_y=freq_x'; %freq_x  ע�⣺���о����ˮƽ����ΪY����ֱ����ΪX
     fz=sqrt((n_refr/Wavelen0).^2-freq_x.^2-freq_y.^2);
     [SpectrumX]=FourrierTrans2D(Ex0, Dx, N, 1);
     [SpectrumY]=FourrierTrans2D(Ey0, Dx, N, 1);
     SpectrumZ=-(freq_x.*SpectrumX+freq_y.*SpectrumY)./fz.*exp(i*2*pi*fz*Z);%freq_y  ע�⣺���о����ˮƽ����ΪY����ֱ����ΪX
     SpectrumX=SpectrumX.*exp(i*2*pi*fz*Z);
     SpectrumY=SpectrumY.*exp(i*2*pi*fz*Z);

     Ex=FourrierTrans2D(SpectrumX, Dx, N, -1);
     Ey=FourrierTrans2D(SpectrumY, Dx, N, -1);
     Ez=FourrierTrans2D(SpectrumZ, Dx, N, -1);  
        
end