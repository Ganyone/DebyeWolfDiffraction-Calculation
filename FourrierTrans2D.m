function [G]=FourrierTrans2D(g, Dx, N, flag)
   num=0:N-1;
   a=exp(i*2*pi/N*(N/2-0.5)*num);
   A=a'*a;
   C=exp(-i*2*pi/N*(N/2-0.5).^2*2)*A ;
   if flag==1
       G=fft2(A.*g); 
       G=[G(2:N,:); G(1,:)];
       G=Dx.^2.*C.*G; 
   end
   if flag==-1
     G=ifft2(conj(A).*g);
     G= (1./(N*Dx)).^2.*N^2*conj(C).*G;
   end
   
end