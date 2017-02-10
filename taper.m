function [yout,cs]=taper(yin,start,endN)
%%%%%%%%%%%%%%%%%%Taper function
% taper record for doing FFT
% SYNTAX: yout=taper(yin,start,end)

      N=length(yin);
      M1 = ceil(N*start);
      M2 = M1 + 1;
      x=[1:N];

      ANG = pi/M1;
      cs(1:M1)=(1.-cos(x(1:M1)*ANG))/2.;      

      M3 = ceil(N*endN);
      M5 = N-M3;
      M4 = M5 + 1;
      ANG =pi/M3;
      cs(M4:N)=fliplr((1.-cos(x(1:M3)*ANG))/2);
      cs(M1+1:M4-1)=1.;
       
      [S1, S2]=size(yin);
      if (S1>S2)
        cs=cs';
      end
      cs=cs.*length(cs)./sum(cs);

      yout=yin.*cs;
end

