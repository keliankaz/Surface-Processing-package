function [fx,PowerStructx]=frequency_spectrum(Z,dx,decimation)
% Runs a fast fourier transform for an N by M array (Z) along all profiles
% along x (along increasing M dimension).

% inputs: 
% Z: N by M detrended array aligned with increasing M (x direction) along slip 
% dx: point spacing (in milimeters)

% iunclude fy and powerstructy for running in both direction

% Decimate the data if necessary
[Nx,Ny]=size(Z); % size of cut array
dec=decimation;  % decimation factor
Z=Z(:,1:dec:Ny); % FOR returning fxi
Z=Z(1:dec:Nx,:); % FOR returning fyi
dx=dx*dec;
dxold=dx;
[Nx,Ny]=size(Z);



%%% take average power spectra over profiles in each direction
% results in average vectors in x and y directions
% for interpolation  

fxi=(0:Ny-1)'/(dx*Ny); 
fxi=fxi(1:floor(end/2));

warning off;  % don't display annoying interpolation warnings
% define noise level = 4 std curvature
%Z2=(del2(Z)/dx.^2)./2;

Z2=del2(Z,dx,dxold)./(dxold.*dx)./2;
Z2=Z2(isfinite(Z2));
NoiseLevel=8*std(Z2);
%NoiseLevel = 10^100;

% set minimum length of row or column for which power will be estimated
Nmin=10;

% initialize structures
maxNumSegment = floor(Ny/Nmin);

% I think this initializtion actually makes things slower 
% PowerStructx(Nx,maxNumSegment).px = [];
% PowerStructx(Nx,maxNumSegment).fx = [];

for i=1:Nx
    
    z   = Z(i,:)';
    Iz  = find(isfinite(z)); 
    
    if length(Iz)>Nmin 
        
        Iz1     = Iz(1);
        Iz2     = Iz(end);
        zfinite = z(Iz1:Iz2); %trim to the actual data
        Iz      = find(isnan(zfinite));
        Iz      =[1 Iz' length(zfinite)]';
        
        for isegment=1:length(Iz)-1  
            
            z       =zfinite(Iz(isegment):Iz(isegment+1)-1);
            
            if (length(z)>Nmin)
                [p,f]   =powerspect(z,dx);
            else
                p       =[NaN NaN];
                f       =[0 1000];
            end
            
        PowerStructx(i,isegment).px = p;
        PowerStructx(i,isegment).fx = f; 
        
        % interpolate onto a consistent f vector
        pxi                             =interp1(f,p,fxi);     
        PowerStructx(i,isegment).pxi    =double(pxi);
        end
    else
        p=[NaN NaN];
        f=[0 1000]; 
        PowerStructx(i,1).px = p;
        PowerStructx(i,1).fx = f;
        
        % interpolate onto a consistent f vector
        pxi                     =interp1(f,p,fxi);         
        PowerStructx(i,1).pxi   =double(pxi);
    end
end

fx=fxi;
end

%%%%%%%%%%%%%%spectral estimation function
function [p,f] = powerspect(z,dx)

% NoiseLevel=maximum reasonable value of 2nd deriv
N=length(z);

%     % check for large 2nd derivs
%     z2=diff(z,2)./dx.^2;
%     I=(abs(z2)>(4*std(abs(z2)))); % SMOOTHED ACCORDING TO TIBO
%     %I=(abs(z2)>NoiseLevel); % ORIGINAL SMOOTHING
%
%     while (sum(I)>0)
%
%         Iz=logical([0 I' 0]');% shift to the z equivalent
%         Ibefore=logical([I' 0 0]'); % index before
%         Iafter=logical([0 Iz(1:end-1)']'); %index after
%         z(Iz)=(z(Ibefore)+z(Iafter))/2;
%         z2=diff(z,2)./dx.^2;
%         I=(abs(z2)>(4*std(abs(z2))));
%     end
%     d
z=z-mean(z);
z=detrend(z);
z=taper(z,0.05,0.05);

y=fft(z);

% power
p=y.*conj(y)./(N*dx); %  /(N.*N);


% put back in dx
p=p.*dx*dx;
f=(0:N-1)'/(dx*N);
p=p(3:N/2);
f=f(3:N/2);

end




