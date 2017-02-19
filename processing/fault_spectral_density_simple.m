function [PXX, fvec] = fault_spectral_density_simple(Z,dx,n,varargin)

% Function returns the 1-D power specrtal density of the fault topography
% for an N by M matrix Z along increasing M for corresponding x  coordinates in X (N
% by M arrray). Spectral estimation uses the Lomb-Scargle periodogram.

% Input is a two dimensional array. The Lomb-Scargle periodogram
% effectively used a least squares astimator to evaluate the spectral
% density similar to a FFT algorithm with the added benifit on not
% requiring evenly spaced data. X in an N by M array with the corresponding
% X coordinates (along M) associated to Z.

% varargin - fvec specify the frequency vector at which the power will be
% computed

% remove entire rows of nan
Z = Z(sum(~isnan(Z),2)>10, :);

[N,M] = size(Z);
X = (1:M)*dx;


for rowCount = 1:N

    %% detrend
    % find best fit to data
    row         = Z(rowCount,:);
    goodDataLoc = ~isnan(row);
    goodData    = row(goodDataLoc);
    goodX       = X(goodDataLoc);

    xbar        = mean(goodX);
    zbar        = mean(goodData);
    numVal      = length(goodX);
    slope       = (sum(goodX.*goodData)-numVal*xbar*zbar)/ ...
                (sum(goodX.^2)-numVal*(xbar^2));
    
    % remove the linear trend in the row
    row(goodDataLoc) = goodData - slope*goodX;
      
    % remove mean value of the row
    row         = row-nanmean(row);
    
    %% Taper
    % create a 5% window
    w           = tukeywin(M,0.1); % built in function that creates a window
    
    % apply the window onto the row
    Z(rowCount,:)      = row.*w';
end

%% run built in Lomb-Scargle periodogram function

% create the query frequency vector

if nargin == 3
    fvec = logspace(log10(1/(M*dx)),log10(1/(2*dx)),n);
    
elseif nargin == 4
    fvec = varargin{1}; 
end
    

% run the function
[pxx, fvec] = plomb(Z',X,fvec);

% average across all rows
PXX = mean(pxx,2);


end 

