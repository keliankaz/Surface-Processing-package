function  parameterStruct = surface_analysis(surfaceGrid, pointSpacing, ...
                                             numberOfScales, decimationFactor)

% the grunt of the analysis is done here...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run surface parameter analysis (std, skew, kurt and asymetry) - set at
% 10 sampled scales

parameterStruct     = surface_parameters(surfaceGrid, pointSpacing, ...
                                         numberOfScales);

% run surface fft analysis (with default decimation factor = 20)

if strcmp(decimationFactor,'default')
    decimationFactor = 1;
end

[fx1, PowerStructx] =  frequency_spectrum(surfaceGrid, pointSpacing, ...
                                          decimationFactor);

% process output - concatonate the structure array  (I think)
[Nx,Ny]             = size(PowerStructx);
PowerStructx        = reshape(PowerStructx,1,Nx*Ny);
[~,Px1]             = FindErr_loop_aniso(PowerStructx);

parameterStruct.FFT = {fx1, Px1};

% run plomb analysis - set at 300 sample frequencies

[Px2, fx2]          = fault_spectral_density_simple(surfaceGrid, ...
    pointSpacing, 300);

parameterStruct.PLOMB = {fx2, Px2};

end
