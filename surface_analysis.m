function  parameterStruct = surface_analysis(surfaceGrid, pointSpacing, ...
                                             numberOfScales, ...
                                             decimationFactor, toDo)

% the grunt of the analysis is done here...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run surface parameter analysis (std, skew, kurt and asymetry) - set at
% 10 sampled scales

if sum(strcmp(toDo,'parameters'))>=1 || strcmp(toDo,'all')

parameterStruct     = surface_parameters(surfaceGrid, pointSpacing, ...
                                         numberOfScales);
end                       

% run surface fft analysis (with default decimation factor = 1)

if strcmp(decimationFactor,'default')
    decimationFactor = 1;
end

if sum(strcmp(toDo,'FFT'))>=1 || strcmp(toDo,'all')
    [fx1, PowerStructx] =  frequency_spectrum(surfaceGrid, pointSpacing, ...
                                              decimationFactor);

    % process output - concatonate the structure array  (I think)
    [Nx,Ny]             = size(PowerStructx);
    PowerStructx        = reshape(PowerStructx,1,Nx*Ny);
    [~,Px1]             = FindErr_loop_aniso(PowerStructx);

    parameterStruct.FFT = {fx1, Px1};

end
    
% run plomb analysis - set at 300 sample frequencies

if sum(strcmp(toDo,'PLOMB'))>=1 || strcmp(toDo,'all')

    [Px2, fx2]          = fault_spectral_density_simple(surfaceGrid, ...
        pointSpacing, 300);

    parameterStruct.PLOMB = {fx2, Px2};

end

end
