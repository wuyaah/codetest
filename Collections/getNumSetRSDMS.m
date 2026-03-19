function numSet = getNumSetRSDMS(dataPack, alphaScaling, numLoadCase, numVert, vertNorm, numTerms)
% Prepare numerical parameters struct for RSDMS
%
% Inputs:
%   dataPack      - Cell array containing element data structs
%   alphaScaling  - Scaling factor for load path amplitudes
%   numLoadCase   - Number of load cases
%   numVert       - Number of vertices used for load path discretization
%   numTerms      - Number of Fourier series terms for decomposition
%
% Output:
%   params        - Struct containing numerical parameter settings

numSet.numVert = numVert;                    % Number of vertices for load path approximation
numSet.vertNorm = vertNorm;
numSet.numElem = numel(dataPack);            % Total number of elements in the model
numSet.numStrComp = dataPack{1}.numStrComp;  % Number of strain components per element

numSet.NGK = dataPack{1}.NGK;                 % Number of Gauss integration points per element
numSet.numTotalGauss = numSet.numElem * numSet.NGK;  % Total Gauss points in the model
numSet.dimension = dataPack{1}.dimension;     % Spatial dimension of the model (2D or 3D)
numSet.Yield=dataPack{1}.Yield;

numSet.numTerms = numTerms;                    % Number of terms in Fourier decomposition
numSet.numLoadCase = numLoadCase;              % Number of different load cases
numSet.numGL = numVert + 1;                    % Number of Gauss-Legendre integration points (vertices + 1)
numSet.alphaScaling = alphaScaling;            % Scaling factor for load factor

% Construct loading path function using trigonometric basis
% pathFunc = LoadingPathFunc(numVert);           % Trigonometric loading path functions
% Alternative: polynomial basis (commented)
% pathFunc = LoadingPathFunc2(numVert);

% numSet.pathFunc = pathFunc;

end
