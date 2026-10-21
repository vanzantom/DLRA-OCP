%%GenerateData permits to regenerate all data analyzed.
%The data will be saved in the folder data_output.
%Use the postprocessing script in data_output to regenerate the figures.

% **CAREFUL**: both ConvergenceOptimization and TimeOptimization can be
% very memory and time consuming

%Figure 1
AnalysisCovariantODE;

%Figure 2
ConvergenceGradients(1/32,500,1);

%Figure 3
GeometricAnalysis(1/32,500,1,1e-3);

%Figure 4
ConvergenceOptimization(1/32,500,1,1e-3);

%Figure 5
TimeOptimization(2^-4,2^7,1,1e-3);
TimeOptimization(2^-5,2^8,1,1e-3);
TimeOptimization(2^-6,2^9,1,1e-3);
