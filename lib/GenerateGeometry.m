function [ geo ] = GenerateGeometry(h, data, para)
%%GENERATE_GEOMETRY assembles the structure geo describing the FE mesh.
%
%Inputs: -h mesh size parameter
%        -data structure describing the problem to be solved (containing e.g., boundary conditions..)
%        -para structure which contains other parameters (e.g., order of FE space)
%Outputs:-geo structure, see package GDGMatlab for a full description

[ geo ] = geometry(h, data, para);

filename = sprintf('data_input/geo_%d.mat', 1/h);
save( filename, 'geo' )

end