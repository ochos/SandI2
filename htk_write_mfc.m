function [status, message] = htk_write_mfc_new(filename, nSamples, sampPeriod, sampSize, parmKind, data)
% function [status, message] = htk_write_mfc_new(filename, nSamples, sampPeriod, sampSize, parmKind, data)
% 
% MATLAB function to write HTK files .mfc files
%
% Inputs:  filename   = name of file to read
%          nSamples   = number of samples in file
%          sampPeriod = sampling period in 100ns units (e.g. 10ms = 100000)
%          sampSize   = 4 * number of vectors used (e.g. 4*14 = 56)
%          parmKind   = type of HTK file (e.g. 8262 for mfc)
%          data       = data stored as nSamples columns, each column being
%                       a feature vector of size (sampSize / 4)
%
% Outputs: status     = status is 2 when an error has occurred
%          message    = reason for error
%
% Jonathan Darch 28/10/03
%
% University of East Anglia, UK
% jonathan.darch@uea.ac.uk
% www.jonathandarch.co.uk

if nargin ~= 6
	help(mfilename)
	return
end

% Obtain dimensions of data:
sizeData = size(data);

% Check that number of ROWS in data is equal to nSamples:
nSamplesData = sizeData(2);

if nSamples ~= nSamplesData
    message = 'ERROR: nSamples does not match data. Terminating...';
    status = 2;
    return
else
    message = '';
    status = 1;
end


% Check that number of COLUMNS in data is a quarter of sampSize:
numVecHead = sampSize / 4;                  % As each vector entry is a 4 byte float.

numVecData = sizeData(1);

if numVecHead ~= numVecData
    message = 'ERROR: sampSize does not match data. Terminating...';
    status = 2;
    return
else
    message = '';
    status = 1;
end


% Open file for writing:
fid = fopen(filename, 'w', 'ieee-be');

% Write the header information
% ==============================
fwrite(fid, nSamples, 'int32');         % number of samples in file (4 byte int):
fwrite(fid, sampPeriod, 'int32');       % sample period in 100ns units (4 byte int):
fwrite(fid, sampSize, 'int16');         % number of bytes per sample (2 byte int):
fwrite(fid, parmKind, 'int16');         % code indicating the sample kind (2 byte int):


% Write the data:
% ================

% Write one sample at a time:
for i = 1:nSamples 
        fwrite(fid, data(:, i), 'float32');
end

% Close the file:
fclose(fid);
