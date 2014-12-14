function [ DCTs ] = featureExtractionNoisy( filename, frameTimeLength, HMMtype)
%featureExtraction Summary of this function goes here
%   Input Arguments:
%             filename = .wav File to be import to matlab.
%                   Fs = Sampleing Frequency of the Audio Signal.
%      frameTimeLength = time of each short time frame.

    % Import speech and noiseFile
    [speechFile,Fs] = wavread(filename);                            % Import the Audio Signal.
    speechFile = speechFile(:,1);                                   % Since the Audio Signal is recorded in 2 channels, it is trimmed to single channel.
    
    speechFile = filter([1,- 0.97],1,speechFile);                   % Pre-emphasis

    % Declare paramteres to create short time frame
    frameSize = frameTimeLength * 10^-3;                            % (20ms)
    frameShift = 10 * 10^-3;                                        % overlapping (10ms)
    frameSizeSample = fix( Fs * frameSize );                        % frame length
    frameShiftSample = fix( Fs * frameShift );                      % frame shift length,  difference is 80 by default
    % calculate number of frame (max)
    numFrames = fix( ( length(speechFile) - (frameSizeSample - frameShiftSample) ) / frameShiftSample ) - 1;
    
    % Create hamming Window.
    hamWindow = hamming(frameSizeSample);                           % Hamming Window in the same length of short time frames.

    % Declare parameters to create Square Wave filters.
    numFilters = 25;                                                % Total Number of filters.
    filterLength = floor(frameSizeSample /numFilters);              % Length of each filters.
    filterShift = filterLength*0.5;                                 % overlapping (50%)
    
    
    numNoiseEstimate = 10;
    
    % pre-allocate 
    result = zeros(1,numFilters);
    filterBank = zeros(numFrames,numFilters);
    DCTs = zeros(numFrames,numFilters);
    logE = zeros(numFrames,1);
    noisyMagEST = zeros(frameShiftSample,numNoiseEstimate);
    
    % create mel-filter 
    [myMel] = melfilter(Fs,frameSizeSample,numFilters);
    
  %-------------------------------------------------------------------------------------------%  
    for frame = 1 : numFrames,
         [ cleanMagEST, energy ] = SpeechProcess(frame, frameSizeSample, frameShiftSample, hamWindow, speechFile);
         
         
        if (frame < 11)
            noisyMagEST ( : , frame ) = cleanMagEST; 
        else
            noiseEST = mean(noisyMagEST,2);
        
            cleanMagEST = cleanMagEST - noiseEST;
        end
        
        for nFilter = 1: numFilters
%             [filteredResult] = myFilter( filterLength, magnitude, nFilter );
           [filteredResult] = myMelFilter(myMel, cleanMagEST, nFilter );
            result(nFilter) = filteredResult;                   % A Matrix storing the filtered Result
        end
        
            filterBank ( frame , :) = result;
            
            DCTs ( frame , :) = dct( log( filterBank(frame,:) ) );
%              DCTs ( frame , :) = dct( ( filterBank(frame,:) ) );

            logE(frame,1) = log(energy);
    end
    
     tc = 14;                                                    % Truncation value
    DCTs = DCTs(1:end, 1:tc);

    if (HMMtype == 1)                                           % create MFCC

        [DCTsRow, DCTsCol] = size(DCTs);
        
        nSamples = DCTsRow;
        sampPeriod  = frameSize*10^3 * 10000/2;
        vecSize = DCTsCol;
        sampSize = vecSize*4;
        parmKind = 6;     

%         HTKwriteMFC( nSamples, sampPeriod, sampSize, parmKind, vecSize, DCTs, mfcName );
        
    elseif (HMMtype == 2)                                       % create MFCC_E 
        
        DCTs = [DCTs logE];                                     % append logE component
   
        [DCTsRow, DCTsCol] = size(DCTs);

        nSamples = DCTsRow;
        sampPeriod  = frameSize*10^3 * 10000/2;
        vecSize = DCTsCol;
        sampSize = vecSize*4;
        parmKind = 6;     
        parmKind = bitset(parmKind, 7);                         % Set the energy bit

%         HTKwriteMFC( nSamples, sampPeriod, sampSize, parmKind, vecSize, DCTs, mfcName );
    
    elseif (HMMtype == 3)                                       % create MFCC_E_D
        
        tc = 14;                                                % Truncation value
        DCTs = DCTs(1:end, 1:tc);
        DCTs = [DCTs logE];                                     % append logE component
    
        % mean variabce normalise
        DCTs = mvnorm(DCTs);
    
        DCTsize = size(DCTs);
        velocity = [DCTs(:,2:end) zeros(DCTsize(1),1)] -[zeros(DCTsize(1),1) DCTs(:,1:end-1)];
        DCTs = [DCTs velocity];
    
        [DCTsRow, DCTsCol] = size(DCTs);

        nSamples = DCTsRow;
        sampPeriod  = frameSize*10^3 * 10000/2;
        vecSize = DCTsCol;
        sampSize = vecSize*4;
        parmKind = 6;     
        parmKind = bitset(parmKind, 7);                         % Set the energy bit
        parmKind = bitset(parmKind, 9);                         % Set the velocity bit
        
%         HTKwriteMFC( nSamples, sampPeriod, sampSize, parmKind, vecSize, DCTs, mfcName );
                
    elseif (HMMtype == 4)                                       % create MFCC_E_D_A
        
        tc = 14;                                                % Truncation value
        DCTs = DCTs(1:end, 1:tc);
        DCTs = [DCTs logE];                                     % append logE component
    
        % mean variabce normalise
        DCTs = mvnorm(DCTs);
    
        DCTsize = size(DCTs);
        velocity = [DCTs(:,2:end) zeros(DCTsize(1),1)] -[zeros(DCTsize(1),1) DCTs(:,1:end-1)];
        acceleratio = [velocity(:,2:end) zeros(DCTsize(1),1)] -[zeros(DCTsize(1),1) velocity(:,1:end-1)];
        DCTs = [DCTs velocity acceleratio];
    
        [DCTsRow, DCTsCol] = size(DCTs);

        nSamples = DCTsRow;
        sampPeriod  = frameSize*10^3 * 10000/2;
        vecSize = DCTsCol;
        sampSize = vecSize*4;
        parmKind = 6;     
        parmKind = bitset(parmKind, 7);                         % Set the energy bit
        parmKind = bitset(parmKind, 9);                         % Set the velocity bit
        parmKind = bitset(parmKind, 10);                        % Set the acceleration bit
        
      %  HTKwriteMFC( nSamples, sampPeriod, sampSize, parmKind, vecSize, DCTs, mfcName );
        
    end

end


function myMel =melfilter(fs1,fsize1,noc1)

    % maximum frequency 
    fmax = fs1/2; 
    % maximum frame index
    Nmax = fsize1/2; 
    % maximum mel frequency
    melmax = hz2mel(fmax); 
    % frequency resolution
    df = fs1/fsize1; 
    
    % frequency increment on mel scale
    dmel = melmax / (noc1 + 1); 
    
    %center frequencies on mel scale
    melcenters = (1:noc1) .* dmel;
    
    %center frequencies in linear scale [Hz]
    fcenters = mel2hz(melcenters);
    
    %on linear scale gives corresponding indices of coefficient of power
    %spectrum ie. fft indices
    indexcenter = round(fcenters ./df);          %center of nc(k)
    
    %compute start indices of windows
    indexstart = [1 , indexcenter(1:noc1-1)];   %start of nc(k)= center of nc(k-1)
    
    %compute stop indices of windows
    indexstop = [indexcenter(2:noc1),Nmax];     %stop of nc(k)=center of nc(k+1)
    
    
    %compute triangle-shaped filter coefficients
    myMel(1:noc1,1:Nmax) =0;
    for c = 1:noc1
        
        %left ramp
        increment = 1.0/(indexcenter(c) - indexstart(c));     %maximum amplitude is 1
            for i = indexstart(c)+1:indexcenter(c)
                myMel(c,i) = (i - indexstart(c))*increment;
            end
            
        %right ramp
        decrement = 1.0/(indexstop(c) - indexcenter(c));
            for i = indexcenter(c):indexstop(c)
                myMel(c,i) = 1.0 - ((i - indexcenter(c))*decrement);
            end
    end
    
myMel=nanclr(myMel); % to replace  'not a number' vlues

% %figure;
% plot(myMel');
end



function mel =  hz2mel(f)
    % convert Hz into mel
    mel = 2595 * log10( 1 + ( f / 700 ) ); 
end


function freq = mel2hz(m)
    % convert mel into Hz
    freq = 700*((10.^(m ./2595)) -1); 
end


% to replace  'not a number' vlues
function w1=nanclr(w1)  
[r1,c1]=size(w1);
for i=1:r1
    temp=isnan(w1(i,1:c1));
    for j=1:c1
        if temp(j)==1
            w1(i,j)=0.0001;
        else
        end;
    end;
end;

end



function [filteredResult] = myMelFilter(melFilter, magnitude, nFilter )

    filteredMag = magnitude .* melFilter(nFilter,:)';    % Apply the mel filter to Magnitude Spectrum.
    
    filteredResult = sum(filteredMag);

end

function [ magnitude, energy ] = SpeechProcess( frame, frameSizeSample,frameShiftSample, hamWindow, speechFile )
%SpeechProcess returns values of magnitude and phase.
%   Input:
%       frame = frame counts from the for-loop.
%       frameSizeSample = number of samples in each frame.
%       frameShiftSample = overlapping value. 0.010 = (10ms)
%       hamWindow = Hamming Window in the same length of short time frames.
%       speechFile = the Audio Signal.

%     startThisFrame = 0.0;
%    
%     if (frame == 1 )
%         % process each frame
%         startThisFrame = 1;                                   % start sampling number 
%         endThisFrame = startThisFrame + frameSizeSample - 1;  % last sampling number
%     else
%         startThisFrame = startThisFrame + frameShiftSample;   %Determine starting point of each short time frame.
%         endThisFrame = startThisFrame + frameSizeSample - 1;  %Determine end point of each short time frame.
%     end

    sample1st = frameSizeSample/2 * (frame - 1 ) + 1;
    sampleEnd = frameSizeSample + sample1st - 1;
    
    shortTimeFrame = speechFile(sample1st : sampleEnd);     %Extract signal from speechFile, 
                                                            %from point sample1st to sampleEnd.
    hamFrame = shortTimeFrame .* hamWindow;              	%Applying Hamming window to the signal.
    result = fft(hamFrame);                                 %Fast Fourier Synthesis.
    mag = abs(result);                                      %Result of magnitude spectrum, by getting rid of imaginary parts.
    magnitude = mag(1 : floor( length(mag)/2) );            %Due to the symmetric characteristic of magnitude spectrum, it is trimmed to first half.
    phase = angle(result);                                  %Phase.
    energy = sum(mag,1);                                    %Calculate log energy
end

function [ mvfeat ] = mvnorm( feat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sfeat = size(feat);
m = sum(feat,2)/sfeat(2);
mfeat = feat-repmat(m,1,sfeat(2));
sd = sqrt(sum(mfeat.^2,2)/sfeat(2));
mvfeat = mfeat./repmat(sd,1,sfeat(2));
end
