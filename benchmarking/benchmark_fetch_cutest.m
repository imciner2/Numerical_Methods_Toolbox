function [ problemDir ] = benchmark_fetch_cutest( problemName, varargin )
%BENCHMARK_FETCH_CUTEST Fetch a problem from the CUTEst benchmark set
%
% This function will fetch a problem from the CUTEst benchmark
% set (ftp://ftp.numerical.rl.ac.uk/pub/cutest/sif/mastsif.html), then
% convert it into the MATLAB problem format from CUTEst and place the
% problem in the active directory.
%
%
% Usage:
%   [ problemDir ] = benchmark_fetch_marosmeszaros( problemName )
%   [ problemDir ] = benchmark_fetch_marosmeszaros( problemName, SIFParams )
%
% Inputs:
%   problemName - The name of the problem in the CUTEst benchmark set to
%                 download.
%   SIFParams   - Parameters that should be passed to the SIF decoder to
%                 create the problem (a character string)
%
% Outputs:
%   problemDir - The file directory (full path) where the mcutest.mex file
%                is located. This directory must be the location the
%                problem is accessed from (otherwise it will fail).
%
% Created by: Ian McInerney
% Created on: November 28, 2017
% Version: 1.0
% Last Modified: November 28, 2017
%
% Revision History:
%   1.0 - Initial Release

disp(['Fetching problem ', upper(problemName)]);

%% Find the directory where the problems will be located, and navigate to it
scriptDir = which('benchmark_fetch_cutest');
scriptDir = strrep(scriptDir, 'benchmark_fetch_cutest.m', '');
originalDir = cd(scriptDir);


%% Call the script to get the problem and build it
disp('Calling problem fetcher');
cd('scripts');

if ( ~isempty(varargin) )
    params = varargin{1};
    % Give parameters to the SIF decoder
    system(['./getSIF.sh ', problemName, ' cutest ftp://ftp.numerical.rl.ac.uk/pub/cutest/sif ', params]);
else
    % No parameters needed
    system(['./getSIF.sh ', problemName, ' cutest ftp://ftp.numerical.rl.ac.uk/pub/cutest/sif']);
end
cd('../');

%% Navigate to the active directory for the problem, get it, and save it
cd('problems/activeCUTEst');


%% Go back to the original directory and save the problem dir
problemDir = cd(originalDir);

end