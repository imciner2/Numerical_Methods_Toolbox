function [ QP ] = benchmark_fetch_marosmeszaros( problemName )
%BENCHMARK_FETCH_MAROSMESZAROS Fetch a problem from the Maros & Meszaros QP benchmark set
%
% This function will fetch a problem from the Maros & Meszaros benchmark
% set (http://www.cuter.rl.ac.uk/Problems/marmes.shtml), then convert it
% into a QP structure in MATLAB. The structure contains fields that map to
% the otpimization problem:
%
%   min x'Qx + c'x
%   s.t. 
%       Aeq*x= beq
%       Ale*x <= ble
%
%       lb <= x <= ub
%
% with intial conditions for the primal variables in x0, and initial
% conditions for the dual variables in mu0.
%
%
% Usage:
%   [ QP ] = benchmark_fetch_marosmeszaros( problemName )
%
% Inputs:
%   problemName - The name of the problem in the Maros & Meszaros QP
%                 benchmark set to download.
%
% Outputs:
%   QP - A structure containing the matrices and vectors that make up the
%        QP problem.
%
% Created by: Ian McInerney
% Created on: November 24, 2017
% Version: 1.0
% Last Modified: November 24, 2017
%
% Revision History:
%   1.0 - Initial Release

disp(['Fetching problem ', upper(problemName)]);

%% Find the directory where the problems will be located, and navigate to it
scriptDir = which('benchmark_fetch_marosmeszaros');
scriptDir = strrep(scriptDir, 'benchmark_fetch_marosmeszaros.m', '');
originalDir = cd(scriptDir);


%% Check to see if the problem exists as a .mat file already
saveFile=['problems/marosmeszaros/', upper(problemName), '.mat'];
if ( exist(saveFile, 'file') ~= 0 )
    % The file exists, load it then leave
    disp('Problem already downloaded and converted, loading .mat file');
    load(saveFile);
    cd(originalDir);
    return;
end


%% Call the script to get the problem
disp('Problem not converted, calling conversion scripts');
cd('scripts');
system(['./getMarosMeszaros.sh ', problemName]);
cd('../');

%% Navigate to the active directory for the problem, get it, and save it
cd('problems/activeCUTEst');
disp('Getting QP form of the problem');
[QP.Q, QP.c, QP.Aeq, QP.beq, QP.Ale, QP.ble, QP.ub, QP.lb, QP.x0, QP.mu0] = cutest_getQP();
cd(scriptDir);
save(saveFile, 'QP');


%% Go back to the original directory
cd(originalDir);

end