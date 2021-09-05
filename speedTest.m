
% To find dependencies
[flist, plist] = matlab.codetools.requiredFilesAndProducts(...
    'generatePoisson2D.m');

%%
tic;
for i = 1:500
    % Simulate the process
    testBoundaryEstimation;
end
toc

