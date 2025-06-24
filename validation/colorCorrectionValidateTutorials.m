function report = colorCorrectionValidateTutorials(varargin)
% Run all tutorials/scripts for a repo and print out which worked and which did not
%
% Syntax:
%    report = colorCorrectionValidateTutorials()
%
% Description:
%   Run all of the tutorials/scripts/validations for a specified repo
%   and print out a report at the end as to whether they threw errors,
%   or not.  The report is also returned.
%
%   The path to the tutorials is setup in the source of this routine, as
%   are various strings that will cause something to be skipped if it is in
%   the pathname (e.g., paths with 'development' are skipped).
%
%   Particular tutorials/scripts will also be skipped if they
%   contain a comment line of the form
%       % UTTBSkip
%   in their source file.
%
%   Also, this script contains a list of directories/tutorials/scripts that
%   can be skipped.
%
% Inputs:
%   None.
%
%   See below.
%
% Outputs
%  report
%
% Optional key/value pairs
%   'saveprint'    - logical, default true.  Save the report in a sensible
%                    place.

% Examples:
%{
   % ETTBSkip
   % 
   % You need to make sure you have the underlying repos on your path
   % these will work.

    colororrectionValidate('tutorials');
%}

% History:
%   06/23/25  dhb, cmd  Wrote the repo specific version.

%% Specify my repo
projectName = 'ColorCorrection';
repo = projectName;
typeToRun = 'tutorials';

%% Specify repos we can test.  
% 
% For each, also need to provide
% the name of the function that gives the repository root path
% and the path to the tutorial directory under that path.
availRepos = {projectName};
repoRootDirFcns = {'ColorCorrectionRootPath'};

p = inputParser;
p.addParameter('saveprint',true,@islogical);
p.parse(varargin{:});

% Figure out where we want to go today
knownRepo = false;
for rr = 1:length(availRepos)
    if (strcmp(repo,availRepos{rr}))
        knownRepo = true;
        selectedRepoNum = rr;
        break;
    end
end
if (~knownRepo)
    error('Unknown repository requested')
end

% Choose the top level directory corresponding to the type of
% validation.  For repo and validation type, we specify the
% top level directory of the scripts to run, and the subdirectory
% of that top level directory where the scripts are.
switch (typeToRun)
    case 'tutorials'
        topLevelDir = eval(repoRootDirFcns{selectedRepoNum});
        subDir = 'tutorials';
        outputFileBase = 'colorCorrection_tutorials';
    % case 'scripts'
    %     topLevelDir = eval(repoRootDirFcns{selectedRepoNum});
    %     subDir = 'scripts';
    %     outputFileBase = 'colorCorrection_scripts';
end

% Set up preferences to work for the selected run
pp = struct(...
    'tutorialsSourceDir',       fullfile(topLevelDir,subDir) ...
    );

%% List of scripts to be skipped from automatic running.
%
% Anything with one of these strings in the path name is skipped.
% library is in ISETBio/scripts.  It creates all the mosaics in a library.
scriptsToSkip = {...
    'Contents', ...
    'data', ...
    'deprecated', ...
    'development', ...
    'Development', ...
    'xNotOnPath', ...
    };

%% Use UnitTestToolbox method to do this.
[~, reportTemp] = UnitTest.runProjectTutorials(pp, scriptsToSkip, 'All');

% This has the effect that if you call the function from the command
% line and don't assign its output to anything and don't put in a
% semi-colon, you don't get the report string dumped out.  The reason we
% want this is that the report string is already printed out with the
% broken ones in color, and if we print it out again that version rolls off
% the screen.
if (nargout > 0)
    report = reportTemp;
end

% Save what happend to a file
if (p.Results.saveprint)
    % Make sure save directory exists and set up save filenme, if saving
    outputDir = fullfile(ColorCorrectionRootPath,'validation','outputfiles');
    if (~exist(outputDir,'dir'))
        mkdir(outputDir);
    end
    outputFile = fullfile(outputDir,[outputFileBase '_' datestr(now,'yy-mm-dd-HH-MM-SS')]);
    outputFID = fopen(outputFile,'w');
    fprintf(outputFID,reportTemp);
    fclose(outputFID);
end

end