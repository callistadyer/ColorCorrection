function result = colorCorrectionValidateExamples(varargin)
% Run all the examples in the specified repo
%
% Syntax:
%     result = colorCorrectionValidateExamples()
%
% Description:
%     Run all the examples in the repository tree, excep those that contain
%     a line of the form 
%          "% ETTBSkip"
%
% Inputs:
%     None.
%
% Outputs:
%    result - describing the outcome
%      result.names   Names of the functions
%      result.status  What happened
%
% Optional key/value pairs
%    'select'      - 'all' or 'one' (one not yet implemented).  Default 'all'.
%    'print'       - Boolean. Print the results to the command line, showing the successes
%                    and failures separately. Default true.
%   'saveprint'    - logical, default true.  Save the report in a sensible
%                    place.

% History:
%   06/24/25  dhb, cmd  Wrote the repo specific version.

% Parse
p = inputParser;
p.addParameter('select','all',@ischar);
p.addParameter('print',true,@islogical);
p.addParameter('saveprint',true,@islogical);
p.parse(varargin{:});

%% Run all the examples
[result.names, result.status ] = ExecuteExamplesInDirectory(ColorCorrectionRootPath,'verbose',false);
outputBaseName = 'colorCorrection_examples';

if (p.Results.print)
    % Make sure save directory exists and set up save filenme, if saving
    if (p.Results.saveprint)
        outputDir = fullfile(ColorCorrectionRootPath,'validation','outputfiles');
        if (~exist(outputDir,'dir'))
            mkdir(outputDir);
        end
        outputFile = fullfile(outputDir,[outputBaseName '_' datestr(now,'yy-mm-dd-HH-MM-SS')]);
        outputFID = fopen(outputFile,'w');
    end

    % Maybe make this a function like examplesResultsPrint
    names = result.names;
    status = result.status;
    goodNames = names(status > 0);
    goodStatus = status(status > 0);
    badNames  = names(status < 0);

    cprintf('Blue','\n\nSuccess cases (%d)\n',numel(goodNames));
    if (p.Results.saveprint)
        fprintf(outputFID,'\n\nSuccess cases (%d)\n',numel(goodNames));
    end
    for ii=1:numel(goodNames)
        fprintf('%d examples suceeded for %s.\n',goodStatus(ii),goodNames{ii});
        if (p.Results.saveprint)
            fprintf(outputFID, '%d examples suceeded for %s.\n',goodStatus(ii),goodNames{ii});
        end
    end

    cprintf('Red','\n\nFailure cases (%d)\n',numel(badNames));
    if (p.Results.saveprint)
        fprintf(outputFID,'\n\nFailure cases (%d)\n',numel(badNames));
    end
    for ii=1:numel(badNames)
        fprintf('At least one example failed for %s.\n',badNames{ii})
         if (p.Results.saveprint)
            fprintf(outputFID, 'At least one example failed for %s.\n',badNames{ii});
        end
    end

    if (p.Results.saveprint)
        fclose(outputFID);
    end
end

%% END