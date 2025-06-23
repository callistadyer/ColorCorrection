% ColorCorrection
%
% Template for setting preferences and other configuration things.

%% Define project
projectName = 'ColorCorrection';

%% Clear out old preferences
if (ispref(projectName))
    rmpref(projectName);
end

%% Specify project location
%
% TbTb can tell us where the project lives on the
% local directory, should we want to know.
projectBaseDir = tbLocateProject(projectName);

% Figure out where baseDir for other kinds of data files is.
%
% Can only do this when we have GetComputerInfo available.
if (exist('GetComputerInfo','file'))
    sysInfo = GetComputerInfo();
    switch (sysInfo.localHostName)
 
        otherwise
            % Some unspecified machine, try user specific customization
            switch(sysInfo.userShortName)
                % Could put user specific things in, but at the moment generic
                % is good enough.
                otherwise
                    % Some unspecified machine, try our generic approach,
                    % which works on a mac to find the dropbox for business
                    % path.
                    if ismac
                        dbJsonConfigFile = '~/.dropbox/info.json';
                        fid = fopen(dbJsonConfigFile);
                        raw = fread(fid,inf);
                        str = char(raw');
                        fclose(fid);
                        val = jsondecode(str);
                        baseDir = val.business.path;
                    end
            end
    end
end

% Set up pointer to where we want to write the output.
setpref(projectName,'outputDir',fullfile(baseDir,'DALT_analysis',projectName));
if (~exist(getpref(projectName,'outputDir'),'dir'))
    mkdir(getpref(projectName,'outputDir'));
end







