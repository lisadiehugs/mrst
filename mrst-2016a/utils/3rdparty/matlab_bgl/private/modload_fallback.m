function modload
   d = fileparts(fileparts(mfilename('fullpath')));

   error(['MatlabBGL software is not available\n\n',    ...
          ' -> Use downloadMBGL script to install ',    ...
          '(copy-paste command below)\n\n',             ...
          '    run ''%s''\n\n',                         ...
          'Then run mrstModule add matlab_bgl again.'], ...
          fullfile(d, 'downloadMBGL'));
end