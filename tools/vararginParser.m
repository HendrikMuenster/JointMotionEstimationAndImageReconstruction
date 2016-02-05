% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 2.0
% Date: 2015-12-09

% All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author and University of Muenster not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.

%iterate through the number of additional input arguments and assign them according to the given input name
% additional arguments should be of the form function(arg1,...,argn,'varname',value);
for i=1:numel(varargin)
    if (mod(i,2)==1)
        assignin('caller',varargin{i}, varargin{i+1})
    end
end
