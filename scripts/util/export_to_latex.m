function [out] = export_to_latex(name, value, format)
% EXPORT_TO_LATEX Export numbers to a newcommand statement for TeX
%
% Parameters:
%   name - name of the command in TeX
%   value - corresponding value
%   format ('raw' | sprintf-compatible)
%       if 'raw', value is printed as is, otherwise passed to num2str
%
% Returns:
%   out - TeX string

    if nargin < 3
        format = 'raw';
    end

    if ~strcmp(format, 'raw')
        value = num2str(value, format);
    end

    out = sprintf("\\newcommand{\\%s}{%s}", name, value);
end

