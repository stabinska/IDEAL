function printParams(s, name, indent)
% PRINTPARAMS  Recursively print all fields of a struct for verification.
%
%   printParams(s)              print struct s with default root label
%   printParams(s, name)        use name as the root label
%   printParams(s, name, indent) (used internally for recursion)
%
%   Example:
%       printParams(P, 'P')

    if nargin < 2 || isempty(name);   name   = inputname(1); end
    if nargin < 3 || isempty(indent); indent = 0;            end

    pad = repmat('  ', 1, indent);

    if indent == 0
        fprintf('=== %s ===\n', name);
    end

    fields = fieldnames(s);
    for i = 1:numel(fields)
        key = fields{i};
        val = s.(key);

        if isstruct(val)
            fprintf('%s%s:\n', pad, key);
            printParams(val, key, indent + 1);

        elseif isnumeric(val) || islogical(val)
            if isscalar(val)
                fprintf('%s%s: %g\n', pad, key, val);
            elseif isvector(val)
                fprintf('%s%s: [%s]\n', pad, key, num2str(val, '%g '));
            else
                % matrix — print dimensions then each row
                [nr, nc] = size(val);
                fprintf('%s%s: [%dx%d matrix]\n', pad, key, nr, nc);
                for r = 1:nr
                    fprintf('%s  row %d: [%s]\n', pad, r, num2str(val(r,:), '%g '));
                end
            end

        elseif ischar(val) || isstring(val)
            fprintf('%s%s: "%s"\n', pad, key, val);

        elseif iscell(val)
            fprintf('%s%s: {%d-element cell}\n', pad, key, numel(val));
            for c = 1:numel(val)
                item = val{c};
                if ischar(item) || isstring(item)
                    fprintf('%s  {%d}: "%s"\n', pad, c, item);
                else
                    fprintf('%s  {%d}: %s\n', pad, c, mat2str(item));
                end
            end

        else
            fprintf('%s%s: <unsupported type: %s>\n', pad, key, class(val));
        end
    end
end
