function loadif(path,ivar)
% Check if a variable is loaded into the workspace already. If not, load
% it. Useful for working with large variables in .mat files that you don't 
% want to load repeatedly. 
%
% path = path to .mat file
% ivar = char, name of variable

narginchk(2,2)

varloaded = evalin('caller', ['exist(''' ivar ''',''var'');']);

    if ~varloaded
        fprintf('Loading variable: %s\n',ivar)
        load(path,ivar)
        assignin('caller',ivar,eval(ivar))
%         varargout{1} = eval(ivar);
    else
        fprintf('Variable: ''%s'' already loaded.\n',ivar)
        return
    end
end