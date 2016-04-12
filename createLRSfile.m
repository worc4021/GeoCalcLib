function createLRSfile(varargin)
% This function creates a text file to be read by the LRS executable.
% createLRSfile('param',param,...) accepts the following parameters:
% 'fname' - Filename of output file, (default lrstest.ine.)
% 'Rep' - {'H','V'}  'H' for vertex enumeration and 'V' facet enumeration,
% default 'V'.
% 'V' - list of vertices/rays, only accepted if Rep = 'V'
% 'Type' - list of 0s for rays and 1s for vertices. V(i,:) is a vertex if
% type(i) == 1.
% 'A' - A matrix for {A*x<=b}, only accepted if Rep = 'H'
% 'b' - b vector, if not provided {A*x<=1}
% 'tol' - tolerance for rational approximation, default 1e-12.


p = inputParser;
expectedReps = {'V','H'};


addOptional(p,'Rep','V',@(x) any(validatestring(x,expectedReps)) );
addOptional(p,'fname','lrstest.ine',@(x) validateattributes(x,{'char'},{'nonempty'}) );
addOptional(p,'V',zeros(0,1),@(x) validateattributes(x,{'numeric'},{'nonempty'}));
addOptional(p,'A',zeros(0,1),@(x) validateattributes(x,{'numeric'},{'nonempty'}));
addOptional(p,'Type',ones(0,1) ,@(x) validateattributes(x,{'numeric'},{'nonempty','ncols',1}));
addOptional(p,'b',ones(0,1) ,@(x) validateattributes(x,{'numeric'},{'nonempty','ncols',1}));
addOptional(p,'tol',1e-12, @isnumeric );

parse(p,varargin{:});


rep = p.Results.Rep;

if strcmp(rep,'V')
    if isempty(p.Results.V)
        error('No vertices passed.')
    else
        dat = p.Results.V;
    end
    
    if isempty(p.Results.Type)
        display('No type passed. Input is taken as vertices.')
        t = ones(size(dat,1),1);
    elseif size(dat,1) ~= size(p.Results.Type,1)
        error('Type must have as many entries as V has rows.');
    else
        t = p.Results.Type;
    end
else
    if isempty(p.Results.A)
        error('No vertices passed.')
    else
        dat = -p.Results.A;
    end
    
    if isempty(p.Results.b)
        display('No b vector passed. Assumed as all ones.');
        t = ones(size(dat,1),1);
    elseif size(dat,1) ~= size(p.Results.b,1)
        error('A and b must have the same number of rows.');
    else
        t = p.Results.b;
    end
end

fname = p.Results.fname;
tol = p.Results.tol;

data = [t,dat];

fid = fopen(fname,'w+');

if strcmp(rep,'V')
    fprintf(fid,'V-representation\n');
else
    fprintf(fid,'H-representation\n');
end

fprintf(fid,'begin\n');
fprintf(fid,'%d %d rational\n', size(data,1), size(data,2));

for i = 1:size(data,1)
    for j = 1:size(data,2)
        [n,d] = rat(data(i,j),tol);
        fprintf(fid,'%s/%s ',num2str(n),num2str(d));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'end');
fclose(fid);

