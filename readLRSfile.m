function varargout = readLRSfile(fname)
% readLRSfile(filename) reads the data from an LRS input file.
% returns [A,b] as {A*x<=b} for H-representations and [V,type] for 
% V-representations, where type(i) = 1 if V(i,:) is a vertex and type(i)=0
% if V(i,:) is a ray.

    validateattributes(fname,{'char'},{'nonempty'});

    fid = fopen(fname,'r');

    while 1
       curStr = returnOneLine(fid);
       if any(strfind(curStr,'H-representation'))
           rep = 'H';
       end

       if any(strfind(curStr,'V-representation'))
           rep = 'V';
       end

       if any(strfind(curStr, 'rational'))
          spacePos = strfind(curStr,',');
          m = str2num(curStr(1:(spacePos(1)-1)));
          n = str2num(curStr((spacePos(1)+1):spacePos(2)-1));
          break;
       end
    end

    dat = zeros(m,n);

    for i = 1:m
        curStr = returnOneLine(fid);
        curStr = breakStringDown(curStr);
        dat(i,:) = convertToDouble(curStr);
    end
    
    fclose(fid);

    if strcmp(rep,'V')
        varargout{1} = dat(:,2:end);
        varargout{2} = dat(:,1);
    else
        varargout{1} = -dat(:,2:end);
        varargout{2} = dat(:,1);
    end
        
end

function string = returnOneLine(fileID)
string = '';
    while 1
        cur = fread(fileID,1,'*char');
        if strcmp(char(10),cur)
            break;
        end
        if strcmp(char(32),cur)
            string = strcat(string,',');
        else
            string = strcat(string,cur);
        end

    end
end

function out = breakStringDown(string)

    spacePos = [0,strfind(string,','),length(string)+1];
    n = numel(spacePos)-1;
    out = cell(n,1);
    for i =1:n
        out{i} = string((spacePos(i)+1):spacePos(i+1)-1);
    end
end

function out = convertToDouble(stringCell)
    n = length(stringCell);
    out = zeros(1,n);
    if strcmp(stringCell{1},'0')
        for i = 1:n
            out(i) = str2num(stringCell{i});
        end
        out = out/norm(out);
    else
        if strcmp(stringCell{1}(1),'-')
            out(1) = -1;
            lcd = stringCell{1}(2:end);
        else
            out(1) = 1;
            lcd = stringCell{1};
        end
        
        for i = 2:n
            out(i) = str2num(strcat(stringCell{i},'/',lcd));
        end
    end
end