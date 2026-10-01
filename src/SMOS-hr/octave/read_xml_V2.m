function [xml att]=read_xml_V2(file)
if ischar(file)
    xmlstr=strsplit_V2(fileread(file));
else
    xmlstr=file;
end
n=1;
if strmatch('<?xml version',xmlstr{n})
    att.version=xmlstr{1};
    n=n+1;
end
att.root=xmlstr{n};

[xml att]=parse([],att,xmlstr(n+1:end-1));

end

function [xst xatt xmlstr]=parse(xst,xatt,xmlstr)

while ~isempty(xmlstr)
    [tag att data xmlstr brk]=next(xmlstr);
    [xst.(tag) xatt.(tag) xmlstr]=parse_next(tag,xst,xatt,xmlstr);
end
end

function [xst xatt xmlstr]=parse_next(stag,xst,xatt,xmlstr)
xst=[];
xatt=[];
if isempty(xmlstr)
    return
end
[tag att data xmlstr brk count]=next(xmlstr);
while ~strcmp(tag,['/' stag]) && ~isempty(xmlstr)
    if ~brk
        xst.(tag)=data;
%        if ~isempty(att)
            xatt.att.(tag)=att;
%        end
    else
        if count==1
            [xst.(tag) xatt.(tag) xmlstr]=parse_next(tag,[],[],xmlstr);
        else
            xatt.att.(tag)=att;
            for k=1:count
               [ctag catt cdata xmlstr brk]=next(xmlstr);
               [xst.(tag)(k).(ctag) xatt.(tag)(k).(ctag) xmlstr]=parse_next(ctag,[],[],xmlstr);
            end
            [tag att data xmlstr brk count]=next(xmlstr);
        end
    end

    [tag att data xmlstr brk count]=next(xmlstr);

end
end


function [tag att data xmlstr brk count]=next(xmlstr)
tag=[]; att=[]; data=[]; brk=0; count=1;
c=true;
ec=true;
while c
    str=xmlstr{1}; xmlstr=xmlstr(2:end);
    if isempty(xmlstr)
        return;
    end
    if ~ec
        if ~isempty(strfind(str,'-->'))
            ec=true;
        end
        continue
    end
    if ~isempty(strfind(str,'<!--'))
        c=true;
        if isempty(strfind(str,'-->'))
            ec=false;
        else
            ec=true;
        end
        continue
    end
    c=false;

    p=find(str=='>',1);
    n=find(str(2:p-1)==' ',1);
    if isempty(n)
        n=p-1;
    end
    tag=str(2:n);
    if tag(end)=='/' %% short cut field
        tag=tag(1:end-1);
        data='';
        return;
    end
    if n+2 < p-1
        att=str(n+2:p-1);
    else
        att=[''];
    end
    q=find(str(p+1:end)=='<',1);
    if ~isempty(q)
        data=str(p+1:p+q-1);
    elseif isempty(find(str(n+1:p)=='/',1))
        brk=1;
    end
    p=strmatch('count=',att);
    if ~isempty(p)
       count=str2num(att(p+7:p+5+find(att(p+7:end)=='"')));
    end
end
end
