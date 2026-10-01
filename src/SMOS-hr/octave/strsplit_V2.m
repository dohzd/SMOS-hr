function cstr = strsplit_V2(strarr, ch)
if ~exist('ch','var')
    ch=char(10);
end
idx=strfind(strarr,ch);
lch=length(ch);
n=1;
idx=[1-lch idx length(strarr)+1];
for k=1:length(idx)-1
    s=strtrim(strarr(idx(k)+lch:idx(k+1)-1));
    if ~isempty(s)
        cstr{n}=s;
        n=n+1;
    end
end
    
    