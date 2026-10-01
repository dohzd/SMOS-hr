function [hdr att]=read_product_header_V2(Prep)
% [hdr att]=read_product_header_V2(Prep)
% Input:
%   Prep: string or cell array of strings: Product folder name or directly
%   an header content as a cell array of strings
%
% Output:
%   hdr: structure; the structred form of the xml header product
%   att: structure; the attribute related information (useful for writing
%   header)
if exist(Prep,'file')
    if isempty(strfind(Prep,'.EEF'))
        [p n]=fileparts(Prep);
        if ~isempty(p)
            p=[p '/'];
        end
        fhdr=[Prep '/' n '.HDR'];
    else
        fhdr=Prep;
    end
else
    fhdr=Prep;
end
% hdr=xml_load(fhdr,'off');
hdr=read_xml_V2(fhdr);
if nargout == 2
    att = read_xml_att_V2(fhdr);
end
end


