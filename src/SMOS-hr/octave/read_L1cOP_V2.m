function [T, SNP, hdr, att, polM, Asc]=read_L1cOP_V2(fL1C,fmt,XY2HV,polM)
% function [T, SNP, hdr, att, polM, Asc]=read_L1cOP(fL1C [,fmt='C', [ XY2HV=false, [ polM=[] ] ] ])
%
% Input:
%   fl1C: string; L1C directory (split DBL/HDR EE)
%   fmt: OPTIONAL char; 'A' T is given as a matrix _n x _p, 'C' T is given as a cell arrays _n x 2
%   XY2HV: OPTIONAL: provide the L1C in Earth TBH/TBV IMPLEMENTED only in fmt=='C' and Full Polarization
%   polM: OPTIONAL string: either 'D' or 'F'; if provided force le
%                        interpreation of the file structure irrespective of its file name. If not
%                        polarization mode is extracted from the filename position 16
%
% Output:
%   For fmt=='A'
%   T: double matrix _n x _p: L1C data block content:
%       Dual Pol DPGS DBL/HDR Format
%           GP_id Lat Lon Alt Mask Count Fl BT RA Th Az Fa Ge Si A1 A2
%           1     2   3   4   5     6    7  8  9  10 11 12 13 14 15 16
%       Full Pol DPGS DBL/HDR Format
%           GP_id Lat Lon Alt Mask Count Fl REBT IMBT RA Th Az Fa Ge Si A1 A2
%           1     2   3   4   5     6    7  8     9   10 11 12 13 14 15 16 17
%
%   For fmt=='C'
%   T: 1 x 2 cell array with L1C data block content in single precision
%      T{1} is _n x 5 matrix = GP_id Lat Lon Alt Mask
%                                1    2   3   4   5
%      T{2} is the following _n x 1 cell array of the following _p x BT_Count L1C record matrix
%           _p = 10 for Dual Pol
%           R(1:10,:) = Fl BT RA Th Az Fa Ge Si A1 A2
%                       1  2  3  4  5  6  7  8  9  10
%           _p = 11 for Full Pol
%           R(1:11,:) = Fl REBT IMBT RA Th Az Fa Ge Si A1 A2
%                       1   2    3   4  5  6  7  8  9  10 11
%
%   if XY2HV is given and not null, then the projected earth surface BTs and their added radiometric accuracy is added right to the R columns
%   XY2HV cand be done only in SMOS Full Polarization mode (which is the normal today conf).
%   T{2} becomes:
%           _p = 19 for Full Pol
%           R(1:19,:) = Fl REBT IMBT RA Th Az Fa Ge Si A1 A2 TBH, TBV, ST3, ST4, RAH, RAV, RA3, RA4
%                       1   2    3   4  5  6  7  8  9  10 11  12   13   14   15   16   17   18   19
%
%   For Dual Polarization for alpha=Geo+Faraday == 45° [90] the MR2 matrix is singular and no inverse exist.
%   Numerically the inversion of MR2(alpha) starts to become unstable in a neighborhood of few degrees around 45°
%   An experimental XY2HV is provided in this version
%   T{2} becomes:
%           _p = 15 for Dual Pol
%           R(1:19,:) = Fl BT RA Th Az Fa Ge Si A1 A2 TBH, TBV, RAH, RAV, Cond
%                       1  2  3  4  5  6  7  8  9  10 11   12   13   14   15
%   Cond is a number close to 1 when Inv(MR2(alpha)) is well condionned, high values indicate that the TBH/TBV were obtained
%   with a numerically close to singular MR2(alpha) maxtrix.
%   It is suggested to filter out these BT when cond is higher than 5 to 10.
%
%   Less DGGs is expected since interpolations in the antenna reference frame are necessary and can be done only with a sufficient number of BTs.
%   The list of of final BTs per DGG is also truncated, the two first and two last BTs are lost due to the interpolation.
%   The REBT IMBT RA associated with a TBH, TBV, ST3, ST4, RAH, RAV, RA3, RA4 is the "real" antenna measurement leading to this Earth BTs vector.
%   while the thre others (no reported) were obtained by linear interpolation
%
%   SNP: double matrix _q x 31
%       SNP(k,:) = Time(1,3) SnapId OBET Xpos Ypos Zpos Xvel YVel Zvel VecSource Q0 Q1 Q2 Q3
%                  1 2 3     4      5    6    7    8    9    10   11   12        13 14 15 16
%                  TEC GeoMagF GeoMagD GeoMagI SunRA SunDEC SunBT Acc RadAcc1 RadAcc2 X_Band
%                  17  18      19      20      21    22     23    24  25      26      27
%                  SoftwareError InsError ADFError CalibFlag
%                  28            29       30       31
%   hdr: structure; L1C product xml header content (read_xml)
%   att: structure; L1C product xml header descriptor content (read_xml_att)
%   polM: char; 'D' if the fL1C is a dual pol product or 'F' if it's full pol one
%
%   Version 2017023: XY2HV added for DP + fix bug shif in stream readi DP data
%
if ~exist('fmt','var') || isempty(fmt)
    fmt='C';
end

if ~exist('XY2HV','var') || isempty(XY2HV)
    XY2HV=false;
end

if ~exist(fL1C,'dir')
    error('read_L1cOP: %s is not existent or not a directory',fL1C);
end
[~,fl1c]=fileparts(fL1C);
if ~exist('polM','var') || isempty(polM)
    polM=fl1c(16);
elseif polM~=fl1c(16)
    warning('read_L1cOP: polM=%s parameter disagrees with filename polarization=%s',polM,fl1c(16));
end
if polM=='D'
    isDP=true;
elseif polM=='F'
    isDP=false;
else
    error('read_L1cOP: polM is not ''D'' or ''F'' - aborting');
end
isgzL1C=false;
## (A supprimer) : code pour décmopresser les fichiers si nécessaire
##[v,~]=system(['ls ' fL1C '/*.gz']);
##if v==0
##    fprintf('====> Gzipped L1C unarchiving to /tmp/Phil/%s\n',fl1c);
##    [v,m]=system(['[ ! -d /tmp/Phil ] && mkdir /tmp/Phil ; cp -Lfr ' fL1C ' /tmp/Phil/. ; gunzip -fr /tmp/Phil/' fl1c ]);
##    if v==1
##        error('======> Error during cp or gunzip: %s',m);
##    end
##    isgzL1C=true;
##    fL1C=['/tmp/Phil/' fl1c];
##end

%% Warning on zero BTR length
WonceV6eq0=true;
%[v l1c_hdr]=system(['ls ' escape_special_char(fL1C) filesep '*SCL[DF]1C*.HDR']);l1c_hdr=l1c_hdr(1:end-1);
%[v l1c_dbl]=system(['ls ' escape_special_char(fL1C) filesep '*SCL[DF]1C*.DBL']);l1c_dbl=l1c_dbl(1:end-1);
[v l1c_hdr]=system(['ls ' escape_special_char(fL1C) '/*.HDR']);l1c_hdr=l1c_hdr(1:end-1);
[v l1c_dbl]=system(['ls ' escape_special_char(fL1C) '/*.DBL']);l1c_dbl=l1c_dbl(1:end-1);

% hdr=xml_load(l1c_hdr,'off');
[hdr,att]=read_product_header_V2(fL1C);

%% Accomodate for 770 dictionary chage Radiometric_Accuracy_Scale => Radiometric_Resolution_Scale
RAname=fieldnames(hdr.Variable_Header.Specific_Product_Header);
RAname=RAname{cellfun(@(x)~isempty(x),strfind(RAname,'Radiometric'))};

v344h=str2num(hdr.Fixed_Header.Source.Creator_Version) >= 344;
Asc=hdr.Variable_Header.Specific_Product_Header.Main_Info.Time_Info.Ascending_Flag;
% att=read_xml_att(l1c_hdr);
lDS_Name=arrayfun(@(x)(x.Data_Set.DS_Name),hdr.Variable_Header.Specific_Product_Header.List_of_Data_Sets,'UniformOutput',0);
iSL=strmatch('Swath_Snapshot_List',lDS_Name);
if isempty(iSL) ; iSL=strmatch('SNAPSHOT_LIST',lDS_Name);end
if isDP
    iTSD=strmatch('Temp_Swath_Dual',lDS_Name);
    if isempty(iTSD) ; iTSD=strmatch('TEMP_SWATH_DUAL',lDS_Name);end
else
    iTSD=strmatch('Temp_Swath_Full',lDS_Name);
    if isempty(iTSD) ; iTSD=strmatch('TEMP_SWATH_FULL',lDS_Name);end
end
if isempty(iTSD)
    error('lit_L1cOP: polM=%s disagree whith product header dataset name - aborting',polM);
end
NumMDR=str2double(hdr.Variable_Header.Specific_Product_Header.List_of_Data_Sets(iTSD).Data_Set.Num_DSR);
SizMDR=str2double(hdr.Variable_Header.Specific_Product_Header.List_of_Data_Sets(iTSD).Data_Set.DSR_Size);
NumSNP=str2double(hdr.Variable_Header.Specific_Product_Header.List_of_Data_Sets(iSL).Data_Set.Num_DSR);
SizSNP=str2double(hdr.Variable_Header.Specific_Product_Header.List_of_Data_Sets(iSL).Data_Set.DSR_Size);

f=fopen(l1c_dbl,'rb','ieee-le');
Num_Snapshots  =fread(f,1,'uint32');

tempSNP = fread(f,[SizSNP,NumSNP],'*uint8')';
if isempty(tempSNP)
     T=[];
     fclose(f);
     return
end

sch_ver=str2double(hdr.Variable_Header.Specific_Product_Header.Main_Info.Datablock_Schema(24:27));
if sch_ver<=400
    SNP = zeros(NumSNP,31);
    SNP(:,    1) =         typecast(reshape(tempSNP(:,  1: 4 )',1,NumSNP*4  ),  'int32')            ; % Day
    SNP(:, 2: 3) = reshape(typecast(reshape(tempSNP(:,  5:12 )',1,NumSNP*4*2), 'uint32'),2,NumSNP)' ; % Second, Microsecond
    SNP(:,    4) =         typecast(reshape(tempSNP(:, 13:16 )',1,NumSNP*2*2), 'uint32')            ; % SnapId
    SNP(:,    5) =         typecast(reshape(tempSNP(:, 17:24 )',1,NumSNP*8  ),  'int64')            ; % OBET
    SNP(:, 6:11) = reshape(typecast(reshape(tempSNP(:, 25:72 )',1,NumSNP*8*6), 'double'),6,NumSNP)' ; % Xpos, Ypos, Zpos, Xvel, Yvel, Zvel
    SNP(:,   12) =                                tempSNP(:,    73 )                                       ; % VecSource
    SNP(:,13:20) = reshape(typecast(reshape(tempSNP(:, 74:137)',1,NumSNP*8*8), 'double'),8,NumSNP)' ; % Q0, Q1, Q2, Q3, TEC, GeomagF, GeomagD, GeomagI
    SNP(:,21:26) = reshape(typecast(reshape(tempSNP(:,138:161)',1,NumSNP*4*6), 'single'),6,NumSNP)' ; % SunRA, SunDec, SunBT, Acc, RadAcc1, Radacc2
    SNP(:,27:31) =                                tempSNP(:,162:166)                                       ; % X_Band, SoftwareError,InsError,ADFError,CalibFlag
else %% Flags 1 byte inserted after OBET for recent versions
    SNP = zeros(NumSNP,32);
    SNP(:,    1) =         typecast(reshape(tempSNP(:,  1: 4 )',1,NumSNP*4  ),  'int32')            ; % Day
    SNP(:, 2: 3) = reshape(typecast(reshape(tempSNP(:,  5:12 )',1,NumSNP*4*2), 'uint32'),2,NumSNP)' ; % Second, Microsecond
    SNP(:,    4) =         typecast(reshape(tempSNP(:, 13:16 )',1,NumSNP*2*2), 'uint32')            ; % SnapId
    SNP(:,    5) =         typecast(reshape(tempSNP(:, 17:24 )',1,NumSNP*8  ),  'int64')            ; % OBET
    SNP(:,    6) =               typecast(reshape(tempSNP(:, 25)',1,NumSNP  ),  'uint8')            ; % Flags
    SNP(:, 7:12) = reshape(typecast(reshape(tempSNP(:, 26:73 )',1,NumSNP*8*6), 'double'),6,NumSNP)' ; % Xpos, Ypos, Zpos, Xvel, Yvel, Zvel
    SNP(:,   13) =                                tempSNP(:,    74 )                                       ; % VecSource
    SNP(:,14:21) = reshape(typecast(reshape(tempSNP(:, 75:138)',1,NumSNP*8*8), 'double'),8,NumSNP)' ; % Q0, Q1, Q2, Q3, TEC, GeomagF, GeomagD, GeomagI
    SNP(:,22:27) = reshape(typecast(reshape(tempSNP(:,139:162)',1,NumSNP*4*6), 'single'),6,NumSNP)' ; % SunRA, SunDec, SunBT, Acc, RadAcc1, Radacc2
    SNP(:,28:32) =                                tempSNP(:,163:167)                                       ; % X_Band, SoftwareError,InsError,ADFError,CalibFlag
end

if XY2HV
    SNPID0=SNP(1,4)-1;
    ts(SNP(:,4)-SNPID0)=SNP(:,1:3)*[24*3600;1;1e-6];
    ts=ts';
end
Num_Grid_Points=fread(f,1,'uint32');

clear tempSNP

tempGP = fread(f,'*uint8');
if isDP %% Dual Pol
    %T=zeros(Num_Grid_Points,1700);
    switch fmt
        case 'A'
            T=zeros(Num_Grid_Points,2805);
            off=0;

            off=0;

            for gp=1:Num_Grid_Points

                T(gp,  1)     =         typecast(reshape(tempGP(off+( 1: 4 ))',1,4  ),  'uint32')       ; off = off+ 4; % GPid
                T(gp, 2: 4)   = reshape(typecast(reshape(tempGP(off+( 1:12 ))',1,4*3),  'single'),3,1)' ; off = off+12; % Lat, Lon, Alt
                T(gp,    5)   =                                tempGP(off+( 1    ))                            ; off = off+ 1; % Mask
                T(gp,    6)   =         typecast(reshape(tempGP(off+( 1: 2 ))',1,2  ),  'uint16')       ; off = off+ 2; % BT_Count

                doff = (1:T(gp,6))*24-24 + off ;
                sn   = (1:T(gp,6)) ;
                doff2= reshape((doff(ones(2,1),:)'+[ones(T(gp,6),1) ones(T(gp,6),1)*2])',1,2*T(gp,6));
                doff4= reshape((doff(ones(4,1),:)'+[ones(T(gp,6),1) ones(T(gp,6),1)*2 ones(T(gp,6),1)*3 ones(T(gp,6),1)*4])',1,4*T(gp,6));

                T(gp, 7+(sn-1)*10) = typecast(reshape(tempGP(doff2 + 0)',1,2 * T(gp,6)  ),  'uint16')      ;  % Flags
                T(gp, 8+(sn-1)*10) = typecast(reshape(tempGP(doff4 + 2)',1,4 * T(gp,6)  ),  'single')      ;  % BT
                T(gp, 9+(sn-1)*10) = typecast(reshape(tempGP(doff2 + 6)',1,2 * T(gp,6)  ),  'uint16')      ;  % RA
                T(gp,10+(sn-1)*10) = typecast(reshape(tempGP(doff2 + 8)',1,2 * T(gp,6)  ),  'uint16')      ;  % Th
                T(gp,11+(sn-1)*10) = typecast(reshape(tempGP(doff2 +10)',1,2 * T(gp,6)  ),  'uint16')      ;  % Az
                T(gp,12+(sn-1)*10) = typecast(reshape(tempGP(doff2 +12)',1,2 * T(gp,6)  ),  'uint16')      ;  % Fa
                T(gp,13+(sn-1)*10) = typecast(reshape(tempGP(doff2 +14)',1,2 * T(gp,6)  ),  'uint16')      ;  % Ge
                T(gp,14+(sn-1)*10) = typecast(reshape(tempGP(doff4 +16)',1,4 * T(gp,6)  ),  'uint32')      ;  % Si
                T(gp,15+(sn-1)*10) = typecast(reshape(tempGP(doff2 +20)',1,2 * T(gp,6)  ),  'uint16')      ;  % A1
                T(gp,16+(sn-1)*10) = typecast(reshape(tempGP(doff2 +22)',1,2 * T(gp,6)  ),  'uint16')      ;  % A2

                off = off + T(gp,6) * 24 ;

            end
            mx=max(T(:,6));
            T=T(:,1:mx*10+6);

            T(:, 9:10:end) = T(:, 9:10:end)* str2double(hdr.Variable_Header.Specific_Product_Header.(RAname))/(2^16);
            T(:,10:10:end) = T(:,10:10:end)* 90/(2^16);
            T(:,11:10:end) = T(:,11:10:end)*360/(2^16);
            T(:,12:10:end) = T(:,12:10:end)*360/(2^16);
            T(:,13:10:end) = T(:,13:10:end)*360/(2^16);
            T(:,15:10:end) = T(:,15:10:end)* str2double(hdr.Variable_Header.Specific_Product_Header.Pixel_Footprint_Scale)/(2^16);
            T(:,16:10:end) = T(:,16:10:end)* str2double(hdr.Variable_Header.Specific_Product_Header.Pixel_Footprint_Scale)/(2^16);
        case 'C'
            T={[],[]};
            T{1}=NaN(Num_Grid_Points,5);
            T{2}=cell(Num_Grid_Points,1);
            off=0;
            V=zeros(1,5);
            C1=str2double(hdr.Variable_Header.Specific_Product_Header.(RAname))/(2^16);
            C2=90/(2^16);
            C345=360/(2^16);
            C6=str2double(hdr.Variable_Header.Specific_Product_Header.Pixel_Footprint_Scale)/(2^16);
            C7=str2double(hdr.Variable_Header.Specific_Product_Header.Pixel_Footprint_Scale)/(2^16);
            nb_dgg_out=0;
            igp=1;
            NS=10;
            NSC=fix(Num_Grid_Points/NS);
            for gp=1:Num_Grid_Points
%                 if mod(gp-1,NSC)==0
%                     fprintf('%02d%%...',fix((gp-1)/NSC)*NS);
%                 end
                V6     = double(typecast(reshape(tempGP(off+17+( 1: 2 ))',1,2  ),  'uint16')); % BT_Count
%                 if XY2HV && V6 <= 6
%                     %% Too few TBs
%                     nb_dgg_out=nb_dgg_out+1;
%                     off = off + 19 + V6 * 24 ; % skip the whole follwoing off + 4 + 12 +1 +2 = 19 + BTCOUNT (V6) * 24 bytes
%                     continue;
%                 end

                V(1)   =         typecast(reshape(tempGP(off+( 1: 4 ))',1,4  ),  'uint32')       ; off = off+ 4; % GPid
                V(2:4) = reshape(typecast(reshape(tempGP(off+( 1:12 ))',1,4*3),  'single'),3,1)' ; off = off+12; % Lat, Lon, Alt
                V(5)   =                                tempGP(off+( 1    ))                            ; off = off+ 1; % Mask
                off = off+ 2; % V6 to skip alread read

                if V6 > 255
                    error('read_L1cOP:Corruted product: TB counter=%d too high',V6);
                end
                doff = (1:V6)*24-24 + off ;
                doff2= reshape((doff(ones(2,1),:)'+[ones(V6,1) ones(V6,1)*2])',1,2*V6);
                doff4= reshape((doff(ones(4,1),:)'+[ones(V6,1) ones(V6,1)*2 ones(V6,1)*3 ones(V6,1)*4])',1,4*V6);
                R=zeros(V6,10);
                R(:, 1) =        typecast(reshape(tempGP(doff2 + 0)',1,2 * V6  ),  'uint16')            ;  % Flags
                R(:, 2) =        typecast(reshape(tempGP(doff4 + 2)',1,4 * V6  ),  'single')            ;  % BT
                R(:, 3) = double(typecast(reshape(tempGP(doff2 + 6)',1,2 * V6  ),  'uint16'))*C1        ;  % RA
                R(:, 4) = double(typecast(reshape(tempGP(doff2 + 8)',1,2 * V6  ),  'uint16'))*C2        ;  % Th
                R(:, 5) = double(typecast(reshape(tempGP(doff2 +10)',1,2 * V6  ),  'uint16'))*C345      ;  % Az
                R(:, 6) = double(typecast(reshape(tempGP(doff2 +12)',1,2 * V6  ),  'uint16'))*C345      ;  % Fa
                R(:, 7) = double(typecast(reshape(tempGP(doff2 +14)',1,2 * V6  ),  'uint16'))*C345      ;  % Ge
                R(:, 8) =        typecast(reshape(tempGP(doff4 +16)',1,4 * V6  ),  'uint32')            ;  % Si
                R(:, 9) = double(typecast(reshape(tempGP(doff2 +20)',1,2 * V6  ),  'uint16'))*C6        ;  % A1
                R(:,10) = double(typecast(reshape(tempGP(doff2 +22)',1,2 * V6  ),  'uint16'))*C7        ;  % A2
                if XY2HV
                    if V6 <= 6
                        T{2}{igp}=[R,NaN(size(R,1),2+2)];
                    else
                        vts=ts(R(:,8)-SNPID0);
                        vpol=bitand(R(:,1),3);
                        flgx=vpol==0; % if sum(flgx) < 2 ; return;end
                        flgy=vpol==1; % if sum(flgy) < 2 ; return;end
                        if sum(flgx) >=2 && sum(flgy)>=2
                            RTB=NaN(V6,2+2);
                            %% too wide gaps in sequence is not handled
                            nflgx=~flgx; % nflgx = nflgx & [true diff(vts(nflgx))<3.8] ;
                            nflgy=~flgy; % nflgy = nflgy & [true diff(vts(nflgy))<3.8] ;

                            RTB(flgx,1)=R(flgx,2);
                            RTB(flgy,2)=R(flgy,2);

                            RTB(flgx,3)=R(flgx,3);
                            RTB(flgy,4)=R(flgy,3);

                            RTB(nflgx,1)=interp1(vts(flgx),R(flgx,2),vts(nflgx));
                            RTB(nflgy,2)=interp1(vts(flgy),R(flgy,2),vts(nflgy));

                            RTB(nflgx,3)=interp1(vts(flgx),R(flgx,3),vts(nflgx));
                            RTB(nflgy,4)=interp1(vts(flgy),R(flgy,3),vts(nflgy));

                            idxok=~any(isnan(RTB),2);
                            RTB=RTB(idxok,:);
                            Rok=R(idxok,:);

                            %% MR2 is a partial rotation matrix in the case of Dual Pol => inv(MR2(a)) != MR2(-a) we have to compute the
                            %% inverse matrix
                            %% a=Far+Geo
                            % InvMR2 = compute_MR2(-sum(R(:,6:7),2)) ;
                            InvMR2 = mat2cell(compute_InvMR2(sum(Rok(:,6:7),2)),2*ones(1,size(Rok,1)));
                            RInvMR2 = cellfun(@(x)cond(x),InvMR2);
                            InvMR2 = blkdiag(InvMR2{:});
                            %                        T{2}{igp}=single([ R(:,1), ...
                            %                                          reshape(InvMR2*reshape(RTB(:,1:4)'.*repmat([1;1;2;-2],1,size(RTB,1)),4*size(RTB,1),1),4,size(RTB,1))', ...
                            %                                          reshape(sqrt(InvMR2.^2*reshape(RTB(:,[5:end end])'.^2.*repmat([1;1;4;4],1,size(RTB,1)),4*size(RTB,1),1)),4,size(RTB,1))'; ...
                            %                                          R(:,5:end)]);
                            RHV=[reshape(InvMR2*reshape(RTB(:,1:2)'.*repmat([1;1],1,size(RTB,1)),2*size(RTB,1),1),2,size(RTB,1))', ...
                                reshape(sqrt(InvMR2.^2*reshape(RTB(:,3:end)'.^2.*repmat([1;1],1,size(RTB,1)),2*size(RTB,1),1)),2,size(RTB,1))', ...
                                RInvMR2];
                            T{2}{igp}=NaN(size(R,1),size(R,2)+size(RHV,2));
                            T{2}{igp}(:,1:size(R,2))=R;
                            T{2}{igp}(idxok,size(R,2)+1:end)=RHV;

                            %                         T{2}{igp}=[ R, ...
                            %                                     reshape(InvMR2*reshape(RTB(:,1:2)'.*repmat([1;1],1,size(RTB,1)),2*size(RTB,1),1),2,size(RTB,1))', ...
                            %                                     reshape(sqrt(InvMR2.^2*reshape(RTB(:,3:end)'.^2.*repmat([1;1],1,size(RTB,1)),2*size(RTB,1),1)),2,size(RTB,1))', ...
                            %                                     RInvMR2
                            %                                   ];
                            % Not that clesr the 2^2ù or only 2*
                            %                         T{2}{igp}=[ R, ...
                            %                                     reshape(InvMR2*reshape(RTB(:,1:4)'.*repmat([1;1;2;-2],1,size(RTB,1)),4*size(RTB,1),1),4,size(RTB,1))', ...
                            %                                     reshape(sqrt(InvMR2.^2*reshape(RTB(:,[5:end end])'.^2.*repmat([1;1;2;2],1,size(RTB,1)),4*size(RTB,1),1)),4,size(RTB,1))'; ...
                            %                                   ];


                        else
                            T{2}{igp}=[R,NaN(size(R,1),2+2)];
                        end
                    end
                else
                    T{2}{igp}=R;
                end
                T{1}(igp,:)=V;
                off = off + V6 * 24 ;
                igp=igp+1;
            end
%             fprintf('100%%\n');
    end

    if XY2HV
        T{1}=T{1}(1:igp-1,:);
        T{2}=T{2}(1:igp-1);
    end

else %% Full Pol
    switch fmt
        case 'A'
            T=zeros(Num_Grid_Points,2805);
            off=0;
            igp=1;
            for gp=1:Num_Grid_Points

                T(gp,  1)     =         typecast(reshape(tempGP(off+( 1: 4 ))',1,4  ),  'uint32')       ; off = off+ 4; % GPid
                T(gp, 2: 4)   = reshape(typecast(reshape(tempGP(off+( 1:12 ))',1,4*3),  'single'),3,1)' ; off = off+12; % Lat, Lon, Alt
                T(gp,    5)   =                                tempGP(off+( 1    ))                            ; off = off+ 1; % Mask
                T(gp,    6)   =         typecast(reshape(tempGP(off+( 1: 2 ))',1,2  ),  'uint16')       ; off = off+ 2; % BT_Count
                if T(gp,6) > 255
                    error('read_L1cOP:Corruted product: TB counter=%d too high',T(gp,6));
               end
                doff = (1:T(gp,6))*28-28 + off ;
                sn   = (1:T(gp,6)) ;
                doff2= reshape((doff(ones(2,1),:)'+[ones(T(gp,6),1) ones(T(gp,6),1)*2])',1,2*T(gp,6));
                doff4= reshape((doff(ones(4,1),:)'+[ones(T(gp,6),1) ones(T(gp,6),1)*2 ones(T(gp,6),1)*3 ones(T(gp,6),1)*4])',1,4*T(gp,6));

                T(gp, 7+(sn-1)*11) = typecast(reshape(tempGP(doff2 + 0)',1,2 * T(gp,6)  ),  'uint16')      ;  % Flags
                T(gp, 8+(sn-1)*11) = typecast(reshape(tempGP(doff4 + 2)',1,4 * T(gp,6)  ),  'single')      ;  % ReBT
                T(gp, 9+(sn-1)*11) = typecast(reshape(tempGP(doff4 + 6)',1,4 * T(gp,6)  ),  'single')      ;  % ImBT
                T(gp,10+(sn-1)*11) = typecast(reshape(tempGP(doff2 +10)',1,2 * T(gp,6)  ),  'uint16')      ;  % RA
                T(gp,11+(sn-1)*11) = typecast(reshape(tempGP(doff2 +12)',1,2 * T(gp,6)  ),  'uint16')      ;  % Th
                T(gp,12+(sn-1)*11) = typecast(reshape(tempGP(doff2 +14)',1,2 * T(gp,6)  ),  'uint16')      ;  % Az
                T(gp,13+(sn-1)*11) = typecast(reshape(tempGP(doff2 +16)',1,2 * T(gp,6)  ),  'uint16')      ;  % Fa
                T(gp,14+(sn-1)*11) = typecast(reshape(tempGP(doff2 +18)',1,2 * T(gp,6)  ),  'uint16')      ;  % Ge
                T(gp,15+(sn-1)*11) = typecast(reshape(tempGP(doff4 +20)',1,4 * T(gp,6)  ),  'uint32')      ;  % Si
                T(gp,16+(sn-1)*11) = typecast(reshape(tempGP(doff2 +24)',1,2 * T(gp,6)  ),  'uint16')      ;  % A1
                T(gp,17+(sn-1)*11) = typecast(reshape(tempGP(doff2 +26)',1,2 * T(gp,6)  ),  'uint16')      ;  % A2

                off = off + T(gp,6) * 28 ;

            end
            mx=max(T(:,6));
            T=T(:,1:mx*11+6);
            T(:,10:11:end) = T(:,10:11:end)* str2double(hdr.Variable_Header.Specific_Product_Header.(RAname))/(2^16);
            T(:,11:11:end) = T(:,11:11:end)* 90/(2^16);
            T(:,12:11:end) = T(:,12:11:end)*360/(2^16);
            T(:,13:11:end) = T(:,13:11:end)*360/(2^16);
            T(:,14:11:end) = T(:,14:11:end)*360/(2^16);
            T(:,16:11:end) = T(:,16:11:end)* str2double(hdr.Variable_Header.Specific_Product_Header.Pixel_Footprint_Scale)/(2^16);
            T(:,17:11:end) = T(:,17:11:end)* str2double(hdr.Variable_Header.Specific_Product_Header.Pixel_Footprint_Scale)/(2^16);
            if XY2HV
                fprintf('In Array mode XY2HV is not yet implemented\n')
            end
        case 'C'
            T={[],[]};
            T{1}=NaN(Num_Grid_Points,5,'single');
            T{2}=cell(Num_Grid_Points,1);
            off=0;
            V=zeros(1,5);
            C1=str2double(hdr.Variable_Header.Specific_Product_Header.(RAname))/(2^16);
            C2=90/(2^16);
            C345=360/(2^16);
            C6=str2double(hdr.Variable_Header.Specific_Product_Header.Pixel_Footprint_Scale)/(2^16);
            C7=str2double(hdr.Variable_Header.Specific_Product_Header.Pixel_Footprint_Scale)/(2^16);
            nb_dgg_out=0;
            igp=1;
            NS=10;
            NSC=fix(Num_Grid_Points/NS);
            for gp=1:Num_Grid_Points
%                 if mod(gp-1,NSC)==0
%                     fprintf('%02d%%...',fix((gp-1)/NSC)*NS);
%                 end
                V6     = double(typecast(reshape(tempGP(off+17+( 1: 2 ))',1,2  ),  'uint16'))       ; % BT_Count off +4+12+1= 17
%                 if XY2HV && V6 <= 6
%                     %% Too few TBs
%                     nb_dgg_out=nb_dgg_out+1;
%                     off = off + 19 + V6 * 28 ; % skip the whole follwoing off + 4 + 12 +1 +2 = 19 + BTCOUNT (V6) * 28 bytes
%                     continue;
%                 end
                V(1)   =         typecast(reshape(tempGP(off+( 1: 4 ))',1,4  ),  'uint32')       ; off = off+ 4; % GPid
                V(2:4) = reshape(typecast(reshape(tempGP(off+( 1:12 ))',1,4*3),  'single'),3,1)' ; off = off+12; % Lat, Lon, Alt
                V(5)   =                                tempGP(off+( 1    ))                            ; off = off+ 1; % Mask
                off = off+ 2;  % V6     = double(typecast(reshape(tempGP(off+( 1: 2 ))',1,2  ),  'uint16'))       ; off = off+ 2; % BT_Count

                if V6 > 255
                    error('read_L1cOP:Corruted product: TB counter=%d too high',V6);
                end
                if V6==0
                    if WonceV6eq0
                        warning('read_L1cOP:zero length BTR detected');
                        WonceV6eq0=false;
                    end
                        nb_dgg_out=nb_dgg_out+1;
%                         off = off + 1;
                        continue;
                end
                doff = (1:V6)*28-28 + off ;
                doff2= reshape((doff(ones(2,1),:)'+[ones(V6,1) ones(V6,1)*2])',1,2*V6);
                doff4= reshape((doff(ones(4,1),:)'+[ones(V6,1) ones(V6,1)*2 ones(V6,1)*3 ones(V6,1)*4])',1,4*V6);
                R=zeros(V6,11);
                R(:,1 ) = typecast(reshape(tempGP(doff2 + 0)',1,2 * V6  ),  'uint16')                   ;  % Flags
                R(:,2 ) = typecast(reshape(tempGP(doff4 + 2)',1,4 * V6  ),  'single')                   ;  % ReBT
                R(:,3 ) = typecast(reshape(tempGP(doff4 + 6)',1,4 * V6  ),  'single')                   ;  % ImBT
                R(:,4 ) = double(typecast(reshape(tempGP(doff2 +10)',1,2 * V6  ),  'uint16'))*C1        ;  % RA
                R(:,5 ) = double(typecast(reshape(tempGP(doff2 +12)',1,2 * V6  ),  'uint16'))*C2        ;  % Th
                R(:,6 ) = double(typecast(reshape(tempGP(doff2 +14)',1,2 * V6  ),  'uint16'))*C345      ;  % Az
                R(:,7 ) = double(typecast(reshape(tempGP(doff2 +16)',1,2 * V6  ),  'uint16'))*C345      ;  % Fa
                R(:,8 ) = double(typecast(reshape(tempGP(doff2 +18)',1,2 * V6  ),  'uint16'))*C345      ;  % Ge
                R(:,9 ) = typecast(reshape(tempGP(doff4 +20)',1,4 * V6  ),  'uint32')                   ;  % Si
                R(:,10) = double(typecast(reshape(tempGP(doff2 +24)',1,2 * V6  ),  'uint16'))*C6        ;  % A1
                R(:,11) = double(typecast(reshape(tempGP(doff2 +26)',1,2 * V6  ),  'uint16'))*C7        ;  % A2
                if XY2HV
                    if V6<=6
                        T{2}{igp}=[R,NaN(size(R,1),4+4)];
                    else
                        vts=ts(R(:,9)-SNPID0);
                        vpol=bitand(R(:,1),3);
                        flgx=vpol==0; % if sum(flgx) < 2 ; return;end
                        flgy=vpol==1; % if sum(flgy) < 2 ; return;end
                        flgxy=vpol>=2; % if sum(flgxy) < 2 ; return;end
                        if sum(flgx) >=2 && sum(flgy)>=2 && sum(flgxy)>=2
                            RTB=NaN(V6,4+3);
                            %% too wide gaps in sequence is not handled
                            nflgx=~flgx; % nflgx = nflgx & [true diff(vts(nflgx))<3.8] ;
                            nflgy=~flgy; % nflgy = nflgy & [true diff(vts(nflgy))<3.8] ;
                            nflgxy=~flgxy; % nflgxy = nflgxy & [true diff(vts(nflgxy))<2.5] ;

                            RTB(flgx,1)=R(flgx,2);
                            RTB(flgy,2)=R(flgy,2);
                            RTB(flgxy,3)=R(flgxy,2);
                            RTB(flgxy,4)=R(flgxy,3);

                            RTB(flgx,5)=R(flgx,4);
                            RTB(flgy,6)=R(flgy,4);
                            RTB(flgxy,7)=R(flgxy,4);

                            RTB(nflgx,1)=interp1(vts(flgx),R(flgx,2),vts(nflgx));
                            RTB(nflgy,2)=interp1(vts(flgy),R(flgy,2),vts(nflgy));
                            RTB(nflgxy,3)=interp1(vts(flgxy),R(flgxy,2),vts(nflgxy));
                            RTB(nflgxy,4)=interp1(vts(flgxy),R(flgxy,3),vts(nflgxy));

                            RTB(nflgx,5)=interp1(vts(flgx),R(flgx,4),vts(nflgx));
                            RTB(nflgy,6)=interp1(vts(flgy),R(flgy,4),vts(nflgy));
                            RTB(nflgxy,7)=interp1(vts(flgxy),R(flgxy,4),vts(nflgxy));

                            idxok=~any(isnan(RTB),2);
                            RTB=RTB(idxok,:);
                            Rok=R(idxok,:);

                            %% MR4 as a true rotation matrix in the case of Full Pol => inv(MR4(a))=MR4(-a) then:
                            %% a=Far+Geo
                            % InvMR4 = compute_MR4(-sum(R(:,7:8),2)) ;
                            InvMR4 = mat2cell(compute_MR4(-sum(Rok(:,7:8),2)),4*ones(1,size(Rok,1)));
                            InvMR4 = blkdiag(InvMR4{:});
                            %                        T{2}{igp}=single([ R(:,1), ...
                            %                                          reshape(InvMR4*reshape(RTB(:,1:4)'.*repmat([1;1;2;-2],1,size(RTB,1)),4*size(RTB,1),1),4,size(RTB,1))', ...
                            %                                          reshape(sqrt(InvMR4.^2*reshape(RTB(:,[5:end end])'.^2.*repmat([1;1;4;4],1,size(RTB,1)),4*size(RTB,1),1)),4,size(RTB,1))'; ...
                            %                                          R(:,5:end)]);
                            RHV=[reshape(InvMR4*reshape(RTB(:,1:4)'.*repmat([1;1;2;-2],1,size(RTB,1)),4*size(RTB,1),1),4,size(RTB,1))', ...
                                reshape(sqrt(InvMR4.^2*reshape(RTB(:,[5:end end])'.^2.*repmat([1;1;4;4],1,size(RTB,1)),4*size(RTB,1),1)),4,size(RTB,1))'];
                            T{2}{igp}=NaN(size(R,1),size(R,2)+size(RHV,2));
                            T{2}{igp}(:,1:size(R,2))=R;
                            T{2}{igp}(idxok,size(R,2)+1:end)=RHV;

                            %                         T{2}{igp}=[ R, ...
                            %                                     reshape(InvMR4*reshape(RTB(:,1:4)'.*repmat([1;1;2;-2],1,size(RTB,1)),4*size(RTB,1),1),4,size(RTB,1))', ...
                            %                                     reshape(sqrt(InvMR4.^2*reshape(RTB(:,[5:end end])'.^2.*repmat([1;1;4;4],1,size(RTB,1)),4*size(RTB,1),1)),4,size(RTB,1))'; ...
                            %                                   ];
                            % Not that clesr the 2^2ù or only 2*
                            %                         T{2}{igp}=[ R, ...
                            %                                     reshape(InvMR4*reshape(RTB(:,1:4)'.*repmat([1;1;2;-2],1,size(RTB,1)),4*size(RTB,1),1),4,size(RTB,1))', ...
                            %                                     reshape(sqrt(InvMR4.^2*reshape(RTB(:,[5:end end])'.^2.*repmat([1;1;2;2],1,size(RTB,1)),4*size(RTB,1),1)),4,size(RTB,1))'; ...
                            %                                   ];


                        else
                            %                             nb_dgg_out=nb_dgg_out+1;
                            %                             off = off + V6 * 28 ;
                            %                             continue;
                            T{2}{igp}=[R,NaN(size(R,1),4+4)];
                        end
                    end
                else
                    T{2}{igp}=R;
                end
                T{1}(igp,:)=V;
                off = off + V6 * 28 ;
                igp=igp+1;
            end
%             fprintf('100%%\n');
    end

    if XY2HV
        T{1}=T{1}(1:igp-1,:);
        T{2}=T{2}(1:igp-1);
    end
end
fclose(f);
if isgzL1C
    [~,~]=system(['rm -f /tmp/Phil/' fl1c '/' fl1c '.DBL /tmp/Phil/' fl1c '/' fl1c '.HDR ; rmdir /tmp/Phil/' fl1c]);
end

end

function MR4 = compute_MR4(alpha)
% MR4
% [ cos(a)^2,  sin(a)^2, -sin(2*a)/2, 0]
% [ sin(a)^2,  cos(a)^2,  sin(2*a)/2, 0]
% [ sin(2*a), -sin(2*a),    cos(2*a), 0]
% [        0,         0,           0, 1]
    c=cosd(alpha);
    cc=c.*c;
    ss=1-cc;
    c2=cosd(2*alpha);
    s2=sind(2*alpha);
    s22=s2/2;
    MR4=zeros(length(alpha)*4,4);
    MR4(1:4:end,1)=cc;MR4(1:4:end,2)=ss; MR4(1:4:end,3)=-s22;
    MR4(2:4:end,1)=ss;MR4(2:4:end,2)=cc; MR4(2:4:end,3)=s22;
    MR4(3:4:end,1)=s2;MR4(3:4:end,2)=-s2;MR4(3:4:end,3)=c2;
    MR4(4:4:end,4)=1;
end
function IMR2 = compute_InvMR2(alpha)
% MR2
% [ cos(a)^2,  sin(a)^2 ]
% [ sin(a)^2,  cos(a)^2 ]
% Inv(MR2)
% [  cc/(cc^2 - ss^2), -ss/(cc^2 - ss^2)]
% [ -ss/(cc^2 - ss^2),  cc/(cc^2 - ss^2)]

    c=cosd(alpha);
    cc=c.*c;
    ss=1-cc;
    d=1./(cc.^2 - ss.^2);
    cc=cc.*d;
    ss=ss.*d;
    IMR2=zeros(length(alpha)*2,2);
    IMR2(1:2:end,1)=cc;IMR2(1:2:end,2)=-ss;
    IMR2(2:2:end,1)=-ss;IMR2(2:2:end,2)=cc;
end
function MR2 = compute_MR2(alpha)
% MR2
% [ cos(a)^2,  sin(a)^2 ]
% [ sin(a)^2,  cos(a)^2 ]
    c=cosd(alpha);
    cc=c.*c;
    ss=1-cc;
    MR2=zeros(length(alpha)*2,2);
    MR2(1:2:end,1)=cc;MR2(1:2:end,2)=ss;
    MR2(2:2:end,1)=ss;MR2(2:2:end,2)=cc;
end

function s=escape_special_char(s)
end
