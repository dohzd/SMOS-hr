function l1c=L1CXY2HV(l1c,SNP)
% function l1c=L1CXY2HV(l1c)
%
% Input:
%   l1c: cell array; typical l1c read from read_L1cOP in cell array mode
%       l1c{1}: double matrix _nx5:l1c{1}(k,:)=[DGGID Lat Lon Alt Mask]
%       l1c{2} if Dual Pol: cell array _nx1 of matrices _p(_n) x 10 l1c{2}{k}(q,:)=[Fl BT RA Th Az Fa Ge Si A1 A2]
%       l1c{2} if Full Pol: cell array _nx1 of matricse _p(_n) x 11 l1c{2}{k}(q,:)=[Fl ReBT ImBT RA Th Az Fa Ge Si A1 A2]
%   SNP: double matrix _q x 31: snapshot information as read by read_L1cOP
%       SNP(k,:) = Time(1,3) SnapId OBET Xpos Ypos Zpos Xvel YVel Zvel VecSource Q0 Q1 Q2 Q3
%                  TEC GeoMagF GeoMagD GeoMagI SunRA SunDEC SunBT Acc RadAcc1 RadAcc2 X_Band
%                  SoftwareError InsError ADFError CalibFlag
% Output:
%   l1c{2} matrices complemented (cols):
%   For Dual Pol with BTX,BTY converted to Earth reference frame BTH,BTV,RAH,RAV
%   For Full Pol with BTX,BTY,ReBTXY,ImBTXY converted to Earth reference frame BTH,BTV,ST3,ST4,RAH,RAV,RA3,RA4
SNPID0=SNP(1,4)-1;
ts(SNP(:,4)-SNPID0)=SNP(:,1:3)*[24*3600;1;1e-6];
ts=ts';
if size(l1c{2}{1},2)==10
    %% Dual Pol
else
    %% Full Pol
    for igp=1:size(l1c{2})
        R=l1c{2}{igp};
        if size(R,1)<=6
            l1c{2}{igp}(:,end+1:end+1+7)=NaN;
        else
            vts=ts(R(:,9)-SNPID0);
            vpol=bitand(R(:,1),3);
            flgx=vpol==0; % if sum(flgx) < 2 ; return;end
            flgy=vpol==1; % if sum(flgy) < 2 ; return;end
            flgxy=vpol>=2; % if sum(flgxy) < 2 ; return;end
            if sum(flgx) >=2 && sum(flgy)>=2 && sum(flgxy)>=2
                RTB=NaN(size(R,1),4+3);
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
                RHV=[reshape(InvMR4*reshape(RTB(:,1:4)'.*repmat([1;1;2;-2],1,size(RTB,1)),4*size(RTB,1),1),4,size(RTB,1))', ...
                    reshape(sqrt(InvMR4.^2*reshape(RTB(:,[5:end end])'.^2.*repmat([1;1;4;4],1,size(RTB,1)),4*size(RTB,1),1)),4,size(RTB,1))'];
                l1c{2}{igp}=NaN(size(R,1),size(R,2)+size(RHV,2));
                l1c{2}{igp}(:,1:size(R,2))=R;
                l1c{2}{igp}(idxok,size(R,2)+1:end)=RHV;
            else
                l1c{2}{igp}(:,end+1:end+1+7)=NaN;
            end
        end
    end
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


