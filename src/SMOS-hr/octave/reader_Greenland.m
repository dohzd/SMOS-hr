
fmt='C';
##xy2hv=true;
##diary example.log
##fprintf('Reading L1C + XY2HV %s ',fL1C);
##t0=tic;
##[l1cdata, l1cSNP, l1chdr, ~, ~, Asc]=read_L1cOP_V2(fL1C,fmt,xy2hv);
##fprintf('%d DGGs\n',size(l1cdata{1},1));
##toc0=toc(t0);
##fprintf('==> Elapsed Time %s\n',datestr(toc0/3600/24,'MM'':SS"'));
##save([fL1C '.mat'],'l1cdata','l1cSNP','l1chdr','Asc');

%% Quicker version; the conversion from XY to HV takes quite a lot of time.
%% When it is activated internally to read_L1cOP_V2() it is made for all DGGs of the products, including those not of interest.
%% In this case that focuses to Antarctica only; read_L1cOP_V2() without conversion will be quicker, then limit the L1C to DGG having Lat < -60°
%% and then apply XY2HV() externally on this limited L1C to complement the DGGs data with HV.

xy2hv=false;
%% Read directory as input argument
%%files = argv ()
%%files=glob("/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/L1C/2020/03/*/SM*300_1")
arg_list = argv ()
files = glob(strcat(arg_list{1},'/SM_*_1'))

for i=1:numel(files)
  fL1C = files{i}
  if ~exist([fL1C "_Greenland.mat"],"file")
    fprintf('\nReading L1C %s ',fL1C);
    t00=tic;
    t0=t00;
    %% Read the L1C as it is without calling internally XY2HV so with no additional HV data
    [l1cdata, l1cSNP, l1chdr, ~, ~, Asc]=read_L1cOP_V2(fL1C,fmt,xy2hv);
    fprintf('%d DGGs\n',size(l1cdata{1},1));
    fprintf('==> Elapsed Time %s\n',datestr(toc(t0)/3600/24,'MM'':SS"'));
    
    %% Then decrease the size of the l1cdata to DGGs falling into a Greenland-centered rectangle
    ixA=l1cdata{1}(:,2)>57 & l1cdata{1}(:,3)<-10 & l1cdata{1}(:,3)>-75;
    l1cdata={l1cdata{1}(ixA,:) l1cdata{2}(ixA)};
    
    %% Then call L1CXY2HV on this shrinked l1cdata to complement l1c BT profiles with additional Earth BTH,BTV,ST3,ST4 and their radiometric uncertainty.
    t0=tic;
    fprintf('XY2HV on %d selected DGGs having lat >57°\n',size(l1cdata{1},1));
    
    if size(l1cdata{1},1)>0
        try 
            l1cdata=L1CXY2HV(l1cdata,l1cSNP);
            
            %% Remove variables that are not used in rSIR image reconstruction
            %% keep: flags (1), incidence (5), azimuth (6), snapshot_ID (9), semi_major_axis (10), semi_minor_axis (11), TB_H (12), TB_V (13), TB_H_acc (16), TB_V_acc (17)
            indexes=[1,5,6,9,10,11,12,13,16,17]
    	    for dgg=1:size(l1cdata{2})
            	R=l1cdata{2}{dgg};
            	l1cdata{2}{dgg}=NaN(size(R,1),size(indexes,2));
                l1cdata{2}{dgg}=R(:,indexes);
            end

            fprintf('==> Elapsed Time %s\n',datestr(toc(t0)/3600/24,'MM'':SS"'));
            fprintf('==> Total Read L1C + select + L1CXY2HV on select %s\n',datestr(toc(t00)/3600/24,'MM'':SS"'));
            save('-v7',[fL1C '_Greenland.mat'],'l1cdata','l1cSNP','l1chdr','Asc');
        catch
            printf("Last error: %s\n",lasterr)
        end
    end

    diary off
  end
endfor


##cd "~/Documents/PZ/programs/L1c_reader"
