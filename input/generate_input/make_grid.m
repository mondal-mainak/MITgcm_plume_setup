clear all
addpath(genpath('/scratch/mm10845/MITgcm_base/verification/MITgcm_10m_fjord/m_files/GSW-Matlab-master'));
addpath(genpath('/scratch/mm10845/MITgcm_base/utils/matlab'));
prec='real*8';
ieee='b';

makeGrid = 1;

%%reading grid data%%%%
filename='/scratch/mm10845/MITgcm_base/verification/MITgcm_5m_fjord/m_files3/grid_des';
T = readtable(filename,'Delimiter',' ');
[r,c] = size(T);
for ii = 1:r
eval([T.Var1{ii} '=' num2str(T.Var2(ii))])
end
%grid
run read_grid.m;
xCells = round(0.45*xdir);
if(botDepth2~=botDepth)
    disp('discrepancy in the vertical gridpts : check grid_des and data');
else
xspan=xdir-xCells;
%horizontal resolution
dx = ones(1,xspan) .* deltaX;
dy = ones(1,ydir) .* deltaX;

nx = length(dx);

%bottom depth
%botDepth = 650;

%r1 = 1; %vertical resolution
h1 = botDepth; 
r1=deltaZ;
%number of z cells
n1 = h1/r1;

dz = [r1*ones(n1,1)];

z = cumsum(dz);
zc = .5*(z(1:end-1)+z(2:end));
zc = [zc; botDepth];

nz = length(z);

%%

%minX = 10;
maxX = maX + deltaX; %100 m to 100km sponge layer

tot = maxX - deltaX;

%xCells = 40; 

c = 1; %scaling parameter
n = 1; %counter

ext_x = zeros(1,xCells);

while(nansum(ext_x) <= tot) %loop until sponge layer distance is reached
      
      for i = 1:xCells
            
            edx = (deltaX + c*i^2) - c; %dx extension
            ext_x(i) = round2even(edx); %only even dx values for cell
            
      end
      
      c = n*0.001;
      n = n+1;
      
end

% for i=1:nx-50
%     
%     dx(i) = 10;
%     
% end

dx = [fliplr(ext_x) dx];

h = zeros(xdir,ydir).*h1;

%%
%T/S initial conditions

%load('Rink_Mean_TS.mat');
load 'sal3.mat';
load 'temp3.mat';
temp=temp3(1:end-1)';
salt=sal3(1:end-1)';
d = 1:botDepth; %depth bins
% spc = 2;
% 
% it = t(1:end-1);
% is = s(1:end-1);
% 
% %fill nan at surface
% it(1) = it(2);
% is(1) = is(2);
% 
% %extend profile to bottom depth
% ind = find(isnan(it),1);
% it(ind:end) = it(ind-1);
% is(ind:end) = is(ind-1);
% 
 temp = interp1(d,temp,z); %interpolate to z cells
 salt = interp1(d,salt,z);
% 
% temp(1) = temp(2);
% salt(1) = salt(2);

sa = gsw_SA_from_SP(salt,z,-53.38,71.313);
ct = gsw_CT_from_t(sa,temp,z);
pt = gsw_pt_from_CT(sa,ct);
[N2, p_mid] = gsw_Nsquared(sa,ct,z,71);
rho = densmdjwf(sa,pt,0); %potential density
plott=0;
if (plott==1)
hFig1 = figure(1);
set(hFig1, 'Position', [1 1 1024 768])
set(gcf, 'color', [1 1 1]);

subplot(121);

plot(temp,-z,'LineWidth',2);
set(gca,'FontSize',22);
xlabel('T (\circC)','FontSize',28);
ylabel('Depth (m)','FontSize',28);

box on
grid on

subplot(122);

plot(salt,-z,'LineWidth',2);
set(gca,'FontSize',22);
xlabel('S (psu)','FontSize',28);

box on
grid on
end
%%
%extend domain

%reshape temp and salt
T = repmat(reshape(temp,[1 1 length(z)]),[length(dx),length(dy),1]); %temperature grid
S = repmat(reshape(salt,[1 1 length(z)]),[length(dx),length(dy),1]); %salinity grid

%west
tw = repmat(temp, [1 length(dy)])';
sw = repmat(salt, [1 length(dy)])';

%north
tn = repmat(temp, [1 length(dx)])';
sn = repmat(salt, [1 length(dx)])';

%east
te = tw * 0; 
te(:) = gsw_CT_freezing(0,botDepth);

se = tw * 0; 

%find the area of the conduit
cst_ind = ydir/2-chan_halfwidth;
cen_ind = ydir/2+chan_halfwidth;
conduit =(cst_ind: cen_ind);
dA=(cen_ind-cst_ind)*deltaX*cond_grid*deltaZ;

ue = (tw*0) - qsg/dA;  

days = 1;

tw_td = zeros(length(dy),length(z),days);
sw_td = zeros(length(dy),length(z),days);

tn_td = zeros(length(dx),length(z),days);
ts_td = zeros(length(dx),length(z),days);

te_td = zeros(length(dy),length(z),days);
se_td = zeros(length(dy),length(z),days);
ue_td = zeros(length(dy),length(z),days);

c = 1;
relaxTime = 4*2; %3 days
useRelax = 0;

for i = 1:(days) %forcing at 21600 s
      
      if (i <= relaxTime && useRelax)
            
            relax = i/relaxTime;
            
      else
            
            relax = 1;
            
      end
      
      te_td(1:length(dy),1:length(z),c) = te;
      se_td(1:length(dy),1:length(z),c) = se;
      ue_td(1:length(dy),1:length(z),c) = ue * relax;
      %ue_td(end/2-1:end/2+1,end,:)=-10;
%       ue_td(end/3-1:end/3+1,end,:)=-10;
       %ue_td(2*end/3-1:2*end/3+1,end,:)=-10;
      tw_td(1:length(dy),1:length(z),c) = tw;
      sw_td(1:length(dy),1:length(z),c) = sw;
      
      tn_td(1:length(dx),1:length(z),c) = tn;
      sn_td(1:length(dx),1:length(z),c) = sn;
      
      c = c + 1;
      
end
disp(ue_td(10,10));
%%
%save grid

if makeGrid
    
    delete *.bin %remove previous files
    
    disp('Saving grid...')

    %bathymetery
    fid = fopen('bathy.bin','w',ieee); fwrite(fid,-h,prec); fclose(fid);

    %dx,dy,dz
    fid = fopen('dx.bin','w',ieee); fwrite(fid,dx,prec); fclose(fid);
    fid = fopen('dy.bin','w',ieee); fwrite(fid,dy,prec); fclose(fid);
    fid = fopen('dz.bin','w',ieee); fwrite(fid,dz,prec); fclose(fid);

    %TS
    fid = fopen('temp.bin','w',ieee); fwrite(fid,T,prec); fclose(fid);
    fid = fopen('salt.bin','w',ieee); fwrite(fid,S,prec); fclose(fid);
  
    %east / west boundary conditions
    fid = fopen('OBET.bin','w',ieee); fwrite(fid,te_td,prec); fclose(fid);
    fid = fopen('OBES.bin','w',ieee); fwrite(fid,se_td,prec); fclose(fid);
    fid = fopen('OBEU.bin','w',ieee); fwrite(fid,ue_td,prec); fclose(fid);
    %fid = fopen('OBWT.bin','w',ieee); fwrite(fid,tw_td,prec); fclose(fid);
    %fid = fopen('OBWS.bin','w',ieee); fwrite(fid,sw_td,prec); fclose(fid);
    
    %north / south boundary conditionss
    %fid = fopen('OBNT.bin','w',ieee); fwrite(fid,tn,prec); fclose(fid);
    %fid = fopen('OBNS.bin','w',ieee); fwrite(fid,sn,prec); fclose(fid);
        
    disp('Done!');

end
end
%%%%%%%%%%%%%%%%

