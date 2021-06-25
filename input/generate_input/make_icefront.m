close
%clear all;
addpath(genpath('/scratch/mm10845/MITgcm_base/verification/MITgcm_10m_fjord/m_files/GSW-Matlab-master'));
%dimensions of grid

%%reading grid data%%%%
filename='/scratch/mm10845/MITgcm_base/verification/MITgcm_5m_fjord/m_files3/grid_des';
T = readtable(filename,'Delimiter',' ');
[r,c] = size(T);
for ii = 1:r
eval([T.Var1{ii} '=' num2str(T.Var2(ii))])
end
nz = 325;
%deltaZ = 5.0;
run read_grid.m;
eos = 'jmd95z';
acc = 'real*8';

%nominal depth of model (meters)
H = -botDepth;

bathy = ones(xdir,ydir)*H;

fid = fopen('bathy.bin','w','b'); fwrite(fid,bathy,acc);fclose(fid);

dz = deltaZ*ones(1,nz);
zgp1 = [0,cumsum(dz)];
zc = .5*(zgp1(1:end-1)+zgp1(2:end));
zg = zgp1(1:end-1);
dz = diff(zgp1);
dx=5.;
%gravity
gravity = 9.81;
rhoConst = 1030;

%compute potential field underneath ice shelf
talpha = 2e-4;
sbeta  = 7.4e-4;
tref = -1.9*ones(nz,1);
t = tref;
sref = 34*ones(nz,1);
s = sref;
gravity = 9.81;
k = 1;
dzm = abs([zg(1)-zc(1) .5*diff(zc)]);
dzp = abs([.5*diff(zc) zc(end)-zg(end)]);
p = abs(zc)*gravity*rhoConst*1e-4;
dp = p;
kp = 0;

while std(dp) > 1e-13
  phiHydF(k) = 0;
  p0 = p;
  kp = kp+1
  for k = 1:nz
    switch eos
     case 'linear'
      drho = rhoConst*(1-talpha*(t(k)-tref(k))+sbeta*(s(k)-sref(k)))-rhoConst;
     case 'jmd95z'
      drho = densjmd95(s(k),t(k),p(k))-rhoConst;
     case 'mdjwf'
      drho = densmdjwf(s(k),t(k),p(k))-rhoConst;
     otherwise
      error(sprintf('unknown EOS: %s',eos))
    end
    phiHydC(k)   = phiHydF(k) + dzm(k)*gravity*drho/rhoConst;
    phiHydF(k+1) = phiHydC(k) + dzp(k)*gravity*drho/rhoConst;
  end
  switch eos
   case 'mdjwf'
    p = (gravity*rhoConst*abs(zc) + phiHydC*rhoConst)/gravity/rhoConst;
  end
  dp = p-p0;
end

icetopo = zeros(xdir,ydir);
runoffVel = icetopo.*0;
runoffRad = icetopo.*0;
plumeMask = icetopo.*0;
icetopo(end,:)= H; %30 m conduit width
buf=1.2;
geo=1;
span=25;
dh=H/span;
for i=1:span
    icetopo(end-i+1,:)= dh*(span-i+1);
    if(geo==2)
        conduit = [ydir/2-chan_halfwidth-round(i/buf):ydir/2+chan_halfwidth+round(i/buf)];
        icetopo(end-i,conduit)=H+cond_grid*deltaZ;
    end
    
end
iceetopo(xdir,conduit)=H+cond_grid*deltaZ;

load mr.mat;
[b,c]=size(mr);
deltaZ=abs(H)/c;
ind=round(dh/deltaZ);
tmp(ydir)=0;
for j=1:ydir
    for x=1:xdir
        iind=xdir+1-x;
        if(icetopo(iind,j)==0)
            tmp(j)=iind;
            break;
        end
    end
end
 
for j=1:ydir
    for i=tmp(j):xdir
        %hh=dh*(i-tmp(j)+1);
        %kind=round(hh/deltaZ)-1;
        icetopo(i,j)=icetopo(i,j)+mr(j,13*(xdir-tmp(j)-1))^2*20;
        if(icetopo(i,j)>0)
            icetopo(i,j)=0;
        end
            
            
    end
end
iceetopo(xdir-1:xdir,conduit)=H+cond_grid*deltaZ;
        
        
        
    
   
    
 %30 m conduit height 
%  for i=xdir-15:xdir-10
%     for j=conduit(1)-5:conduit(end)+5
%        icetopo(i,j)=0.3*(icetopo(i,j+1)+icetopo(i,j)+icetopo(i,j-1));
%     end
%  end
% conduit = [ydir/2-chan_halfwidth:ydir/2+chan_halfwidth];
% icetopo(end,conduit) = H+cond_grid*deltaZ;
fid=fopen('icetopo.bin','w','b'); fwrite(fid,icetopo,acc);fclose(fid);
%fid=fopen('pload.bin','w','b'); fwrite(fid,-icetopo,acc);fclose(fid);

phi0surf = zeros(xdir,ydir);
for ix=1:xdir
      for iy=1:ydir
            k=max(find(abs(zg)<abs(icetopo(ix,iy))));
            if isempty(k)
                  k=0;
            end
            if k>0
                  phi0surf(ix,iy) = phiHydF(k)*rhoConst;
            end
      end
end
fid=fopen('phi0surf.bin','w','b'); fwrite(fid,-icetopo,acc);fclose(fid);
%icefront files

%surf(icetopo);
%movefile *.bin ../run4/
surf(icetopo);

 frontdepth=zeros(xdir,ydir);
%frontdepth(end-1,:)=H+deltaZ;
for j=1:ydir
    for i=1:xdir-1
       
        frontdepth(i,j)=icetopo(i+1,j);
    end
end
load mr.mat;
frontdepth(end,:)=0;
frontlength=zeros(xdir,ydir);
k=1;
for j=1:ydir
    for i=tmp(j):xdir
      if(k < 325)  
        frontlength(i,j)=(deltaX-mr(j,k))/(deltaX)^2;
        k=k+13;
      end
    end
end

%frontlength(end-1,:)=1./deltaX;

heff=zeros(xdir,ydir);
for j=1:ydir
    for i=1:tmp(j)-1
        heff(i,j)=0.5;
    end
end

writebin('frontdepth.bin',frontdepth,1,acc);
writebin('frontlen.bin',frontlength,1,acc);  
writebin('heff.bin',frontdepth,1,acc);
