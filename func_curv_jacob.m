%Constructs curvilinear mesh and its Cartesian analogus. Calculates
%Jacobian 
% J=[dksi_dx dksi_dy;
%    deta_dx deta_dy];
% Ji=[dx_dksi dx_deta;
%     dy_dksi dy_deta];

% Input arguments:
% nx - number of nx grid points
% ny - number of ny grid points
% xmin, xmax - min and max values over OX
% ymin, ymax - min and max values over OY
% argument - string, for sin(argument), as function of x, that depend on xmin and xmax
% dxx - x spacing for Cartesian grid
% dyy - y spacing for Cartesian grid

function [xx,yy,ksi,eta,J] = func_curv_jacob(nx,ny,xmin,xmax,ymin,ymax,argument,dxx,dyy,curvature,showornot)
%npml=10;

%Set rang
ymid=(ymax+ymin)/2;

% dksi=xmax/nx;
% deta=ymax/ny;
% ksipmls=dksi*2*npml;
% etapmls=deta*2*npml;
% xmax=xmax+ksipmls;
% ymax=ymax+etapmls;

x=linspace(xmin,xmax,nx+1);
%argument of sin
%phi=2*pi*x/max(x);
% phi=-(1.25*pi*x/max(x)+0.25*pi);

phi=inline(argument,'x');


%Bottom part
y=curvature*ymax*sin(phi(x));
y=ymid+y;
%Top part, shifted by pi to be in phase
y1=-ymid-curvature*ymax*sin(phi(x)-pi);

%Phys. domain
xx=zeros(nx+1,ny+1);
yy=zeros(nx+1,ny+1);

%Comp. domain
ksi=zeros(nx+1,ny+1);
eta=zeros(nx+1,ny+1);
% dxx=1.0;
% dyy=1.0;

for i=1:nx+1
    for j=1:ny+1
        ksi(i,j)=dxx*(i-1);
        eta(i,j)=dyy*(j-1);
    end
end

%Calculate x-ksi in phys domain
for i=1:nx+1
        xx(i,:)=x(1,i);
end

%subplot(2,2,2);
%cksi=x;  %x in real space
nyr=round((ny)/2);
for i=1:ny+1
    if i<=nyr
        ceta=(i-1)*y/nyr; %y in real space
        yy(:,i)=ceta(1,:);
        %plot(cksi,ceta); hold on  %plot horizontal y lines
    else
        %disp('>');
        ceta=ymax+(ny-i+1)*y1/(ny-nyr); %y in real space
        yy(:,i)=ceta(1,:);
        %plot(cksi,ceta); hold on  %plot horizontal y lines
    end
%     if i==ny+1 %plot vertical lines for x
%         for j=1:nx
%             line([x(j),x(j)],[0,ceta(j)]); hold on
%         endend
%     end
end

% pmln=zeros(npml,ny+2*npml+1);
% pmls=zeros(npml,ny+2*npml+1);
% pmlw=zeros(nx+1, npml);
% pmle=zeros(nx+1, npml);
% 
% %fill the pmls with it's coordinates
% 
% 
% %Left and right PMLs
% ksi=[pmlw ksi pmle];
% eta=[pmlw eta pmle];
% 
% % %Top and bottom PMLs
%  ksi=[pmls; ksi; pmln];
%  eta=[pmls; eta; pmln];

% for i=1:nx+1
%     for j=1:ny+1
%         scatter(ksi(i,j),eta(i,j))
%     end
% end

%-------------------------------------------------------------------------
if showornot
    %Plot Cartesian grid
    subplot(1,2,2)
    title('Cartesian grid')
    %plot grid lines
    for i=2:nx+1
        for j=2:ny+1
            line([ksi(i-1,j-1) ksi(i-1,j-1)],[eta(i-1,j-1) eta(i-1,j)]); hold on; %vertical
            line([ksi(i-1,j-1) ksi(i,j-1)],[eta(i-1,j-1) eta(i-1,j-1)]); hold on; %horizontal
        end
    end
    line([ksi(1,1) ksi(end,1)],[max(eta(1,:)) max(eta(1,:))]);
    line([ksi(end,1) ksi(end,1)],[eta(1,1) eta(1,end)]);
    %--------------------------------------------------------------------------
    %Plot curvilinear coordinates
    subplot(1,2,1)
    %plot grid lines
    for i=1:ny+1
        plot(xx(:,i),yy(:,i)); hold on
    end
    for i=1:nx+1
        line([xx(i,1), xx(i,1)],[0, yy(i,end)]); hold on
    end
    title('Curvilinear grid')
end

%--------------------------------------------------------------------------
%Derivatives and Jacobian
J=cell(nx,ny);
for i=2:nx
    for j=2:ny
        dksi_dx=(ksi(i+1,j)-ksi(i-1,j))/(xx(i+1,j)-xx(i-1,j));
        dksi_dy=(ksi(i,j+1)-ksi(i,j-1))/(yy(i,j+1)-yy(i,j-1));
        deta_dx=(eta(i+1,j)-eta(i-1,j))/(xx(i+1,j)-xx(i-1,j));
        deta_dy=(eta(i,j+1)-eta(i,j-1))/(yy(i,j+1)-yy(i,j-1));
        J{i-1,j-1}=[dksi_dy deta_dy; dksi_dx deta_dx];
        
%         dx_dksi=(xx(i+1,j)-xx(i-1,j))/(ksi(i+1,j)-ksi(i-1,j));
%         dx_deta=(xx(i,j+1)-xx(i,j-1))/(eta(i,j+1)-eta(i,j-1));
%         dy_dksi=(yy(i+1,j)-yy(i-1,j))/(ksi(i+1,j)-ksi(i-1,j));
%         dy_deta=(yy(i,j+1)-yy(i,j-1))/(eta(i,j+1)-eta(i,j-1));
%         J{i-1,j-1}=[dx_deta dy_deta; dx_dksi dy_dksi];

    end
end
fprintf('Jacobian cell-array %d x %d created\n',size(J,1),size(J,2));
disp('Finished');
end