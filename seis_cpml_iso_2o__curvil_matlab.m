%
% SEISMIC_CPML Version 1.1.1, November 2009.
%
% Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
% Contributor: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
%
% This software is a computer program whose purpose is to solve
% the two-dimensional isotropic elastic wave equation
% using a finite-difference method with Convolutional Perfectly Matched
% Layer (C-PML) conditions.
%
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software. You can use,
% modify and/or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info".
%
% As a counterpart to the access to the source code and rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty and the software's author, the holder of the
% economic rights, and the successive licensors have only limited
% liability.
%
% In this respect, the user's attention is drawn to the risks associated
% with loading, using, modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean that it is complicated to manipulate, and that also
% therefore means that it is reserved for developers and experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or
% data to be ensured and, more generally, to use and operate it in the
% same conditions as regards security.
%
% The full text of the license is available at the end of this program
% and in file "LICENSE".

% 2D elastic finite-difference code in velocity and stress formulation
% with Convolutional-PML (C-PML) absorbing conditions for an isotropic medium

% Dimitri Komatitsch, University of Pau, France, April 2007.

% The second-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
%
%            ^ y
%            |
%            |
%
%            +-------------------+
%            |                   |
%            |                   |
%            |                   |
%            |                   |
%            |        v_y        |
%   sigma_xy +---------+         |
%            |         |         |
%            |         |         |
%            |         |         |
%            |         |         |
%            |         |         |
%            +---------+---------+  ---> x
%           v_x    sigma_xx
%                  sigma_yy
%

% The C-PML implementation is based in part on formulas given in Roden and Gedney (2000).
% If you use this code for your own research, please cite some (or all) of these
% articles:
%
% @ARTICLE{MaKoEz08,
% author = {Roland Martin and Dimitri Komatitsch and Abdela\^aziz Ezziani},
% title = {An unsplit convolutional perfectly matched layer improved at grazing
% incidence for seismic wave equation in poroelastic media},
% journal = {Geophysics},
% year = {2008},
% volume = {73},
% pages = {T51-T61},
% number = {4},
% doi = {10.1190/1.2939484}}
%
% @ARTICLE{MaKo09,
% author = {Roland Martin and Dimitri Komatitsch},
% title = {An unsplit convolutional perfectly matched layer technique improved
% at grazing incidence for the viscoelastic wave equation},
% journal = {Geophysical Journal International},
% year = {2009},
% volume = {179},
% pages = {333-344},
% number = {1},
% doi = {10.1111/j.1365-246X.2009.04278.x}}
%
% @ARTICLE{MaKoGe08,
% author = {Roland Martin and Dimitri Komatitsch and Stephen D. Gedney},
% title = {A variational formulation of a stabilized unsplit convolutional perfectly
% matched layer for the isotropic or anisotropic seismic wave equation},
% journal = {Computer Modeling in Engineering and Sciences},
% year = {2008},
% volume = {37},
% pages = {274-304},
% number = {3}}
%
% @ARTICLE{KoMa07,
% author = {Dimitri Komatitsch and Roland Martin},
% title = {An unsplit convolutional {P}erfectly {M}atched {L}ayer improved
%          at grazing incidence for the seismic wave equation},
% journal = {Geophysics},
% year = {2007},
% volume = {72},
% number = {5},
% pages = {SM155-SM167},
% doi = {10.1190/1.2757586}}
%
% The original CPML technique for Maxwell's equations is described in:
%
% @ARTICLE{RoGe00,
% author = {J. A. Roden and S. D. Gedney},
% title = {Convolution {PML} ({CPML}): {A}n Efficient {FDTD} Implementation
%          of the {CFS}-{PML} for Arbitrary Media},
% journal = {Microwave and Optical Technology Letters},
% year = {2000},
% volume = {27},
% number = {5},
% pages = {334-339},
% doi = {10.1002/1098-2760(20001205)27:5<334::AID-MOP14>3.0.CO;2-A}}

%
% To display the 2D results as color images, use:
%
%   " display image*.gif " or " gimp image*.gif "
%
% or
%
%   " montage -geometry +0+3 -rotate 90 -tile 1x21 image*Vx*.gif allfiles_Vx.gif "
%   " montage -geometry +0+3 -rotate 90 -tile 1x21 image*Vy*.gif allfiles_Vy.gif "
%   then " display allfiles_Vx.gif " or " gimp allfiles_Vx.gif "
%   then " display allfiles_Vy.gif " or " gimp allfiles_Vy.gif "
%

% IMPORTANT : all our CPML codes work fine in single precision as well (which is significantly faster).
%             If you want you can thus force automatic conversion to single precision at compile time
%             or change all the declarations and constants in the code from double precision to single.


%Translated to MATLAB from FORTRAN by O.Ovcharenko, 2015

  tic;  %start timer 
  close all;  %close all extra windows
  clc;  %clear console
  clear all; %clear all variables
  
% total number of grid points in each direction of the grid

%Fine
 NX =282;  %X
 NY =141;  %Y
 
 %NX=621;
 %NY=342;
 
 %NX=121;
 %NY=102;
 
%Rough, but not so bad
 %NX=86;
 %NY=42;

%Very rough
%NX=43;
%NY=21;
 

%--------------------------------------------------------------------------
%---------------------- FLAGS ---------------------------------------------
%Take instant snapshot
    SNAPSHOT=false;
    snapshot_time=150;

%Use explosive source or gaussian?
    EXPLOSIVE_SOURCE=true;
    
%Show source position on the graph?
    SHOW_SOURCE_POSITION=true;

%Pause a little bit each iteration
   PAUSE_ON=false;
   pause_time=0.01; %[sec]
   
% To show or don't show wavefield 
  SAVE_VX_JPG =false;
  SAVE_VY_JPG =true;

%Record video - corresponding SAVE_VX or VY must be turned on
%because video is being created by capturing of current frame
%Matlab 2012 + required, saves video to a current folder
  MAKE_MOVIE_VX=false;
  MAKE_MOVIE_VY=false;
  
%Apply flat Free surface
  FREE_SURFACE=false;
  %Position of flat horizontal free surface
  NFS = NY-round(NY/2)+round(NY/3);
  
%Apply flat Flat elastic boundary
  FE_BOUNDARY=true;
  cp_above_eb=1800.d0;
  cp_below_eb=3300.d0;
  rho_above_eb=2400.d0;
  rho_below_eb=3200.d0;
  
% flags to add PML layers to the edges of the grid
  USE_PML_XMIN = true;
  USE_PML_XMAX = true;
  USE_PML_YMIN = true;
  USE_PML_YMAX = true;
  
  
  
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
  
  %Check if it is possible to save video
  if SAVE_VY_JPG
    fig_vy=figure;
  else
      %to get rid of errors that can occur
      disp('Error. It is necesary to have SAVE_VY_JPG=true.');
      MAKE_MOVIE_VY=false;
  end
  
  if SAVE_VX_JPG
    fig_vx=figure;
  else
      %to get rid of errors that can occur
      disp('Error. It is necesary to have SAVE_VX_JPG=true.');
      MAKE_MOVIE_VX=false;
  end
  

 % size of a grid cell
 DELTAX = 10.d0;	%[m]
 DELTAY = DELTAX;
    
 %vectors for visualisation using imagesec
 nx_vec=1:NX*DELTAX;	%[m]
 ny_vec=1:NY*DELTAY;

 
 % thickness of the PML layer in grid points
 NPOINTS_PML = 10;

 % P-velocity, S-velocity and density
 cp = 3300.d0;	%[km/s]
 cs = cp / 1.732d0;	%[km/s]
 density = 2800.d0;	%[kg/m3]

% total number of time steps
 NSTEP = 500;

% time step in seconds
  DELTAT = 2.d-3;	%[sec]

% parameters for the source
  f0 = 20.d0;
  t0 = 1.20d0 / f0;
  factor = 1.d7;

% source
  %ISOURCE = NX - 2*NPOINTS_PML - 1-round(NX/3);
  ISOURCE = NX - round(NX/3);
  JSOURCE = (NY / 3) + 1;
   
  xsource = (ISOURCE - 1) * DELTAX;
  ysource = (JSOURCE - 1) * DELTAY;
% angle of source force clockwise with respect to vertical (Y) axis
  ANGLE_FORCE = 45.d0;

% receivers
  NREC = 2;
  xdeb = xsource - 100.d0;   % first receiver x in meters
  ydeb = 2300.d0;            % first receiver y in meters
  xfin = xsource;            % last receiver x in meters
  yfin =  300.d0;           % last receiver y in meters

% display information on the screen from time to time
  IT_DISPLAY = 5;

% value of PI
  PI = 3.141592653589793238462643d0;

% conversion from degrees to radians
  DEGREES_TO_RADIANS = PI / 180.d0;

% zero
  ZERO = 0.d0;

% large value for maximum
  HUGEVAL = 1.d+30;

% velocity threshold above which we consider that the code became unstable
  STABILITY_THRESHOLD = 1.d+25;

% main arrays
  vx=zeros(NX,NY);
  vy=zeros(NX,NY);
  sigmaxx=zeros(NX,NY);
  sigmayy=zeros(NX,NY);
  sigmaxy=zeros(NX,NY);
  lambda=zeros(NX,NY);
  mu=zeros(NX,NY);
  rho=zeros(NX,NY);

  total_energy_kinetic=zeros(NSTEP);
  total_energy_potential=zeros(NSTEP);

% power to compute d0 profile
  NPOWER = 2.d0;

  K_MAX_PML = 1.d0; % from Gedney page 8.11
  ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0); % from Festa and Vilotte

% arrays for the memory variables
% could declare these arrays in PML only to save a lot of memory, but proof of concept only here
  memory_dvx_dx=zeros(NX,NY);
  memory_dvx_dy=zeros(NX,NY);
  memory_dvy_dx=zeros(NX,NY);
  memory_dvy_dy=zeros(NX,NY);
  memory_dsigmaxx_dx=zeros(NX,NY);
  memory_dsigmayy_dy=zeros(NX,NY);
  memory_dsigmaxy_dx=zeros(NX,NY);
  memory_dsigmaxy_dy=zeros(NX,NY);


 % 1D arrays for the damping profiles
 d_x=zeros(NX,1);
 K_x=zeros(NX,1);
 alpha_x=zeros(NX,1);
 a_x=zeros(NX,1);
 b_x=zeros(NX,1);
 d_x_half=zeros(NX,1);
 K_x_half=zeros(NX,1);
 alpha_x_half=zeros(NX,1);
 a_x_half=zeros(NX,1);
 b_x_half=zeros(NX,1);

 d_y=zeros(NY,1);
 K_y=zeros(NY,1);
 alpha_y=zeros(NY,1);
 a_y=zeros(NY,1);
 b_y=zeros(NY,1);
 d_y_half=zeros(NY,1);
 K_y_half=zeros(NY,1);
 alpha_y_half=zeros(NY,1);
 a_y_half=zeros(NY,1);
 b_y_half=zeros(NY,1);

 % for receivers
 ix_rec=zeros(NREC,1);
 iy_rec=zeros(NREC,1);
 xrec=zeros(NREC,1);
 yrec=zeros(NREC,1); 

 % for seismograms
 sisvx=zeros(NSTEP,NREC);
 sisvy=zeros(NSTEP,NREC);

 %character(len=100) :: filename
 
   %Initiate video object for vx
  if MAKE_MOVIE_VX
	  movie_name_vx=['vy_video_' num2str(NX) '_' num2str(NY) '_' num2str(DELTAX) '_' num2str(f0) '.avi'];
	  vidObj_vx=VideoWriter(movie_name_vx);
	  open(vidObj_vx);
  end

  %Initiate video object for vy
  if MAKE_MOVIE_VY
	  movie_name_vy=['vy_video_' num2str(NX) '_' num2str(NY) '_' num2str(DELTAX) '_' num2str(f0) '.avi'];
	  vidObj_vy=VideoWriter(movie_name_vy);
	  open(vidObj_vy);
  end
 
 %----------------------------------------
 %--- program starts here ----------------
 %----------------------------------------

    fprintf('2D elastic finite-difference code in velocity and stress formulation with C-PML\n\n');
    
    %  display size of the model
    fprintf('NX = %d\n',NX);
    fprintf('NY = %d\n\n',NY);
    fprintf('size of the model along X = %.2f\n',(NX - 1) * DELTAX);
    fprintf('size of the model along Y = %.2f\n\n',(NY - 1) * DELTAY);
    fprintf('Total number of grid points = %.2f\n\n',NX * NY);
    
%--------------------------------------------------------------------------
%---------Constructing of Curvilinear and Cartesian meshes-----------------
%--------------------------------------------------------------------------
if FE_BOUNDARY
    fprintf('func_curv_jacob started...\n');
    [ksi,eta, xx,yy,J, Ji] = func_curv_jacob(NX,NY,0,NX*DELTAX, 0, NY*DELTAY,'-(1.25*pi*x/max(x)+0.25*pi)',DELTAX,DELTAY,false);
    fprintf('func_curv_jacob finished\n\n');
end
    
    

%--- define profile of absorption in PML region ---

% thickness of the PML layer in meters
  thickness_PML_x = NPOINTS_PML * DELTAX;
  thickness_PML_y = NPOINTS_PML * DELTAY;

% reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.001d0;

% check that NPOWER is okaymarkers=zeros(nx+1,ny+1);
  if(NPOWER < 1)       
      disp('NPOWER must be greater than 1');
      break;
  end

% compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - (NPOWER + 1) * cp * log(Rcoef) / (2.d0 * thickness_PML_x);
  d0_y = - (NPOWER + 1) * cp * log(Rcoef) / (2.d0 * thickness_PML_y);
 
  fprintf('d0_x = %.2f\n',d0_x);
  fprintf('d0_y = %.2f\n\n',d0_y);

  d_x(:) = ZERO;
  d_x_half(:) = ZERO;
  K_x(:) = 1.d0;
  K_x_half(:) = 1.d0;
  alpha_x(:) = ZERO;
  alpha_x_half(:) = ZERO;
  a_x(:) = ZERO;
  a_x_half(:) = ZERO;

  d_y(:) = ZERO;
  d_y_half(:) = ZERO;
  K_y(:) = 1.d0;
  K_y_half(:) = 1.d0;
  alpha_y(:) = ZERO;
  alpha_y_half(:) = ZERO;
  a_y(:) = ZERO;
  a_y_half(:) = ZERO;
  %break;
%--------------------------------------------------------------------------
% damping in the X direction

% origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x;
  xoriginright = (NX-1)*DELTAX - thickness_PML_x;

for i = 1:NX
    % abscissa of current grid point along the damping profile
    xval = DELTAX * double(i-1);
    %---------- left edge
    if(USE_PML_XMIN)
        % define damping profile at the grid points
        abscissa_in_PML = xoriginleft - xval;
        if(abscissa_in_PML >= ZERO)
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            d_x(i) = d0_x * abscissa_normalized^NPOWER;
            % this taken from Gedney page 8.2
            K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
            alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
        end
        % define damping profile at half the grid points
        abscissa_in_PML = xoriginleft - (xval + DELTAX/2.d0);
        if(abscissa_in_PML >= ZERO)
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            d_x_half(i) = d0_x * abscissa_normalized^NPOWER;
            % this taken from Gedney page 8.2
            K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
            alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
        end
    end

%---------- right edge
   if(USE_PML_XMAX)
        % define damping profile at the grid points
        abscissa_in_PML = xval - xoriginright;
        if(abscissa_in_PML >= ZERO)
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            d_x(i) = d0_x * abscissa_normalized^NPOWER;
            % this taken from Gedney page 8.2
            K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
            alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
        end

        % define damping profile at half the grid points
        abscissa_in_PML = xval + DELTAX/2.d0 - xoriginright;
        if(abscissa_in_PML >= ZERO)
            abscissa_normalized = abscissa_in_PML / thickness_PML_x;
            d_x_half(i) = d0_x * abscissa_normalized^NPOWER;
            % this taken from Gedney page 8.2
            K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
            alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
        end
   end

    % just in case, for -5 at the end
    if(alpha_x(i) < ZERO) 
        alpha_x(i) = ZERO;
    end
    if(alpha_x_half(i) < ZERO) 
        alpha_x_half(i) = ZERO;
    end
    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT);
    b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DELTAT);

    % this to avoid division by zero outside the PML
    if(abs(d_x(i)) > 1.d-6) 
        a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)));
    end    
    if(abs(d_x_half(i)) > 1.d-6) 
        a_x_half(i)=d_x_half(i)*(b_x_half(i) - 1.d0) / (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_x_half(i)));
    end
end

%--------------------------------------------------------------------------
% damping in the Y direction

% origin of the PML layer (position of right edge minus thickness, in meters)
  yoriginbottom = thickness_PML_y;
  yorigintop = (NY-1)*DELTAY - thickness_PML_y;

  for j = 1:NY
    % abscissa of current grid point along the damping profile
    yval = DELTAY * double(j-1);
    %---------- bottom edge
    if(USE_PML_YMIN)
      % define damping profile at the grid points
      abscissa_in_PML = yoriginbottom - yval;
      if(abscissa_in_PML >= ZERO)
        abscissa_normalized = abscissa_in_PML / thickness_PML_y;
        d_y(j) = d0_y * abscissa_normalized^NPOWER;
        % this taken from Gedney page 8.2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
      end

      % define damping profile at half the grid points
      abscissa_in_PML = yoriginbottom - (yval + DELTAY/2.d0);
      if(abscissa_in_PML >= ZERO)
        abscissa_normalized = abscissa_in_PML / thickness_PML_y;
        d_y_half(j) = d0_y * abscissa_normalized^NPOWER;
        % this taken from Gedney page 8.2
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
      end
    end

%---------- top edge
    if(USE_PML_YMAX)
      % define damping profile at the grid points
      abscissa_in_PML = yval - yorigintop;
      if(abscissa_in_PML >= ZERO)
        abscissa_normalized = abscissa_in_PML / thickness_PML_y;
        d_y(j) = d0_y * abscissa_normalized^NPOWER;
        % this taken from Gedney page 8.2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
      end

      % define damping profile at half the grid points
      abscissa_in_PML = yval + DELTAY/2.d0 - yorigintop;
      if(abscissa_in_PML >= ZERO)
        abscissa_normalized = abscissa_in_PML / thickness_PML_y;
        d_y_half(j) = d0_y * abscissa_normalized^NPOWER;
        % this taken from Gedney page 8.2
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized^NPOWER;
        alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML;
      end
    end

    b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_y(j)) * DELTAT);
    b_y_half(j) = exp(- (d_y_half(j) / K_y_half(j) + alpha_y_half(j)) * DELTAT);

    % this to avoid division by zero outside the PML
    if(abs(d_y(j)) > 1.d-6) 
        a_y(j) = d_y(j) * (b_y(j) - 1.d0) / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)));
    end  
    
    if(abs(d_y_half(j)) > 1.d-6) 
        a_y_half(j) = d_y_half(j)*(b_y_half(j) - 1.d0) / (K_y_half(j) * (d_y_half(j) + K_y_half(j) * alpha_y_half(j)));
    end 
  end
  
  
  ymtr=NY*DELTAY; %ymax [m]
  xmtr=NX*DELTAX; %xmax [m]
  
  %xdscr=[0:NX]*DELTAX;
  xdscr=linspace(0,xmtr,20*NX);
  phi=-(1.25*PI*xdscr/max(xdscr)+0.25*PI);
  ydscr=sin(phi);
  ydscr=ymtr*ydscr/4;
  ydscr=abs(min(ydscr))+ydscr+ymtr/4;

  %------------------------------------------------------------------------
  if FE_BOUNDARY
        %grid point coordinates in physical domain
        gr_x=zeros(NX+1,NY+1);
        gr_y=zeros(NX+1,NY+1);

        %Cartesian grid
        for i=1:NX+1
            for j=1:NY+1
                gr_x(i,j)=(i-1)*DELTAX;
                gr_y(i,j)=(j-1)*DELTAY;    
            end    
        end
  end 
  
  % compute the Lame parameters and density

  for j = 1:NY
    y_trial=j*DELTAY;
    for i = 1:NX
        x_trial = i*DELTAX;
        rho(i,j) = density;
        mu(i,j) = density*cs*cs;
        lambda(i,j) = density*(cp*cp - 2.d0*cs*cs);
        %If we use flast elastic inner boundary
        
        if FE_BOUNDARY
            temp_dn=round(length(ydscr)/NX); %-number of topography points in dx sample
            if y_trial>=ydscr(1,i*temp_dn)
                 rho(i,j) = rho_above_eb;
                 cp = cp_above_eb;	%[km/s]
                 cs = cp / 1.732d0;	%[km/s]
                 lambda(i,j) =density*(cp*cp - 2.d0*cs*cs);
                 mu(i,j) = density*cs*cs;
            else
                 rho(i,j) = rho_below_eb;
                 cp = cp_below_eb;	%[km/s]
                 cs = cp / 1.732d0;	%[km/s]
                 lambda(i,j) =density*(cp*cp - 2.d0*cs*cs);
                 mu(i,j) = density*cs*cs;
            end
        end
    end
  end
  
% % Create Cijkl matrix
% fprintf('\nStart Cijkl');
% if FE_BOUNDARY
%     order=2;
%     cijkla=zeros(order,order,order,order);
%     cijklb=zeros(order,order,order,order);
%     densitya = rho_above_eb;
%     cpa = cp_above_eb;	%[km/s]
%     csa = cp / 1.732d0;	%[km/s]
%     lambdaa =density*(cp*cp - 2.d0*cs*cs);
%     mua = density*cs*cs;
%     densityb = rho_below_eb;
%     cpb = cp_below_eb;	%[km/s]
%     csb = cp / 1.732d0;	%[km/s]
%     lambdab =density*(cp*cp - 2.d0*cs*cs);
%     mub = density*cs*cs;
%     for i=1:order
%         for j=1:order
%             for k=1:order
%                 for l=1:order
%                      cijkla(i,j,k,l)=lambdaa*dkr(i,j)*dkr(k,l)+mua*(dkr(i,k)*dkr(j,l)+dkr(i,l)*dkr(j,k));
%                      cijklb(i,j,k,l)=lambdab*dkr(i,j)*dkr(k,l)+mub*(dkr(i,k)*dkr(j,l)+dkr(i,l)*dkr(j,k));
%                      fprintf('.');
%                 end
%             end
%             fprintf('\n');
%         end
%     end
%     clearvars densitya cpa csa lambdaa mua densityb cpb csb lambdab mub;
%     fprintf('Cijkl dim %d created\n',order);
% end
%------------------------------------------------------------------------
  
  % print position of the source
   fprintf('Position of the source:\n');
   fprintf('x = %.2f\n',xsource);
   fprintf('y = %.2f\n\n',ysource);
 
%  define location of receivers
   fprintf('There are %d receivers\n',NREC);

   xspacerec = (xfin-xdeb) / double(NREC-1);
   yspacerec = (yfin-ydeb) / double(NREC-1);
   for irec=1:NREC
     xrec(irec) = xdeb + double(irec-1)*xspacerec;
     yrec(irec) = ydeb + double(irec-1)*yspacerec;
   end

% find closest grid point for each receiver
   for irec=1:NREC
   dist = HUGEVAL;
   for j = 1:NY
    for i = 1:NX
      distval = sqrt((DELTAX*double(i-1) - xrec(irec))^2 + (DELTAY*double(j-1) - yrec(irec))^2);
      if(distval < dist)
        dist = distval;
        ix_rec(irec) = i;
        iy_rec(irec) = j;
      end
    end
   end
   fprintf('receiver %d x_target,y_target = %.2f  %.2f\n',irec,xrec(irec),yrec(irec))
   fprintf('closest grid point found at distance %.2f in i,j = %d  %d\n\n',dist,ix_rec(irec),iy_rec(irec));
   end

%--------------------------------------------------------------------------
% check the Courant stability condition for the explicit time scheme
% R. Courant et K. O. Friedrichs et H. Lewy (1928)
  Courant_number = cp * DELTAT * sqrt(1.d0/DELTAX^2 + 1.d0/DELTAY^2);
  fprintf('Courant number = %.4f\n\n',Courant_number);

  if(Courant_number > 1.d0) 
      disp('time step is too large, simulation will be unstable');
      break;
  end
  
  if FE_BOUNDARY
      %calculate involved grid points, descritized coordinates of curve,
      %normal vectors, coordinates of middles of the descritized samples.
      %All the output variables are vectors
      [markers, xt_dis, yt_dis, nvec, xmn, ymn] = func_find_closest_grid_nodes(NX,NY,1,gr_x,gr_y ,xdscr, ydscr);
      %[markers, xt_dis, yt_dis, nvec, xmn, ymn] = func_c_find_closest_grid_nodes(NX,NY,1,gr_x,gr_y ,xdscr, ydscr, DELTAX);
      nvecx=nvec*[1 0]';
      nvecy=nvec*[0 1]';
  end

%--------------------------------------------------------------------------
% initialize arrays
  vx(:,:) = ZERO;
  vy(:,:) = ZERO;
  sigmaxx(:,:) = ZERO;
  sigmayy(:,:) = ZERO;
  sigmaxy(:,:) = ZERO;

% PML
  memory_dvx_dx(:,:) = ZERO;
  memory_dvx_dy(:,:) = ZERO;
  memory_dvy_dx(:,:) = ZERO;
  memory_dvy_dy(:,:) = ZERO;
  memory_dsigmaxx_dx(:,:) = ZERO;
  memory_dsigmayy_dy(:,:) = ZERO;
  memory_dsigmaxy_dx(:,:) = ZERO;
  memory_dsigmaxy_dy(:,:) = ZERO;
  memory_dsigmayy_dx(:,:) = ZERO;

% initialize seismograms
  sisvx(:,:) = ZERO;
  sisvy(:,:) = ZERO;

% initialize total energy
  total_energy_kinetic(:) = ZERO;
  total_energy_potential(:) = ZERO;

  input('Press Enter to start time loop ...');
    
  %Set red-blue colormap for images
  CMAP=zeros(256,3);
  c1=[0 0 1]; %blue
  c2=[1 1 1]; %white
  c3=[1 0 0]; %red
  for nc=1:128
	  f=(nc-1)/128;
	  c=(1-sqrt(f))*c1+sqrt(f)*c2;
	  CMAP(nc,:)=c;
	  c=(1-f^2)*c2+f^2*c3;
	  CMAP(128+nc,:)=c;
	end
  colormap(CMAP);
  set(gca,'YDir','normal');
%---------------------------------
%---  beginning of time loop -----
%---------------------------------
  for it = 1:NSTEP
%------------------------------------------------------------
% compute stress sigma and update memory variables for C-PML
%------------------------------------------------------------
        for j = 2:NY
            for i = 1:NX-1
              % interpolate material parameters at the right location in the staggered grid cell
              lambda_half_x = 0.5d0 * (lambda(i+1,j) + lambda(i,j));
              mu_half_x = 0.5d0 * (mu(i+1,j) + mu(i,j));
              lambda_plus_two_mu_half_x = lambda_half_x + 2.d0 * mu_half_x;

              value_dvx_dx = (vx(i+1,j) - vx(i,j)) / DELTAX;
              value_dvy_dy = (vy(i,j) - vy(i,j-1)) / DELTAY;

              memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + a_x_half(i) * value_dvx_dx;
              memory_dvy_dy(i,j) = b_y(j) * memory_dvy_dy(i,j) + a_y(j) * value_dvy_dy;

              value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j);
              value_dvy_dy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j);

              sigmaxx(i,j) = sigmaxx(i,j) + ...
                 (lambda_plus_two_mu_half_x * value_dvx_dx + lambda_half_x * value_dvy_dy) * DELTAT;

              sigmayy(i,j) = sigmayy(i,j) + ...
                 (lambda_half_x * value_dvx_dx + lambda_plus_two_mu_half_x * value_dvy_dy) * DELTAT;
               %--------------------------------------------------------------------------------------------------
               if(FREE_SURFACE)
                    if(j==NFS+1)
                        sigmayy(i,j)=-sigmayy(i,j-1);
                        sigmaxx(i,j)=ZERO;
                    end
               end    
               %--------------------------------------------------------------------------------------------------
            end 
        end
    

    for j = 1:NY-1
        for i = 2:NX
              % interpolate material parameters at the right location in the staggered grid cell
              mu_half_y = 0.5d0 * (mu(i,j+1) + mu(i,j));

              value_dvy_dx = (vy(i,j) - vy(i-1,j)) / DELTAX;
              value_dvx_dy = (vx(i,j+1) - vx(i,j)) / DELTAY;

              memory_dvy_dx(i,j) = b_x(i) * memory_dvy_dx(i,j) + a_x(i) * value_dvy_dx;
              memory_dvx_dy(i,j) = b_y_half(j) * memory_dvx_dy(i,j) + a_y_half(j) * value_dvx_dy;

              value_dvy_dx = value_dvy_dx / K_x(i) + memory_dvy_dx(i,j);
              value_dvx_dy = value_dvx_dy / K_y_half(j) + memory_dvx_dy(i,j);

              sigmaxy(i,j) = sigmaxy(i,j) + mu_half_y * (value_dvy_dx + value_dvx_dy) * DELTAT;
       end
    end
    if(FREE_SURFACE)
        sigmaxy(:,NFS)=ZERO;
    end

    %--------------------------------------------------------
    % compute velocity and update memory variables for C-PML
    %--------------------------------------------------------

    for j = 2:NY
        for i = 2:NX
              value_dsigmaxx_dx = (sigmaxx(i,j) - sigmaxx(i-1,j)) / DELTAX;
              value_dsigmaxy_dy = (sigmaxy(i,j) - sigmaxy(i,j-1)) / DELTAY;

              memory_dsigmaxx_dx(i,j) = b_x(i) * memory_dsigmaxx_dx(i,j) + a_x(i) * value_dsigmaxx_dx;
              memory_dsigmaxy_dy(i,j) = b_y(j) * memory_dsigmaxy_dy(i,j) + a_y(j) * value_dsigmaxy_dy;

              value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + memory_dsigmaxx_dx(i,j);
              value_dsigmaxy_dy = value_dsigmaxy_dy / K_y(j) + memory_dsigmaxy_dy(i,j);

              
              if (FE_BOUNDARY && markers(i,j)==1)
                  vx(i,j) = vx(i,j) + (value_dsigmaxx_dx + value_dsigmaxy_dy) * DELTAT / rho(i,j);
                  %fprintf('Vx FE_B: %d %d\n',i,j);
%                       Bpx=[1 0 0 0 0 0; ...
%                           0 0 nvx*cijklb(:,1,1,2)+nvy*cijklb(:,2,1,2) 0 0 0; ...
%                           0 nvy -nvx 0 0 0; ...
%                           0 0 0 0 -nvx^2*cijklb(:,1,1,2) nvx*nvy*cijklb(:,1,1,2); ...
%                           0 0 0 0 cijklb(:,2,1,2)/rhob cijklb(:,1,1,2)/rhob;...
%                           0 0 0 nvy^2 nvx^2 -2*nvx*nvy]
              else
                  vx(i,j) = vx(i,j) + (value_dsigmaxx_dx + value_dsigmaxy_dy) * DELTAT / rho(i,j);
              end
        end
    end

    for j = 1:NY-1
        for i = 1:NX-1
          % interpolate density at the right location in the staggered grid cell
          rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1));
          
          value_dsigmaxy_dx = (sigmaxy(i+1,j) - sigmaxy(i,j)) / DELTAX;
          value_dsigmayy_dy = (sigmayy(i,j+1) - sigmayy(i,j)) / DELTAY;
       
          memory_dsigmaxy_dx(i,j) = b_x_half(i) * memory_dsigmaxy_dx(i,j) + a_x_half(i) * value_dsigmaxy_dx;
          memory_dsigmayy_dy(i,j) = b_y_half(j) * memory_dsigmayy_dy(i,j) + a_y_half(j) * value_dsigmayy_dy;

          value_dsigmaxy_dx = value_dsigmaxy_dx / K_x_half(i) + memory_dsigmaxy_dx(i,j);
          value_dsigmayy_dy = value_dsigmayy_dy / K_y_half(j) + memory_dsigmayy_dy(i,j);  

          vy(i,j) = vy(i,j) + (value_dsigmaxy_dx + value_dsigmayy_dy) * DELTAT / rho_half_x_half_y;
        end
        if(FREE_SURFACE)
            if(j==NFS)
                sigmayy(:,NFS+1)=ZERO; %Once mirrored sigmayy(:,NFS+1) was used, we have to set it equal to zero
            end                              % to get rid of artificial values above the free surface
        end
    end

    % add the source (force vector located at a given grid point)
    a = pi*pi*f0*f0;
    t = double(it-1)*DELTAT;

    % Gaussian
     %source_term = factor * exp(-a*(t-t0)^2);

     %source_term = factor * (t-t0);
    
     % first derivative of a Gaussian
        source_term =  -factor*2.d0*a*(t-t0)*exp(-a*(t-t0)^2);
    

    % Ricker source time function (second derivative of a Gaussian)
     %source_term = factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);

    force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;
    force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term;

    % define location of the source
    i = ISOURCE;
    j = JSOURCE;

    % interpolate density at the right location in the staggered grid cell
    rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1));
    
    if EXPLOSIVE_SOURCE
        source_term = factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);
        force_x = source_term;
        force_y = source_term;
        sigmaxx(i,j)=sigmaxx(i,j)+force_x;
        sigmayy(i,j)=sigmayy(i,j)+force_y;
    else
        vx(i,j) = vx(i,j) + force_x * DELTAT / rho(i,j);
        vy(i,j) = vy(i,j) + force_y * DELTAT / rho_half_x_half_y;
    end
    
    % Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
    vx(1,:) = ZERO;
    vx(NX,:) = ZERO;

    vx(:,1) = ZERO;
    vx(:,NY) = ZERO;

    vy(1,:) = ZERO;
    vy(NX,:) = ZERO;

    vy(:,1) = ZERO;
    vy(:,NY) = ZERO;

    % store seismograms
    for irec = 1:NREC
        sisvx(it,irec) = vx(ix_rec(irec),iy_rec(irec));
        sisvy(it,irec) = vy(ix_rec(irec),iy_rec(irec));
    end

    % compute total energy in the medium (without the PML layers)

    % compute kinetic energy first, defined as 1/2 rho ||v||^2
    % in principle we should use rho_half_x_half_y instead of rho for vy
    % in order to interpolate density at the right location in the staggered grid cell
    % but in a homogeneous medium we can safely ignore it
    
    %total_energy_kinetic(it) = 0.5d0 .*sum(rho((NPOINTS_PML+1):(NX-NPOINTS_PML),(NPOINTS_PML+1):(NY-NPOINTS_PML))*( ...
    %    vx((NPOINTS_PML+1):(NX-NPOINTS_PML),(NPOINTS_PML+1):(NY-NPOINTS_PML)).^2 +  ...
    %    vy((NPOINTS_PML+1):(NX-NPOINTS_PML),(NPOINTS_PML+1):(NY-NPOINTS_PML)).^2));

    % add potential energy, defined as 1/2 epsilon_ij sigma_ij
    % in principle we should interpolate the medium parameters at the right lo thencation
    % in the staggered grid cell but in a homogeneous medium we can safely ignore it
    total_energy_potential(it) = ZERO;
    for j = NPOINTS_PML+1: NY-NPOINTS_PML
        for i = NPOINTS_PML+1: NX-NPOINTS_PML
            epsilon_xx = ((lambda(i,j) + 2.d0*mu(i,j)) * sigmaxx(i,j) - lambda(i,j) * ...
                sigmayy(i,j)) / (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)));
            epsilon_yy = ((lambda(i,j) + 2.d0*mu(i,j)) * sigmayy(i,j) - lambda(i,j) * ...
                sigmaxx(i,j)) / (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)));
            epsilon_xy = sigmaxy(i,j) / (2.d0 * mu(i,j));
            total_energy_potential(it) = total_energy_potential(it) + ...
                0.5d0 * (epsilon_xx * sigmaxx(i,j) + epsilon_yy * sigmayy(i,j) + 2.d0 * epsilon_xy * sigmaxy(i,j));
        end
    end

    % output information
    if(mod(it,IT_DISPLAY) == 0 || it == 5) 
        % print maximum of norm of velocity
        velocnorm = max(sqrt(vx.^2 + vy.^2));
        fprintf('Time step: %d\n',it)
        fprintf('Time: %.2f seconds\n',single((it-1)*DELTAT));
        %fprintf('Max norm velocity vector V (m/s) = %.2f\n',velocnorm);
        %     print *,'total energy = ',total_energy_kinetic(it) + total_energy_potential(it)
        %     print *
        % check stability of the code, exit if unstable
        if(velocnorm > STABILITY_THRESHOLD)
            break 
            disp('code became unstable and blew up');
        end    

        if(SAVE_VX_JPG)
            clf;
            %fig_vx=figure;
            imagesc(nx_vec,ny_vec,vy');
            colorbar();
            xlabel('m');
            ylabel('m');
            %fig_vx_name=['Snapshot Vx. t= ' num2str(it) 'sec'];
            %set(fig_vx,'Name',fig_vx_name); % - set name of a figure
            %set(fig_vx,'NumberTitle','off'); %- hide "Figure #" in name of figure
            set(gca,'YDir','normal');
            drawnow;
            if MAKE_MOVIE_VX
                F_x=getframe(gcf);  %-  capture figure or use gcf to get current figure. Or get current
                writeVideo(vidObj_vx,F_x);  %- add frame to the movie
                fprintf('Frame for %s captured\n',movie_name_vx);
            end
        end
        
    
        if(SAVE_VY_JPG)
            clf;	%clear current frame
            %fig_vy=figure;
 
            imagesc(nx_vec,ny_vec,vy'); 
            xlabel('m');
            ylabel('m');
            set(gca,'YDir','normal');
            drawnow; hold on;
            
            if FE_BOUNDARY
                plot(xdscr,ydscr); 
                xlabel('m');
                ylabel('m');
                set(gca,'YDir','normal');
                drawnow;
            end
            
            if SHOW_SOURCE_POSITION
                scatter(xsource, ysource,'g','filled'); drawnow;
            end
           
            if SNAPSHOT
                if it==snapshot_time
                    snapshat = getframe(gcf);
                    imgg = frame2im(snapshat);
                    scrsht_name=['im' num2str(it) '.png'];
                    imwrite(imgg,scrsht_name);
                    fprintf('Screenshot %s saved to %s\n', scrsht_name, pwd);
                    clearvars scrsht_name imgg snapshat
                    fprintf('Unused variables removed\n');
                end  
            end
            
            %input('Next');
            

            if MAKE_MOVIE_VY
                F_y=getframe(gcf);  %-  capture figure or use gcf to get current figure. Or get current
                writeVideo(vidObj_vy,F_y);  %- add frame to the movie
                fprintf('Frame for %s captured\n',movie_name_vy);
            end
        end
        fprintf('\n'); 
    end
    if PAUSE_ON
        pause(pause_time);
    end
  end
  % end of time loop
  
  current_folder=pwd;	%current path
  if MAKE_MOVIE_VX
	  close(vidObj_vx);     %- close video file
      printf('Video %s saved in %s\n',movie_name_vx,current_folder);
  end
  
  if MAKE_MOVIE_VY
	  close(vidObj_vy);     %- close video file
      fprintf('Video %s saved in %s\n',movie_name_vy, current_folder);
  end
  
  toc;	%stop timer
  disp('End');
  
  %
% CeCILL FREE SOFTWARE LICENSE AGREEMENT
%
%     Notice
%
% This Agreement is a Free Software license agreement that is the result
% of discussions between its authors in order to ensure compliance with
% the two main principles guiding its drafting:
%
%     * firstly, compliance with the principles governing the distribution
%       of Free Software: access to source code, broad rights granted to
%       users,
%     * secondly, the election of a governing law, French law, with which
%       it is conformant, both as regards the law of torts and
%       intellectual property law, and the protection that it offers to
%       both authors and holders of the economic rights over software.
%
% The authors of the CeCILL (for Ce[a] C[nrs] I[nria] L[ogiciel] L[ibre])
% license are:
%
% Commissariat a l'Energie Atomique - CEA, a public scientific, technical
% and industrial research establishment, having its principal place of
% business at 25 rue Leblanc, immeuble Le Ponant D, 75015 Paris, France.
%
% Centre National de la Recherche Scientifique - CNRS, a public scientific
% and technological establishment, having its principal place of business
% at 3 rue Michel-Ange, 75794 Paris cedex 16, France.
%
% Institut National de Recherche en Informatique et en Automatique -
% INRIA, a public scientific and technological establishment, having its
% principal place of business at Domaine de Voluceau, Rocquencourt, BP
% 105, 78153 Le Chesnay cedex, France.
%
%     Preamble
%
% The purpose of this Free Software license agreement is to grant users
% the right to modify and redistribute the software governed by this
% license within the framework of an open source distribution model.
%
% The exercising of these rights is conditional upon certain obligations
% for users so as to preserve this status for all subsequent redistributions.
%
% In consideration of access to the source code and the rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty and the software's author, the holder of the
% economic rights, and the successive licensors only have limited liability.
%
% In this respect, the risks associated with loading, using, modifying
% and/or developing or reproducing the software by the user are brought to
% the user's attention, given its Free Software status, which may make it
% complicated to use, with the result that its use is reserved for
% developers and experienced professionals having in-depth computer
% knowledge. Users are therefore encouraged to load and test the
% suitability of the software as regards their requirements in conditions
% enabling the security of their systems and/or data to be ensured and,
% more generally, to use and operate it in the same conditions of
% security. This Agreement may be freely reproduced and published,
% provided it is not altered, and that no provisions are either added or
% removed herefrom.
%
% This Agreement may apply to any or all software for which the holder of
% the economic rights decides to submit the use thereof to its provisions.
%
%     Article 1 - DEFINITIONS
%
% For the purpose of this Agreement, when the following expressions
% commence with a capital letter, they shall have the following meaning:
%
% Agreement: means this license agreement, and its possible subsequent
% versions and annexes.
%
% Software: means the software in its Object Code and/or Source Code form
% and, where applicable, its documentation, "as is" when the Licensee
% accepts the Agreement.
%
% Initial Software: means the Software in its Source Code and possibly its
% Object Code form and, where applicable, its documentation, "as is" when
% it is first distributed under the terms and conditions of the Agreement.
%
% Modified Software: means the Software modified by at least one
% Contribution.
%
% Source Code: means all the Software's instructions and program lines to
% which access is required so as to modify the Software.
%
% Object Code: means the binary files originating from the compilation of
% the Source Code.
%
% Holder: means the holder(s) of the economic rights over the Initial
% Software.
%
% Licensee: means the Software user(s) having accepted the Agreement.
%
% Contributor: means a Licensee having made at least one Contribution.
%
% Licensor: means the Holder, or any other individual or legal entity, who
% distributes the Software under the Agreement.
%
% Contribution: means any or all modifications, corrections, translations,
% adaptations and/or new functions integrated into the Software by any or
% all Contributors, as well as any or all Internal Modules.
%
% Module: means a set of sources files including their documentation that
% enables supplementary functions or services in addition to those offered
% by the Software.
%
% External Module: means any or all Modules, not derived from the
% Software, so that this Module and the Software run in separate address
% spaces, with one calling the other when they are run.
%
% Internal Module: means any or all Module, connected to the Software so
% that they both execute in the same address space.
%
% GNU GPL: means the GNU General Public License version 2 or any
% subsequent version, as published by the Free Software Foundation Inc.
%
% Parties: mean both the Licensee and the Licensor.
%
% These expressions may be used both in singular and plural form.
%
%     Article 2 - PURPOSE
%
% The purpose of the Agreement is the grant by the Licensor to the
% Licensee of a non-exclusive, transferable and worldwide license for the
% Software as set forth in Article 5 hereinafter for the whole term of the
% protection granted by the rights over said Software.
%
%     Article 3 - ACCEPTANCE
%
% 3.1 The Licensee shall be deemed as having accepted the terms and
% conditions of this Agreement upon the occurrence of the first of the
% following events:
%
%     * (i) loading the Software by any or all means, notably, by
%       downloading from a remote server, or by loading from a physical
%       medium;
%     * (ii) the first time the Licensee exercises any of the rights
%       granted hereunder.
%
% 3.2 One copy of the Agreement, containing a notice relating to the
% characteristics of the Software, to the limited warranty, and to the
% fact that its use is restricted to experienced users has been provided
% to the Licensee prior to its acceptance as set forth in Article 3.1
% hereinabove, and the Licensee hereby acknowledges that it has read and
% understood it.
%
%     Article 4 - EFFECTIVE DATE AND TERM
%
%       4.1 EFFECTIVE DATE
%
% The Agreement shall become effective on the date when it is accepted by
% the Licensee as set forth in Article 3.1.
%
%       4.2 TERM
%
% The Agreement shall remain in force for the entire legal term of
% protection of the economic rights over the Software.
%
%     Article 5 - SCOPE OF RIGHTS GRANTED
%
% The Licensor hereby grants to the Licensee, who accepts, the following
% rights over the Software for any or all use, and for the term of the
% Agreement, on the basis of the terms and conditions set forth hereinafter.
%
% Besides, if the Licensor owns or comes to own one or more patents
% protecting all or part of the functions of the Software or of its
% components, the Licensor undertakes not to enforce the rights granted by
% these patents against successive Licensees using, exploiting or
% modifying the Software. If these patents are transferred, the Licensor
% undertakes to have the transferees subscribe to the obligations set
% forth in this paragraph.
%
%       5.1 RIGHT OF USE
%
% The Licensee is authorized to use the Software, without any limitation
% as to its fields of application, with it being hereinafter specified
% that this comprises:
%
%    1. permanent or temporary reproduction of all or part of the Software
%       by any or all means and in any or all form.
%
%    2. loading, displaying, running, or storing the Software on any or
%       all medium.
%
%    3. entitlement to observe, study or test its operation so as to
%       determine the ideas and principles behind any or all constituent
%       elements of said Software. This shall apply when the Licensee
%       carries out any or all loading, displaying, running, transmission
%       or storage operation as regards the Software, that it is entitled
%       to carry out hereunder.
%
%       5.2 ENTITLEMENT TO MAKE CONTRIBUTIONS
%
% The right to make Contributions includes the right to translate, adapt,
% arrange, or make any or all modifications to the Software, and the right
% to reproduce the resulting software.
%
% The Licensee is authorized to make any or all Contributions to the
% Software provided that it includes an explicit notice that it is the
% author of said Contribution and indicates the date of the creation thereof.
%
%       5.3 RIGHT OF DISTRIBUTION
%
% In particular, the right of distribution includes the right to publish,
% transmit and communicate the Software to the general public on any or
% all medium, and by any or all means, and the right to market, either in
% consideration of a fee, or free of charge, one or more copies of the
% Software by any means.
%
% The Licensee is further authorized to distribute copies of the modified
% or unmodified Software to third parties according to the terms and
% conditions set forth hereinafter.
%
%         5.3.1 DISTRIBUTION OF SOFTWARE WITHOUT MODIFICATION
%
% The Licensee is authorized to distribute true copies of the Software in
% Source Code or Object Code form, provided that said distribution
% complies with all the provisions of the Agreement and is accompanied by:
%
%    1. a copy of the Agreement,
%
%    2. a notice relating to the limitation of both the Licensor's
%       warranty and liability as set forth in Articles 8 and 9,
%
% and that, in the event that only the Object Code of the Software is
% redistributed, the Licensee allows future Licensees unhindered access to
% the full Source Code of the Software by indicating how to access it, it
% being understood that the additional cost of acquiring the Source Code
% shall not exceed the cost of transferring the data.
%
%         5.3.2 DISTRIBUTION OF MODIFIED SOFTWARE
%
% When the Licensee makes a Contribution to the Software, the terms and
% conditions for the distribution of the resulting Modified Software
% become subject to all the provisions of this Agreement.
%
% The Licensee is authorized to distribute the Modified Software, in
% source code or object code form, provided that said distribution
% complies with all the provisions of the Agreement and is accompanied by:
%
%    1. a copy of the Agreement,
%
%    2. a notice relating to the limitation of both the Licensor's
%       warranty and liability as set forth in Articles 8 and 9,
%
% and that, in the event that only the object code of the Modified
% Software is redistributed, the Licensee allows future Licensees
% unhindered access to the full source code of the Modified Software by
% indicating how to access it, it being understood that the additional
% cost of acquiring the source code shall not exceed the cost of
% transferring the data.
%
%         5.3.3 DISTRIBUTION OF EXTERNAL MODULES
%
% When the Licensee has developed an External Module, the terms and
% conditions of this Agreement do not apply to said External Module, that
% may be distributed under a separate license agreement.
%
%         5.3.4 COMPATIBILITY WITH THE GNU GPL
%
% The Licensee can include a code that is subject to the provisions of one
% of the versions of the GNU GPL in the Modified or unmodified Software,
% and distribute that entire code under the terms of the same version of
% the GNU GPL.
%
% The Licensee can include the Modified or unmodified Software in a code
% that is subject to the provisions of one of the versions of the GNU GPL,
% and distribute that entire code under the terms of the same version of
% the GNU GPL.
%
%     Article 6 - INTELLECTUAL PROPERTY
%
%       6.1 OVER THE INITIAL SOFTWARE
%
% The Holder owns the economic rights over the Initial Software. Any or
% all use of the Initial Software is subject to compliance with the terms
% and conditions under which the Holder has elected to distribute its work
% and no one shall be entitled to modify the terms and conditions for the
% distribution of said Initial Software.
%
% The Holder undertakes that the Initial Software will remain ruled at
% least by this Agreement, for the duration set forth in Article 4.2.
%
%       6.2 OVER THE CONTRIBUTIONS
%
% The Licensee who develops a Contribution is the owner of the
% intellectual property rights over this Contribution as defined by
% applicable law.
%
%       6.3 OVER THE EXTERNAL MODULES
%
% The Licensee who develops an External Module is the owner of the
% intellectual property rights over this External Module as defined by
% applicable law and is free to choose the type of agreement that shall
% govern its distribution.
%
%       6.4 JOINT PROVISIONS
%
% The Licensee expressly undertakes:
%
%    1. not to remove, or modify, in any manner, the intellectual property
%       notices attached to the Software;
%
%    2. to reproduce said notices, in an identical manner, in the copies
%       of the Software modified or not.
%
% The Licensee undertakes not to directly or indirectly infringe the
% intellectual property rights of the Holder and/or Contributors on the
% Software and to take, where applicable, vis-a-vis its staff, any and all
% measures required to ensure respect of said intellectual property rights
% of the Holder and/or Contributors.
%
%     Article 7 - RELATED SERVICES
%
% 7.1 Under no circumstances shall the Agreement oblige the Licensor to
% provide technical assistance or maintenance services for the Software.
%
% However, the Licensor is entitled to offer this type of services. The
% terms and conditions of such technical assistance, and/or such
% maintenance, shall be set forth in a separate instrument. Only the
% Licensor offering said maintenance and/or technical assistance services
% shall incur liability therefor.
%
% 7.2 Similarly, any Licensor is entitled to offer to its licensees, under
% its sole responsibility, a warranty, that shall only be binding upon
% itself, for the redistribution of the Software and/or the Modified
% Software, under terms and conditions that it is free to decide. Said
% warranty, and the financial terms and conditions of its application,
% shall be subject of a separate instrument executed between the Licensor
% and the Licensee.
%
%     Article 8 - LIABILITY
%
% 8.1 Subject to the provisions of Article 8.2, the Licensee shall be
% entitled to claim compensation for any direct loss it may have suffered
% from the Software as a result of a fault on the part of the relevant
% Licensor, subject to providing evidence thereof.
%
% 8.2 The Licensor's liability is limited to the commitments made under
% this Agreement and shall not be incurred as a result of in particular:
% (i) loss due the Licensee's total or partial failure to fulfill its
% obligations, (ii) direct or consequential loss that is suffered by the
% Licensee due to the use or performance of the Software, and (iii) more
% generally, any consequential loss. In particular the Parties expressly
% agree that any or all pecuniary or business loss (i.e. loss of data,
% loss of profits, operating loss, loss of customers or orders,
% opportunity cost, any disturbance to business activities) or any or all
% legal proceedings instituted against the Licensee by a third party,
% shall constitute consequential loss and shall not provide entitlement to
% any or all compensation from the Licensor.
%
%     Article 9 - WARRANTY
%
% 9.1 The Licensee acknowledges that the scientific and technical
% state-of-the-art when the Software was distributed did not enable all
% possible uses to be tested and verified, nor for the presence of
% possible defects to be detected. In this respect, the Licensee's
% attention has been drawn to the risks associated with loading, using,
% modifying and/or developing and reproducing the Software which are
% reserved for experienced users.
%
% The Licensee shall be responsible for verifying, by any or all means,
% the suitability of the product for its requirements, its good working
% order, and for ensuring that it shall not cause damage to either persons
% or properties.
%
% 9.2 The Licensor hereby represents, in good faith, that it is entitled
% to grant all the rights over the Software (including in particular the
% rights set forth in Article 5).
%
% 9.3 The Licensee acknowledges that the Software is supplied "as is" by
% the Licensor without any other express or tacit warranty, other than
% that provided for in Article 9.2 and, in particular, without any warranty
% as to its commercial value, its secured, safe, innovative or relevant
% nature.
%
% Specifically, the Licensor does not warrant that the Software is free
% from any error, that it will operate without interruption, that it will
% be compatible with the Licensee's own equipment and software
% configuration, nor that it will meet the Licensee's requirements.
%
% 9.4 The Licensor does not either expressly or tacitly warrant that the
% Software does not infringe any third party intellectual property right
% relating to a patent, software or any other property right. Therefore,
% the Licensor disclaims any and all liability towards the Licensee
% arising out of any or all proceedings for infringement that may be
% instituted in respect of the use, modification and redistribution of the
% Software. Nevertheless, should such proceedings be instituted against
% the Licensee, the Licensor shall provide it with technical and legal
% assistance for its defense. Such technical and legal assistance shall be
% decided on a case-by-case basis between the relevant Licensor and the
% Licensee pursuant to a memorandum of understanding. The Licensor
% disclaims any and all liability as regards the Licensee's use of the
% name of the Software. No warranty is given as regards the existence of
% prior rights over the name of the Software or as regards the existence
% of a trademark.
%
%     Article 10 - TERMINATION
%
% 10.1 In the event of a breach by the Licensee of its obligations
% hereunder, the Licensor may automatically terminate this Agreement
% thirty (30) days after notice has been sent to the Licensee and has
% remained ineffective.
%
% 10.2 A Licensee whose Agreement is terminated shall no longer be
% authorized to use, modify or distribute the Software. However, any
% licenses that it may have granted prior to termination of the Agreement
% shall remain valid subject to their having been granted in compliance
% with the terms and conditions hereof.
%
%     Article 11 - MISCELLANEOUS
%
%       11.1 EXCUSABLE EVENTS
%
% Neither Party shall be liable for any or all delay, or failure to
% perform the Agreement, that may be attributable to an event of force
% majeure, an act of God or an outside cause, such as defective
% functioning or interruptions of the electricity or telecommunications
% networks, network paralysis following a virus attack, intervention by
% government authorities, natural disasters, water damage, earthquakes,
% fire, explosions, strikes and labor unrest, war, etc.
%
% 11.2 Any failure by either Party, on one or more occasions, to invoke
% one or more of the provisions hereof, shall under no circumstances be
% interpreted as being a waiver by the interested Party of its right to
% invoke said provision(s) subsequently.
%
% 11.3 The Agreement cancels and replaces any or all previous agreements,
% whether written or oral, between the Parties and having the same
% purpose, and constitutes the entirety of the agreement between said
% Parties concerning said purpose. No supplement or modification to the
% terms and conditions hereof shall be effective as between the Parties
% unless it is made in writing and signed by their duly authorized
% representatives.
%
% 11.4 In the event that one or more of the provisions hereof were to
% conflict with a current or future applicable act or legislative text,
% said act or legislative text shall prevail, and the Parties shall make
% the necessary amendments so as to comply with said act or legislative
% text. All other provisions shall remain effective. Similarly, invalidity
% of a provision of the Agreement, for any reason whatsoever, shall not
% cause the Agreement as a whole to be invalid.
%
%       11.5 LANGUAGE
%
% The Agreement is drafted in both French and English and both versions
% are deemed authentic.
%
%     Article 12 - NEW VERSIONS OF THE AGREEMENT
%
% 12.1 Any person is authorized to duplicate and distribute copies of this
% Agreement.
%
% 12.2 So as to ensure coherence, the wording of this Agreement is
% protected and may only be modified by the authors of the License, who
% reserve the right to periodically publish updates or new versions of the
% Agreement, each with a separate number. These subsequent versions may
% address new issues encountered by Free Software.
%
% 12.3 Any Software distributed under a given version of the Agreement may
% only be subsequently distributed under the same version of the Agreement
% or a subsequent version, subject to the provisions of Article 5.3.4.
%
%     Article 13 - GOVERNING LAW AND JURISDICTION
%
% 13.1 The Agreement is governed by French law. The Parties agree to
% endeavor to seek an amicable solution to any disagreements or disputes
% that may arise during the performance of the Agreement.
%
% 13.2 Failing an amicable solution within two (2) months as from their
% occurrence, and unless emergency proceedings are necessary, the
% disagreements or disputes shall be referred to the Paris Courts having
% jurisdiction, by the more diligent Party.
%
% Version 2.0 dated 2006-09-05.
