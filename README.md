# matlab_seismic_cpml_iso_2d_curvil
seis_cpml_iso_2o_mtb_disp_curvil.m - disp formulation + curvilinear grid + fe_boundary - without multiple ifs
ksi - curvil coord x
eta - curvil coord y


curveintersect.m - func to find x0,y0 intersection of curves
func_find_closest_grid_nodes.m - finds nearest gridpoints, normals and so on
seis_cpml_iso_2o__curvil_matlab.m - main file

func_curv_jacob_pml - creates curvilinear grid. to test use the following expression
%[ksi,eta, gr_x,gr_y,J, Ji] = func_curv_jacob_pml(20,30,5,0,100, 0, 100,sphi,10,10,true)
%Constructs curvilinear mesh with regular regions for pmls and its Cartesian analogus. Calculates
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

