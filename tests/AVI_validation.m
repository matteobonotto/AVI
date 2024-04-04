%%
close all;
clc;
clearvars
restoredefaultpath;

if ismac
    set(0,'DefaultFigureWindowStyle','docked')
end

addpath(genpath('../src/'))

%% Define geometry
case_ = 1;

switch case_
    case 1
        C = [5 0];
        load('../data/coil_geo_5m.mat')
    case 2
        C = [10 0];
        load('../data/coil_geo_10m.mat')
    case 3
        C = [20 0];
        load('../data/coil_geo_20m.mat')
end
delta_r = 0.05;
delta_z = 0.05;

coil = C + [
    -delta_r -delta_z;
    delta_r -delta_z;
    delta_r delta_z;
    -delta_r delta_z;
    ];


%% Mesh geometry
meshing = false;

if meshing
    clear GEOMETRY
    GEOMETRY.shape = {coil};
    GEOMETRY.spacing = {.01};
    GEOMETRY.surface = {1};
    filename='GMSH_GEOMETRY';
    order = 1;
    dirname = '../src/AutoMESH_VI_ver_1.0/';
    [p,e,t]=autoMESH_complex(filename,GEOMETRY,order,dirname);
    save('../data/coil_geo_10m.mat', 'p', 'e', 't')
% else
%     load('../data/coil_geo_5m.mat')
% %     load('../data/coil_geo_20m.mat')
end

%% Build meshData_pas
N_order= 1;

meshData_pas.t = t;
meshData_pas.n = p;
meshData_pas.type = ones(size(meshData_pas.t,1),1);
meshData_pas.nn = size(meshData_pas.n,1);
meshData_pas.nt = size(meshData_pas.t,1);
meshData_pas.N_order = N_order;

tri = meshData_pas.t;
nodes = meshData_pas.n;
type_pas = meshData_pas.type;

nodes_matrix = [nodes(tri(:,1),:) ...
    nodes(tri(:,2),:) ...
    nodes(tri(:,3),:)];

P1 = nodes_matrix(:,[1 2]);
P2 = nodes_matrix(:,[3 4]);
P3 = nodes_matrix(:,[5 6]);

degree_G_target = 8;
degree_G_source = 12;
[w_G_target,P_G_target,n_G_target] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_target);
[w_G_source,P_G_source,n_G_source] = fun_Gauss_points_triangle_Dunavant(P1,P2,P3,degree_G_source);

meshData_pas.n_G_source = n_G_source;
meshData_pas.P_G_source = P_G_source;
meshData_pas.w_G_source = w_G_source;

[MatInd_nodes_tri] = fun_MatInd_nodes_tri(meshData_pas.N_order,meshData_pas.t,meshData_pas.nn);

degree_G_source_ext = degree_G_target;


%% Compute AVI matrices
rho = 1.68e-5; % copper
eta = rho*ones(meshData_pas.nt,1);
RUN_MEX = 0;


tic
L_VI = fun_L_VI_stable_symmetric_fast(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    degree_G_target, ...
    n_G_target, ....
    P_G_target, ...
    RUN_MEX);
toc

tic
R_VI = fun_R_VI_stable(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    eta);
toc

tic
U_VI = fun_U_VI_stable_fast(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source);
toc

tic
V_VI = fun_V_VI_stable_fast(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source);
toc


tic
ind_t_map = [1 1];
D_VI = fun_D_VI(tri,...
    type_pas, ...
    ind_t_map);
toc


% % tic
% % [M_VI,~] = fun_M_act_VI_stable_fast(tri, ...
% %     nodes, ...
% %     N_order, ...
% %     degree_G_source_ext, ...
% %     degree_G_target, ...
% %     n_G_target, ....
% %     P_G_target,...
% %     n_act, ...
% %     tri_act, ...
% %     nodes_act, ...
% %     ind_act, ...
% %     keyreg_act);
% % toc

R_fluxloop = mean(coil(:,1));
Z_fluxloop = mean(coil(:,2));

%%% Passive contribution
tic
G_flux_VI_pas = fun_G_flux_passive_VI_stable_fast(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source, ...
    n_G_source, ...
    P_G_source,...
    R_fluxloop, ...
    Z_fluxloop, ...
    RUN_MEX);
toc

ind_pas = 1;
vec_A = fun_vec_Area_element_VI_stable(tri, ...
    nodes, ...
    N_order, ...
    MatInd_nodes_tri,...
    degree_G_source, ...
    n_G_source, ...
    ind_pas, ...
    meshData_pas.type);


%%% for source coil use the passive coils shifted outwards by 1m 
n_act = 1;
tri_act = tri;
nodes_act = nodes + [10 0];
ind_act = 1;
keyreg_act = ones(size(tri_act(:,1)));

tic
[M_VI,~] = fun_M_act_VI_stable_fast(tri, ...
    nodes, ...
    N_order, ...
    degree_G_source_ext, ...
    degree_G_target, ...
    n_G_target, ....
    P_G_target,...
    n_act, ...
    tri_act, ...
    nodes_act, ...
    ind_act, ...
    keyreg_act, ...
    RUN_MEX);
toc

tic
G_flux_VI_act =fun_G_flux_active_VI_stable_fast(n_act, ...
    tri_act, ...
    nodes_act, ...
    ind_act, ...
    keyreg_act, ...
    degree_G_source, ...
    R_fluxloop, ...
    Z_fluxloop, ...
    RUN_MEX);
toc

figure; hold on, axis equal;
triplot(t,p(:,1),p(:,2));
triplot(tri_act,nodes_act(:,1),nodes_act(:,2));
plot(R_fluxloop,Z_fluxloop,'o')
% xlim([0,1.3])

%% State-space matrices
L_VI_inv = L_VI\eye(size(L_VI));

% A_SS_ODE = -R_VI*L_VI_inv;
% A = A_SS_ODE;
% B = -A_SS_ODE;
% C = L_VI_inv;
% D = -L_VI_inv;


%%
I0 = 1e+4;
f=50;
omega = 2*pi*f;
time = linspace(0,15*2*pi/omega,10000);
R0 = mean(coil(:,1));
source = I0*sin(omega*time);

addpath '/Users/matte/Library/CloudStorage/Dropbox/PhD/RESEARCH_ACTIVITY/Inductance_Coefficient/fun_Inductance_ver_1.01'

L_tot = fun_M_Green_Axi_Quad_Quad(coil,coil,'F','F');
R_tot = fun_R_quad(coil,rho);

%%% Reference solution
% sinusoidal current on the source coil
R_tot = rho*2*pi*R0/(2*delta_r*2*delta_z);
L_tot = G_flux_VI_pas*(ones(size(vec_A)).'/sum(vec_A)).'; %to get total unit current = 1
% L_tot = fun_M_Green_Axi_Quad_Quad(coil,coil,'F','F');
M_tot = G_flux_VI_act;

% quad_1 = coil;
% quad_1 = fun_ordinapunti(quad_1);
% Area_coil = polyarea(quad_1(:,1),quad_1(:,2));
% center_coil = sum(quad_1)/size(quad_1,1);
% R = rho_mat*2*pi*center_coil(1)/Area_coil;

A = -R_tot/L_tot;
B = R_tot/L_tot*M_tot;
C = 1/L_tot;
D = -M_tot/L_tot;

x_t0 = -A\B*source(:,1);
tic
phi = fun_ODE_Crank_Nicolson_State_Space(time,x_t0,source,A,B);
Ic = C*phi + D*source;
toc

Vc = R_tot*Ic;



%%% VI solution
% sinusoidal current on the source coil
L_VI_inv = L_VI\eye(size(L_VI));


A = -R_VI*L_VI_inv;
B = R_VI*L_VI_inv*M_VI;
C = L_VI_inv;
D = -L_VI_inv*M_VI;

x_t0 = -A\B*source(:,1);
tic
aa = fun_ODE_Crank_Nicolson_State_Space(time,x_t0,source,A,B);
jc = C*aa + D*source;
toc

Ic_VI = jc.'*vec_A;

max(Ic_VI)/max(Ic)
max(Ic)/max(Ic_VI)

Vc_VI = mean(((inv(U_VI)*R_VI*2*pi.*meshData_pas.n(:,1))*jc),1);

Vc_VI_lumped = R_tot*Ic_VI;

ind = floor(1000*rand(1));
figure(ind)
subplot(2,1,1); hold on
plot(time,Ic, 'color', 'k')
subplot(2,1,2)
plot(time,Vc, 'color', 'k'); hold on


figure(ind)
subplot(2,1,1)
plot(time(1:10:end),Ic_VI(1:10:end),'o')
subplot(2,1,2)
plot(time(1:10:end),Vc_VI(1:10:end),'ro')
plot(time(1:10:end),Vc_VI_lumped(1:10:end),'go')













here = 1;





















