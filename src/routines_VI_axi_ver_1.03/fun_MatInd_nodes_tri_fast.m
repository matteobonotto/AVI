function [MatInd_nodes_tri] = fun_MatInd_nodes_tri_fast(N_order,tri,nn) %#codegen

MatInd_nodes_tri = [];

if N_order == 1
    
    tri1 = tri(:,1);
    tri2 = tri(:,2);
    tri3 = tri(:,3);
    
    num_el_max = 100;
    
    num_el = zeros(nn,1);
    MatInd_nodes_tri = zeros(nn,num_el_max);
    parfor ii=1:nn
        qq = unique([find(tri1 == ii); ...
            find(tri2 == ii); ...
            find(tri3 == ii)]);
        
        tmp = numel(qq);
        num_el(ii) = tmp;
        
        ind = 1:numel(qq);
        tmp = zeros(1,num_el_max);
        tmp(ind) = qq;
        MatInd_nodes_tri(ii,:) = tmp;
        
    end
    
    num_el_max = max(num_el);
    MatInd_nodes_tri = MatInd_nodes_tri(:,1:num_el_max);
    
    
elseif N_order == 2
    
    tri1 = tri(:,1);
    tri2 = tri(:,2);
    tri3 = tri(:,3);
    tri4 = tri(:,4);
    tri5 = tri(:,5);
    tri6 = tri(:,6);
    
    num_el_max = 100;
    
    num_el = zeros(nn,1);
    MatInd_nodes_tri = zeros(nn,num_el_max);
    parfor ii=1:nn
        qq = unique([find(tri1 == ii); ...
            find(tri2 == ii); ...
            find(tri3 == ii); ...
            find(tri4 == ii); ...
            find(tri5 == ii); ...
            find(tri6 == ii)]);
        
        tmp = numel(qq);
        num_el(ii) = tmp;
        
        ind = 1:numel(qq);
        tmp = zeros(1,num_el_max);
        tmp(ind) = qq;
        MatInd_nodes_tri(ii,:) = tmp;
        
    end
    
    num_el_max = max(num_el);
    MatInd_nodes_tri = MatInd_nodes_tri(:,1:num_el_max);
    
end

