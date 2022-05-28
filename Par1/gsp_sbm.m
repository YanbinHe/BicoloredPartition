function [G] = gsp_sbm(N, P,XCoords, YCoords, param)
if nargin < 5
    param = struct;
end

if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'N_try'), param.N_try = 10; end
if ~isfield(param, 'distribute'), param.distribute = 0; end
if ~isfield(param, 'connected'), param.connected = 1; end
if ~isfield(param, 'nnparam'), param.nnparam = {}; end
if ~isfield(param.nnparam, 'k'), param.nnparam.k = 4; end


% generate adjacency matrix W
% intra_clu = 0.9;
% inter_clu = 0.2;
% P = [intra_clu inter_clu 0 inter_clu-0.1;
%     inter_clu intra_clu inter_clu-0.1 0;
%     0 inter_clu-0.1 intra_clu inter_clu;
%     inter_clu-0.1 0 inter_clu intra_clu];
    
c = [ones(N/4,1);2*ones(N/4,1);3*ones(N/4,1);4*ones(N/4,1)];
W = generateSbm(c,P);


for n=1:param.N_try
%     [XCoords, YCoords] = create_coords_m(N,param.distribute);
    % sort rows for plotting reasons
    G = gsp_graph(W, ([XCoords, YCoords]),param.nnparam);
    if gsp_check_connectivity(G)
        break;
    elseif n == param.N_try
        fprintf('Warning! Graph is not connected\n');
    end
end

G.type = 'sensor';


G = gsp_graph_default_parameters(G);
end