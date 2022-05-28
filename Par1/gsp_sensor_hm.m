
function [G, XCoords, YCoords] = gsp_sensor_hm(N, param)
if nargin < 2
    param = struct;
end

if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'N_try'), param.N_try = 10; end
if ~isfield(param, 'distribute'), param.distribute = 0; end
if ~isfield(param, 'connected'), param.connected = 1; end
if ~isfield(param, 'nnparam'), param.nnparam = {}; end
if ~isfield(param.nnparam, 'k'), param.nnparam.k = 4; end

for n=1:param.N_try
    [XCoords, YCoords] = create_coords_m(N,param.distribute);
    % sort rows for plotting reasons
    G = gsp_nn_graph(sortrows([XCoords, YCoords]),param.nnparam);
    if gsp_check_connectivity(G)
        break;
    elseif n == param.N_try
        fprintf('Warning! Graph is not connected\n');
    end
end

% Return the values

G.type = 'sensor';


G = gsp_graph_default_parameters(G);
end