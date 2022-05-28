
function [XCoords, YCoords] = create_coords_m(N, distribute)


N = (N/4);
dis = 2;


XCoords1 = rand(N,1);
YCoords1 = rand(N,1)+dis;

XCoords2 = rand(N,1);
YCoords2 = rand(N,1);

XCoords3 = rand(N,1)+dis;
YCoords3 = rand(N,1);

XCoords4 = rand(N,1)+dis;
YCoords4 = rand(N,1)+dis;

XCoords = [XCoords1;XCoords2;XCoords3;XCoords4];
YCoords = [YCoords1;YCoords2;YCoords3;YCoords4];

end