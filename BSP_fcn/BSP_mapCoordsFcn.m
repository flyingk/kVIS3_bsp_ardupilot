function [lat, lon, alt, t] = BSP_mapCoordsFcn(fds)

group = 'GPS';

t = kVIS_fdsGetChannel(fds, group, 'Time');

lon = kVIS_fdsGetChannel(fds, group, 'Lng');

lat = kVIS_fdsGetChannel(fds, group, 'Lat');

alt = kVIS_fdsGetChannel(fds, group, 'Alt');
alt(alt < 0) = 0;

end