landareas = shaperead('landareas.shp','UseGeoCoords',true);
axesm ('eqdazim', 'Frame', 'on', 'Grid', 'on');
geoshow(landareas,'FaceColor',[1 1 .5],'EdgeColor',[.6 .6 .6]);
tissot;
