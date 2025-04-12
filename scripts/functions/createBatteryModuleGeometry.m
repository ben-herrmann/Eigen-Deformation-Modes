function [geomModule, domainIDs, boundaryIDs, volume, boundaryArea, ReferencePoint] = createBatteryModuleGeometry(numCellsInModule, cellWidth,cellThickness,tabThickness,tabWidth,cellHeight,tabHeight, connectorHeight )
% Make 2-D geometry to extured
[geom, ReferencePoint] = getCellRectangels(cellWidth, cellThickness, tabThickness, tabWidth);

% Extrude cell body and merge internal domains.
geom1=fegeometry(geom);
geom2=extrude (geom1,cellHeight);

% Merge unwanted internal cells
geom2 = geom2.mergeCells([1,2]);
geom2 = geom2.mergeCells([1,2]);

% Extrude tab faces
tabTopFaces  = geom2.nearestFace([ReferencePoint.TabLeft;ReferencePoint.TabRight] + [0,0,cellHeight]);
geom3=extrude(geom2,tabTopFaces,[tabHeight,connectorHeight]); % Cell 1, 2 and 3

%Extrude connector faces
connectorFaces1 = geom3.nearestFace(ReferencePoint.TabLeft+[ tabThickness/2,0,cellHeight + tabHeight + connectorHeight/2]);
connectorFaces2 = geom3.nearestFace(ReferencePoint.TabRight+[-tabThickness/2,0,cellHeight + tabHeight + connectorHeight/2]);
geom4=extrude(geom3,connectorFaces1,(cellThickness-tabThickness)/2);
geom5=extrude(geom4,connectorFaces2,(cellThickness-tabThickness)/2);

ReferencePoint.Cell = ReferencePoint.Cell + [0,0,cellHeight/2];
ReferencePoint.TabLeft = ReferencePoint.TabLeft + [0,0,cellHeight+tabHeight/2];
ReferencePoint.ConnectorLeft = ReferencePoint.TabLeft + [0,0,connectorHeight];
ReferencePoint.ConnectorLeftExtension = ReferencePoint.ConnectorLeft + [(cellThickness-tabThickness)/4+tabThickness/2,0,0];
ReferencePoint.TabRight = ReferencePoint.TabRight + [0,0,cellHeight+tabHeight/2];
ReferencePoint.ConnectorRight = ReferencePoint.TabRight + [0,0,connectorHeight];
ReferencePoint.ConnectorRightExtension = ReferencePoint.ConnectorRight + [-(cellThickness-tabThickness)/4-tabThickness/2,0,0];

geom5=generateMesh(geom5,'GeometricOrder','linear');
domains.Cell = findCellID(geom5,ReferencePoint.Cell);
domains.TabLeft = findCellID(geom5,ReferencePoint.TabLeft);
domains.TabRight = findCellID(geom5,ReferencePoint.TabRight);
domains.ConnectorLeft = findCellID(geom5,[ReferencePoint.ConnectorLeft; ReferencePoint.ConnectorLeftExtension]);
domains.ConnectorRight = findCellID(geom5,[ReferencePoint.ConnectorRight; ReferencePoint.ConnectorRightExtension]);

% Construct one cell geometry using mesh.
elems = geom5.Mesh.Elements(1:4,:)';
nodes = geom5.Mesh.Nodes';
elemToRegion = zeros(size(elems,1),1);
c = 1;
renumberedDomainIDs = zeros(5,1); % 1 cell, 2 tabs and two connector ends per cell.
for i = fields(domains)'
    cellElems = findElements(geom5.Mesh,'region','Cell',domains.(i{1}));
    elemToRegion(cellElems) = c;
    renumberedDomainIDs(c) = c;
    c = c+1;
end

nodesModue = nodes;
elemsModule = elems;
elemToRegionModule = elemToRegion;

domainIDs(1:numCellsInModule) = struct("Cell",1, ...
    "TabLeft",2,"TabRight",3, ...
    "ConnectorLeft",4,"ConnectorRight",5);

numNodes = size(nodes,1);
numElems = size(elems,1);

% Add the remaining copies of cell geometry. Mirror alternate cells so that
% connectors are in touch such that the whole pack is in series connection.
for n = 1:numCellsInModule-1
    newDominIDs = renumberedDomainIDs+(numel(renumberedDomainIDs))*n;
    nodesToAppend = nodes;
    elemsAppend = elems + size(nodesModue,1); % offset indices
    if mod(n,2) ~= 0
        nodesToAppend(:,1) = -1*nodesToAppend(:,1); % mirror
        elemsAppend = elemsAppend(:,[2 1 3 4]); % fix node-order
    else
        maxX = max(nodesToAppend(:,1));
        nodesToAppend(:,1) = nodesToAppend(:,1)-maxX; % translate
    end
    nodesModue(end+1:end+numNodes,:) = nodesToAppend;
    elemsModule(end+1:end+numElems,:) = elemsAppend;
    elemToRegionModule(end+1:end+numElems) = newDominIDs(elemToRegion);
    minX = min(nodesModue(:,1));
    nodesModue(:,1) = nodesModue (:,1)- minX;
    
    cellIDi = struct("Cell",newDominIDs(1),...
        "TabLeft",newDominIDs(2),"TabRight",newDominIDs(3),...
        "ConnectorLeft",newDominIDs(4),"ConnectorRight",newDominIDs(5));
    domainIDs(n+1) = cellIDi;
end

% Remove duplicate nodes
[nodesModue,~,ic] = uniquetol(nodesModue,1e-14,'ByRows',true);
elemsModule = ic(elemsModule);

% Construct geometry from mesh for the whole module
model = createpde;
geomModule = geometryFromMesh(model,nodesModue',elemsModule',elemToRegionModule);

% Offset x-coordinate of reference points the the center of the first cell
% of the pack.
xOffset = max(geomModule.Vertices(:,1))-cellThickness/2;
for f = fields(ReferencePoint)'
ReferencePoint.(f{1})(1) = xOffset;
end

% Helper function handles
getFaceID = @(offsetVal,offsetDirection,cellNumber) nearestFace(geomModule,...
    ReferencePoint.Cell + offsetVal/2 .*offsetDirection ... % offset ref. point to face
    - cellThickness*(cellNumber-1)*[1,0,0]);        % offset to cell

getVolume = @(geomCellID) model.Mesh.volume(findElements(model.Mesh,'region',"Cell",geomCellID));

% Initilize strcut arrays for boundary IDs, boundary area, and volume of
% cells.
boundaryIDs(1:numCellsInModule) = struct("FrontFace",[],"BackFace",[], ...
    "RightFace",[],"LeftFace",[], ...
    "TopFace",[],"BottomFace",[]);
boundaryArea(1:numCellsInModule) = struct("FrontFace",[],"BackFace",[], ...
    "RightFace",[],"LeftFace",[], ...
    "TopFace",[],"BottomFace",[]);
volume(1:numCellsInModule) = struct("Cell",0, ...
            "TabLeft",0,"TabRight",0, ...
            "ConnectorLeft",0,"ConnectorRight",0);

for n = 1:numCellsInModule
    boundaryIDs(n).FrontFace =  getFaceID(cellThickness, [ 1, 0, 0],n);
    boundaryIDs(n).BackFace =   getFaceID(cellThickness, [-1, 0, 0],n);
    boundaryIDs(n).RightFace =  getFaceID(cellWidth,     [ 0, 1, 0],n);
    boundaryIDs(n).LeftFace =   getFaceID(cellWidth,     [ 0,-1, 0],n);
    boundaryIDs(n).TopFace =    getFaceID(cellHeight,    [ 0, 0, 1],n);
    boundaryIDs(n).BottomFace = getFaceID(cellHeight,    [ 0, 0,-1],n);
    
    boundaryArea(n).FrontFace = cellHeight*cellWidth;
    boundaryArea(n).BackFace =  cellHeight*cellWidth;
    boundaryArea(n).RightFace = cellThickness*cellHeight;
    boundaryArea(n).LeftFace = cellThickness*cellHeight;
    boundaryArea(n).TopFace =  cellThickness*cellWidth - tabThickness*tabWidth;
    boundaryArea(n).BottomFace = cellThickness*cellWidth - tabThickness*tabWidth;

    volume(n).Cell =  getVolume(domainIDs(n).Cell);
    volume(n).TabLeft =  getVolume(domainIDs(n).TabLeft);
    volume(n).TabRight =  getVolume(domainIDs(n).TabRight);
    volume(n).ConnectorLeft =  getVolume(domainIDs(n).ConnectorLeft);
    volume(n).ConnectorRight =  getVolume(domainIDs(n).ConnectorRight);
end

end

% Helper functions
function [gm, ReferencePoint] = getCellRectangels(width, thickness, tabThickness, tabWidth)
xCoords = [0,thickness,thickness,0];
yCoords = [0,0,width,width];

xCoordsT = [0,tabThickness,tabThickness,0];
yCoordsT = [0,0,tabWidth,tabWidth];

Rcell = [3, 4, xCoords, yCoords]';
ReferencePoint.Cell = [thickness/2,width/2,0];

RT1 = [3,4, xCoordsT + (thickness-tabThickness)/2,yCoordsT +  tabWidth/2]';
RT2 = [3,4, xCoordsT + (thickness-tabThickness)/2,yCoordsT +  width - tabWidth*3/2]';
ReferencePoint.TabLeft = [thickness/2,tabWidth,0];
ReferencePoint.TabRight = [thickness/2,width - tabWidth,0];

gd = [Rcell, RT1, RT2];
sf = join("R"+ (1:3),"+");
ns = "R" + (1:3)';
gm = decsg(gd,sf,ns');
end

% Helper function
function fID = findCellID(gm,pt)
% A utility function to identify the cell ID given a point inside the
% geometric cell.

msh = gm.Mesh;
doms = cell(1,gm.NumCells);

for k=1:gm.NumCells
    doms{k} = findElements(msh,'region','Cell',k);
end
% Convert mesh to triangulation and find containing element
tr = triangulation(msh.Elements',msh.Nodes');
fID = zeros(1,size(pt,1));
for i = 1:size(pt,1)
    zz = pointLocation(tr,pt(i,:));
    mem = cellfun(@(x)ismember(zz,x),doms);

    fID(i) = find(mem,1,'first');
end
if isempty(fID)
    error("Outside of all domains");
end
end