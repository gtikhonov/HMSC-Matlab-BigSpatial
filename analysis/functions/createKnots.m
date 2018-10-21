function [ xyKnot ] = createKnots( hull, knotN )
%CREATEKNONTS Summary of this function goes here
%   Detailed explanation goes here

bounds = [min(hull), max(hull)];
origin = (min(hull) + max(hull))/2;
edgeSise = sqrt(prod(max(hull)-min(hull))/knotN);
p = 1;
while p < 10
   edgeSiseL = edgeSise*(1-2^-p);
   edgeSiseU = edgeSise*(1+2^-p);
   xyKnot = hexagonalGrid(bounds, origin, edgeSise);
   xyKnotL = hexagonalGrid(bounds, origin, edgeSiseL);
   xyKnotU = hexagonalGrid(bounds, origin, edgeSiseU);
   in = inpolygon(xyKnot(:,1),xyKnot(:,2),hull(:,1),hull(:,2));
   inL = inpolygon(xyKnotL(:,1),xyKnotL(:,2),hull(:,1),hull(:,2));
   inU = inpolygon(xyKnotU(:,1),xyKnotU(:,2),hull(:,1),hull(:,2));
   kN = [sum(inL),sum(in),sum(inU)];
   [~,i] = min(abs(kN-knotN));
   switch i
      case 1
         edgeSise = edgeSiseL;
      case 2
         p = p+1;
      case 3
         edgeSise = edgeSiseU;
   end
end
xyKnot = hexagonalGrid(bounds, origin, edgeSise);
in = inpolygon(xyKnot(:,1),xyKnot(:,2),hull(:,1),hull(:,2));
xyKnot = xyKnot(in,:);

end

