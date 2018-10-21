function [ hull ] = createHull( xy, hullMult )
%CREATEHULL Summary of this function goes here
%   Detailed explanation goes here
p = 1000;
K = convhull(xy(:,1),xy(:,2));
hull = xy(K,:);
hull = 1.03*(hull-repmat(mean(hull,1),size(hull,1),1)) + repmat(mean(hull,1),size(hull,1),1);
hull = insertRow(hull, p*[2460,2340], 14) ;
hull = insertRow(hull, p*[2760,2380], 17) ;
hull = insertRow(hull, p*[2270,2330], 13) ;
hull = insertRow(hull, p*[2890,2480], 21) ;
hull = insertRow(hull, p*[2788,2530], 22) ;
hull = insertRow(hull, p*[2720,2630], 2) ;
hull = insertRow(hull, p*[2460,2660], 3) ;
hull = insertRow(hull, p*[2380,2700], 4) ;
hull = insertRow(hull, p*[2350,2750], 5) ;
hull = insertRow(hull, p*[2310,2770], 6) ;
hull = insertRow(hull, p*[2280,2775], 7) ;
hull = insertRow(hull, p*[2248,2815], 8) ;

hull = hullMult/1.03*(hull-repmat(mean(hull,1),size(hull,1),1)) + repmat(mean(hull,1),size(hull,1),1);

end

