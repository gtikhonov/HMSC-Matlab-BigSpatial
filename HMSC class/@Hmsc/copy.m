function new = copy(m)
% Instantiate new object of the same class.
new = Hmsc(m.folder, m.traits, m.speciesX, m.phylogeny, m.includeXs, m.includeXv, m.outlierSpecies, m.spatial, m.factorCov);

% Copy all non-hidden properties.
p = properties(m);
for i = 1:length(p)
   new.(p{i}) = m.(p{i});
end
end