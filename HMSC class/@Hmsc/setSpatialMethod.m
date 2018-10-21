function setSpatialMethod(m, spatialMethod)

if length(spatialMethod) ~= m.nr
	error('HMSC: Length of vector must be equal to the number of levels');
end;
m.spatialMethod = spatialMethod;

end