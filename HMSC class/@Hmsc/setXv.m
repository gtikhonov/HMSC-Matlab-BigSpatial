function setXv(m, Xv)

% check for correctness


if m.includeXv == 0 && ~isempty(Xv)
	error('HMSC: model was defined without dimension reduction part');
end


if m.includeXv == 2
	if length(Xv)~=m.ns
		error('HMSC: Length of Xv does not match mumber of species, defined by Y');
	end
	m.ncv = size(Xv{1}, 2);
	for i = 1:m.ns
		if any(size(Xv{i}) ~= [m.ny, m.ncv])
			error(paste('HMSC: Xv{', int2str(i), '} dimensions are not compatible with specified'));
		end
	end
else
	m.ncv = size(Xv, 2);
	if any(size(Xv) ~= [m.ny, m.ncv])
		error('HMSC: Xv dimensions are not compatible with specified');
	end
end
m.Xv = Xv;
m.vCovScaleFlag = zeros(1, m.ncv);
m.vCovScale = zeros(2, m.ncv);

end