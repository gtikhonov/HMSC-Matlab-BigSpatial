function setKnotNumber(m, knotNumber)

if length(knotNumber) ~= m.nr
	error('HMSC: Length of vector with number of knots  must be equal to the number of levels');
end;
m.knotNumber = knotNumber;
m.xyKnot = cell(1,m.nr);

end