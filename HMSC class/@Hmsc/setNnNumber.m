function setNnNumber(m, nnNumber)

if length(nnNumber) ~= m.nr
	error('HMSC: Length of vector with number of knots  must be equal to the number of levels');
end;
m.nnNumber = nnNumber;

end