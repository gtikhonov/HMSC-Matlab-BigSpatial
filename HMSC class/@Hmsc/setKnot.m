function setKnot(m, xyKnot)

if length(xyKnot) ~= m.nr
	error('HMSC: Length of cell array with knots  must be equal to the number of levels');
end;
knotNumber = nan(1,m.nr);
for r = 1:m.nr
   if ~isempty(xyKnot{r})
      if m.spatial(r)
         knotNumber(r) = size(xyKnot{r},1);
      else
         error('HMSC: knots were provided for non spatial level %d', r);
      end
   end
m.knotNumber = knotNumber;
m.xyKnot = xyKnot;

end