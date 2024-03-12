function s = fill_area(t, xt1, xt2, color,txt)
t2 = [t, fliplr(t)];
inBetween = [xt1, fliplr(xt2)];
s = fill(t2, inBetween, color,'DisplayName',txt);
s.EdgeColor = color;
alpha(s,.3);
end