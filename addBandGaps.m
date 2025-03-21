function addBandGaps(bg,U)
val = [min(abs(U))*0.9 max(abs(U))*1.1];
for i1 = 1:size(bg,1)
    hold on; patch( [bg(i1,1) bg(i1,1) bg(i1,2) bg(i1,2)], ...
                    [val(1) val(end) val(end) val(1)], ...                 
                    .8*[1 1 1],'EdgeColor','none' ); alpha(0.3)
end
