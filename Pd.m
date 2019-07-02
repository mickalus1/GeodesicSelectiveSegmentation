function Pd1 = Pd( mask )
%Pd calculates the distance matrix Pd as in the Spencer-Chen 2015 paper
%   Based on user input polygon this function determines the diatance from
%   each pixel to the selected polygon.

Pd=bwdist(mask);
Pd1=Pd/max(Pd(:));
Pd1=double(Pd1);

%imagesc(Pd1);colormap(gray);

end

