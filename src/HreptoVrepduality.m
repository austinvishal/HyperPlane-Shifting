% source: % https://stackoverflow.com/questions/26809630/how-to-convert-the-half-spaces-that-constitute-a-convex-hull-to-a-set-of-extreme
function HreptoVrepduality(A,b)
% A= C';
% b=d';
hv=A\b %least square solution, backlash invokes linear solver
b=b-A*hv 
D=A./repmat(b,[1 size(A,2)])
  [k,vol]=convhulln(D);
% trisurf(k,D(:,1),D(:,2),D(:,3),'FaceColor','cyan')
for ix = 1:size(k,1)
    F = D(k(ix,:),:);
    G(ix,:)=F\ones(size(F,1),1); %least square solution, backlash invokes linear solver
end
V = G + repmat(hv',[size(G,1),1]);
[kA,vol]=convhulln(V);
trisurf(kA,V(:,1),V(:,2),V(:,3),'FaceColor','none')

end