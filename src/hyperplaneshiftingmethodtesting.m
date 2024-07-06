% Algorithm based on
%Characterization of Parallel Manipulator AvailableWrench Set Facets,
%Gouttefarde M
% function [outs]=hyperplaneshiftingmethod(q,tmin,tmax)
% output C,d
% C(list): matrix of half-space representation `Cx<d`
% d(list): vector of half-space representation `Cx<d`
   
clear all
tmin=2;
tmax=120;
% tmin=[2 2 2 2 2 2];
% tmax=[120 120 120 120 120 ];
Wrench_shomain =[ 0.1046   -0.0231    0.0050   -0.0157   -0.0555    0.0201
    0.0045   -0.0197   -0.0403   -0.0087    0.0043   -0.0064
   -0.0386    0.0464   -0.0311   -0.0020    0.0151   -0.0048];
% Wrench_shomain=[0.0532    0.0532   -0.0266   -0.0266   -0.0266   -0.0266
%          0         0    0.0460    0.0460   -0.0460   -0.0460
%    -0.0442    0.0442   -0.0442    0.0442   -0.0442    0.0442]; %single shoulder
% Wrench_shomain=[0.0393   -0.0000    0.0327  %doubleshoulder
%     0.0393    0.0000   -0.0327
%    -0.0197    0.0341    0.0327
%    -0.0197    0.0341   -0.0327
%    -0.0197   -0.0341    0.0327
%    -0.0197   -0.0341   -0.0327
%    -0.0473   -0.0192    0.0485
%    -0.0473    0.0192   -0.0485
%     0.0402   -0.0313    0.0485
%     0.0070   -0.0505   -0.0485
%     0.0070    0.0505    0.0485
%     0.0402    0.0313   -0.0485]';
%   Wrench_shomain = 1.0e-04 *[ -0.2457   -0.2457   -0.2444    0.4901    0.4901   -0.2444
%     0.4237   -0.4237   -0.4245   -0.0008    0.0008    0.4245
%     0.0008   -0.0008    0.0008   -0.0008    0.0008   -0.0008];
% Wrench_shomain= [-0.0470    0.0552   -0.0485   -0.0463    0.0783    0.0302
%     0.0343    0.0158   -0.0057    0.0258    0.0141   -0.0083
%    -0.0137    0.0029   -0.0128    0.0265   -0.0386    0.0319];
% Wrench_shomain =[  0         0         0         0         0         0
%          0         0    0.0460    0.0460   -0.0460   -0.0460
%    -0.0442    0.0442   -0.0442    0.0442   -0.0442    0.0442]; % singular case, rank drop, failign hsm
%            Wrench_shomain=[0 1 1 0 0 1 1 0
%                            0 0 0 1 0 1 1 1
%                            0 0 1 1 1 1 0 0];
W=Wrench_shomain;
% W=wrench_matrix(q); % here q is the configuration, unit vectors and cross prods,
%W is cable space Jacobian for multi link and for single body it is the
%Jacobian
m=size(W,2); %no. of cables %linearly independent directions
n=size(W,1); %no. of dofs 
% m-n+1 remaining unit wrenches w_i if 3x6 n=3,m=6 n-1=2, m-n+1=4, 2
% linearly independent wrench and 4 remaining unit wrenches
% nbcomb= nchoosek(m,n-1); %number of combinations of m,n-1 ie choose n-1 linearly independent cable unit wrenches elements from set of m distinct objects
% I= nchoosek(W,nbcomb);%list of combinations of (m n-1) indices
% C=zeros(2*nbcomb,n);
% d=zeros(2*nbcomb,1);
% What = ObtainLinearlyIndependentSet(W);
%         w_t             =   -What*ones(n,1); % A pose is WCW if the positive span contains this vector
%         combinations    =   nchoosek(1:m,n);
%         nbcomb2= nchoosek(1:m,n-1);
%         size(nbcomb2,1)
% G = A(:,combinations(i,:));
% nbcomb2= nchoosek(1:m,n-1);
nbcomb=size(nchoosek(1:m,n-1),1);
%  C=zeros(2*nbcomb,n);
%   d=zeros(2*nbcomb,1);
I= nchoosek(1:m,n-1);
I_rem=nchoosek(1:m,m-n+1);
d=[];
dx=[];
% nbcombrem=size(nchoosek(1:m,n-1),1); % this will be same size as nbcomb in the first dimension, second
%dimension will be m-n+1
% for i= 1:nbcomb-1
for i= 1:nbcomb
    I_0=I(i,:);
    V=W(:,I_0); %matrix formed by the n-1 columns (wi) of V
%     c=cross_product_normalized(V(:,2),V(:,1))
     c= null(V');% basis of the nullspace of V
    nullityOfc = size(c, 2);
    if nullityOfc==1
%         I1_pos= 
%         C(:,2*i)=c;
%         C(:,2*i+1)=-c;
        C(:,2*i-1)=c;
         C(:,2*i)=-c;
%          I_main=c'*W(:,I_rem(i,:))
          I_main=c'*W  %either take remaining wrenches(m-n+1) or take full projection matrix(which might end up with zero values 
%          ... after multiplication with orthogonal vector)
%if we take full matrix, the 0 value comes as 0.0000 which actually is a small value and not fully zero
% so find a better way to eliminate this decimal 0 value as it causes different array sizes
% as of now remaining wrenches is safer option and works fine
% we need to keep an eps/tol value tol = 1e-15, eps epsilon relative
% accuracy of floating point number
%           I_positive=find(I_main>1e-15); %instead of zero, we put the tol value
%           I_negative=find(I_main<-1e-15);
            I_positive=find(I_main>0); %instead of zero, we put the tol value
          I_negative=find(I_main<0);
%          for j=1:m
%              if I_main>0.0  %note equality is not checked where I is closed halfspace bounded by Hi
%                 I_positive=j;
%              elseif I_main<0.0
%                 I_negative=j;
%              else
%                  disp('I_main is null')
%              end
         d1=sum(tmax*I_main(I_positive))+sum(tmin*I_main(I_positive)) %d_positive
         d2=-sum(tmax*I_main(I_negative))-sum(tmin*I_main(I_positive))  %d_negative
         end
%         d(2*i)=d1;
%         d(2*i+1)=d2;
%            d(:,2*i-1)=d1;
%          d(:,2*i)=d2;
         d=[d,[d1,d2]];
         
%          dx=[dx,d]
end
    % A=Cx<d in H-rep the zonotope shaped available wrench set A is
    % expressed
%  end
% f= %wrench induced by the weight of the platform at configuration q
% dq= C*f;
% for i=0: 2*nbcomb
%    if dq(i)>d(i)
%        outs = False;      
%    end
% outs = true; 
% % return 1
% end
% end
%% computing vertices V from H rep is vertex enumeration problem
%% computing H from V rep is facet enumeration problem
% HreptoVrepduality(C',d')
%% ref: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/48509/versions/3/previews/COMP_GEOM_TLBX/html/Vertex_enumeration_3D.html
% a convex polytope can be specified as convex hull  of vertex set V of P or
% intersection of set H of its facet including half spaces
% using duality principle convert the planes to equivalent points
% https://stackoverflow.com/questions/26809630/how-to-convert-the-half-spaces-that-constitute-a-convex-hull-to-a-set-of-extreme
A= C';
b=d';
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

% axis equal
% axis vis3d
% axis off
% h=camlight(0,90);