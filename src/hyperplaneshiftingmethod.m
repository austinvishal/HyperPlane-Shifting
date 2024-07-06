% Algorithm based on
%Characterization of Parallel Manipulator AvailableWrench Set Facets,
%Gouttefarde M
function [C,d]=hyperplaneshiftingmethod(J,tmin,tmax,tol)
W=J;
% W=wrench_matrix(q); % here q is the configuration, unit vectors and cross prods,
%W is cable space Jacobian for multi link and for single body it is the
m=size(W,2); %no. of cables %linearly independent directions
n=size(W,1); %no. of dofs 
% m-n+1 remaining unit wrenches w_i if 3x6 n=3,m=6 n-1=2, m-n+1=4, 2
% linearly independent wrench and 4 remaining unit wrenches
% nbcomb= nchoosek(m,n-1); %number of combinations of m,n-1 ie choose n-1 linearly independent cable unit wrenches elements from set of m distinct objects
nbcomb=size(nchoosek(1:m,n-1),1);
I= nchoosek(1:m,n-1);
I_rem=nchoosek(1:m,m-n+1);
d=[];
% nbcombrem=size(nchoosek(1:m,n-1),1); % this will be same size as nbcomb in the first dimension, second
%dimension will be m-n+1
for i= 1:nbcomb
    I_0=I(i,:);
    V=W(:,I_0); %matrix formed by the n-1 columns (wi) of V
    c= null(V');% basis of the nullspace of V
    nullityOfc = size(c, 2);
    if nullityOfc==1
        C(:,2*i-1)=c;
        C(:,2*i)=-c;
        I_main=c'*W;
        I_positive=find(I_main>0); %instead of zero, we put the tol value also
        I_negative=find(I_main<0);
        d1=sum(tmax*I_main(I_positive))+sum(tmin*I_main(I_positive)); %d_positive
        d2=-sum(tmax*I_main(I_negative))-sum(tmin*I_main(I_positive));  %d_negative
    end
    d=[d,[d1,d2]];
    HreptoVrepduality(C',d') % for plotting 
end