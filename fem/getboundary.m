function [geo,geosquare]=getboundary(Pb,geo,lsquare)

geosquare=geo;

if lsquare==0
    v=zeros(size(Pb,2),1);
    tol=10^(-10);
    for k=1:length(v)
        if(abs(Pb(2,k))<tol || abs(Pb(2,k)-1)<tol || abs(Pb(1,k))<tol || (abs(Pb(1,k)-1)<tol)) 
            v(k)=-1;
        end
    end
    geo.int=find(v==0);
    geo.out=find(v==-1);
end


if  lsquare==1
    geosquare=geo;
    
    Nh=sqrt(size(geo.Pb,2))-2;
    v=ones(size(Pb,2),1); %get matrix Pb containining also boundary nodes for lsquare
    vec_new_to_old=zeros(size(Pb,2),1); %vec(jj)=k is the global_index k of node jj in the modified mesh.
    vec_old_to_new=zeros(size(Pb,2),1);
    jj=0;
    for k=1:size(Pb,2)
        if Pb(1,k)>0.5 &&  Pb(2,k)>0.5
            v(k)=0; %point to remove
        else
          jj=jj+1; %if it is not to remove, keep track of the new index
          vec_new_to_old(jj)=k;
          vec_old_to_new(k)=jj;
      end
  end
  vec_new_to_old=vec_new_to_old(1:jj);
  toremove=find(v==0);
  T=geo.T;
  w=ones(size(T,2),1); %update triangle structure
  for k=1:size(T,2)
      if sum(find(T(1,k)==toremove))+sum(find(T(2,k)==toremove))+sum(find(T(3,k)==toremove))>0
          w(k)=0;
      else
          T(1,k)=vec_old_to_new(T(1,k));
          T(2,k)=vec_old_to_new(T(2,k));
          T(3,k)=vec_old_to_new(T(3,k));
      end
  end
  T=T(:,find(w==1));%remove triangle that are outside lsquare
  P=Pb(:,find(v==1));
  Pb=P;
  Tb=T(1:3,:);
  geo.P=P;geo.Pb=P;geo.T=T;geo.Tb=Tb; %update structures
  v=zeros(size(P,2),1); %get interior nodes in the new lsquare mesh
  tol=10^(-10);
  for k=1:length(v)
    if(abs(P(2,k))<tol || abs(P(1,k))<tol || (abs(P(1,k)-1)<tol) || (abs(P(2,k)-1)<tol) ||  (abs(P(2,k))>0.5-tol && abs(P(1,k))>0.5-tol))
        v(k)=-1;
    end
  end
  geo.int=find(v==0);
  geo.out=find(v==-1); 
  
  
end