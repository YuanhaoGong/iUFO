function [u,energy] = iUFO(im,iteration)
%Image Filtering with Generic Geometric Prior

% @ARTICLE{8470940, 
% author={Y. Gong and X. Hou and F. Li and G. Qiu}, 
% journal={IEEE Access}, 
% title={Image Filtering with Generic Geometric Prior}, 
% year={2018},  
% doi={10.1109/ACCESS.2018.2871829}, 
% ISSN={2169-3536},}

u = single(im);
f = u;
energy = zeros(iteration,size(u,3));
offset=uint32(reshape(1:size(u,1)*size(u,2),size(u,1),size(u,2)));
for ch=1:size(u,3)
    %main loop
    for i=1:iteration
        tmp = filter(u(:,:,ch),offset);
        u(:,:,ch) = u(:,:,ch) + tmp;
        e = u(:,:,ch) - f(:,:,ch);
        energy(i,ch)=sum(sum(e.^2));
    end
    %boundary
    u(1,:,ch)=u(2,:,ch);u(end,:,ch)=u(end-1,:,ch);u(:,1,ch)=u(:,2,ch);u(end,:,ch)=u(end-1,:,ch);
end

function dm=filter(u,offset)
dist = zeros([size(u),16],'single');
%line kernels
k1=[0,0,0;
    1,-2,1;
    0,0,0]/3;
k2=[1,0,0;
    0,-2,0;
    0,0,1]/3;
%edge kernels
k3=[1,1,0;
    1,-5,0;
    1,1,0]/6;
%mc kernels
k4=[1,2,0;
	4,-10,0;
	1,2,0]/10;
k5=[1,4,1/2;
    4,-10,0;
    1/2,0,0]/10;
%line
dist(:,:,1) = conv2(u,k1,'same');
dist(:,:,2) = conv2(u,k1','same');
dist(:,:,3) = conv2(u,k2,'same');
dist(:,:,4) = conv2(u,flipud(k2),'same');
%edge
dist(:,:,5) = conv2(u,k3,'same');
dist(:,:,6) = conv2(u,fliplr(k3),'same');
dist(:,:,7) = conv2(u,k3','same');
dist(:,:,8) = conv2(u,flipud(k3'),'same');
%wmc
dist(:,:,9) = conv2(u,k4,'same');
dist(:,:,10) = conv2(u,fliplr(k4),'same');
dist(:,:,11) = conv2(u,k4','same');
dist(:,:,12) = conv2(u,flipud(k4'),'same');
dist(:,:,13) = conv2(u,k5,'same');
dist(:,:,14) = conv2(u,fliplr(k5),'same');
dist(:,:,15) = conv2(u,flipud(k5),'same');
dist(:,:,16) = conv2(u,rot90(k5,2),'same');

%minimum absolute distance
tmp = abs(dist); 
[~,ind] = min(tmp,[],3);

%turn sub to index, but faster than sub2ind
dim2 = uint32(size(dist,1)*size(dist,2));
index = offset+uint32(ind-1)*dim2;
dm = dist(index); 