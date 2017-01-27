function faces=findfracboundaryfaces2d(G,tol)
% Generates a list of global faces belonging to fractures. These global
% faces intersect with the external boundary with the domain. This is to
% facilitate assigning boundary conditions to these faces. The list of
% faces can be plugged into addBC. The inputs are the global grid which
% needs to contain fracture grid data. Fractures crossing the domain
% boundary are not considered.
%
% The output is a struct class object containing the fields West, East,
% South, North with each field containing a list of
% fracture faces.
%
% Note that this function only works for 2D rectangular domains.

faces=struct('West',[],'East',[],'South',[],'North', []);

%% Extract external faces
mincoords=min(G.Matrix.nodes.coords);
maxcoords=max(G.Matrix.nodes.coords);

% West face
westpoint=mincoords;
westnormal=[1 0];

% East face
eastpoint=maxcoords;
eastnormal=[-1 0];

% South face
southpoint=mincoords;
southnormal=[0 1];

% North face
northpoint=maxcoords;
northnormal=[0 -1];

%% Extract boundary faces
numfracs=length(fieldnames(G.FracGrid));
for i=1:numfracs
    fieldname=['Frac',num2str(i)];
    F=G.FracGrid.(fieldname);
    facestart=F.faces.start;
    numnodes = F.nodes.num;
        
    % first treat faces at beginning and end of fracture, faces are part of 
    % a boundary if the face centroid lies on the respective boundary
    % get indices of these faces (start and end point of fracture)
    n = uniquetol(F.faces.normals(all(F.faces.neighbors,2),:),'ByRows',true); % get normal
    nfaces = all(abs(F.faces.normals-n) < 1e4*eps(min(abs(F.faces.normals),abs(n))),2); % get faces with that normal
    idx = find(any(nfaces & (~all(F.faces.neighbors,2)),2));
    for j = idx'
        % special treatment for fractures starting/ending in corner points
        % that are not parallel to either x- or y-axis
        if (any(F.faces.centroids(j,:) == [mincoords; maxcoords; mincoords(1), maxcoords(2); maxcoords(1), mincoords(2)],2))
            if ((find(n)==find(westnormal) && (pointplanedistance(F.faces.centroids(j,:),westnormal,westpoint)<tol))) % West
                faces.West(end+1)=facestart-1+j;
            elseif (find(n)==find(eastnormal) && pointplanedistance(F.faces.centroids(j,:),eastnormal,eastpoint)<tol) % East
                faces.East(end+1)=facestart-1+j;
            elseif (find(n)==find(southnormal) && pointplanedistance(F.faces.centroids(j,:),southnormal,southpoint)<tol) % South
                faces.South(end+1)=facestart-1+j;
            elseif (find(n)==find(northnormal) && pointplanedistance(F.faces.centroids(j,:),northnormal,northpoint)<tol) % North
                faces.North(end+1)=facestart-1+j;
            end
        else
            if (pointplanedistance(F.faces.centroids(j,:),westnormal,westpoint)<tol) % West
                faces.West(end+1)=facestart-1+j;
            elseif (pointplanedistance(F.faces.centroids(j,:),eastnormal,eastpoint)<tol) % East
                faces.East(end+1)=facestart-1+j;
            elseif (pointplanedistance(F.faces.centroids(j,:),southnormal,southpoint)<tol) % South
                faces.South(end+1)=facestart-1+j;
            elseif (pointplanedistance(F.faces.centroids(j,:),northnormal,northpoint)<tol) % North
                faces.North(end+1)=facestart-1+j;
            end 
        end
    end
    
    % as fracture bc are not considered for fractures intersecting the
    % boundary, the only other possible situation is a fracture aligning 
    % with one of the boundaries
    if (size(uniquetol(F.cells.centroids(:,1), 'ByRows',true),1)==1) % // to y-axis
        fline = [accumarray(kron(1:numnodes/2,[1 1])', F.nodes.coords(:,1), [numnodes/2 1])/2, accumarray(kron(1:numnodes/2,[1 1])', F.nodes.coords(:,2), [numnodes/2 1])/2];
        if pointplanedistance(uniquetol(fline(:,1)),eastnormal,eastpoint)<tol % East
        faces.East = [faces.East, facestart-1+find(F.faces.centroids(:,1) == max(F.faces.centroids(:,1)))];
        elseif pointplanedistance(uniquetol(fline(:,1)),westnormal,westpoint)<tol % West  
        faces.West= [faces.West, facestart-1+find(F.faces.centroids(:,1) == min(F.faces.centroids(:,1)))];
        end
    elseif (size(uniquetol(F.cells.centroids(:,2), 'ByRows',true), 1)==1) % // to x-axis
        fline = [accumarray([1:numnodes/2, 1:numnodes/2]', F.nodes.coords(:,1), [numnodes/2 1])/2, accumarray([1:numnodes/2, 1:numnodes/2]', F.nodes.coords(:,2), [numnodes/2 1])/2];
        % if north boundary, assign upper fracture faces
        if pointplanedistance(uniquetol(fline(:,2)),northnormal,northpoint)<tol % North
        faces.North= [faces.North, facestart-1+find(F.faces.centroids(:,2) == min(F.faces.centroids(:,2)))];
        % if south boundary, assign lower fracture faces  
        elseif pointplanedistance(uniquetol(fline(:,2)),southnormal,southpoint)<tol % South
        faces.South= [faces.South, facestart-1+find(F.faces.centroids(:,2) == min(F.faces.centroids(:,2)))];
        end
    end
end 
end