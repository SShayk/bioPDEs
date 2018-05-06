% based on neuron drawn in .pdf

% note - should make this generalizable (can load # junctions & branch
% lengths from files)

clear

%% set up main junction structure
N = 1:9;                    % indices of junctions
n_junc = length(N);         % number of nodes
R = [0,1,1,1,3,3,2,2,8];    % rootward neighbor of each junction

C = {1,n_junc};             % children of each junction
FC = zeros(1,n_junc);       % first child of each junction
NS = zeros(1,n_junc);       % next sibling of each junction

level = zeros(1,n_junc);         % # of junctions between each junction and the root junction

for ii = N
    
    % find all children of current junction
    C{ii} = find(R==N(ii));
    % set first child of current junction
    if C{ii} 
        FC(ii) = C{ii}(1);
    end
        
end

% find distance between each junction and the root junction
for ii = N
temp_R = R(ii);
level_ct = 0;
    while(temp_R)
        temp_R = R(temp_R);
        level_ct = level_ct +1;
    end
    level(ii) = level_ct;
end
% set next sibling of each junction
for ii = 1:max(level)
    % get junctions at current level
    juncs_level = find(level==ii);
    % set next siblings of all junctions at current level
    NS(juncs_level(1:end-1)) = juncs_level(2:end);
end

%% set up node structure

L = zeros(1,n_junc-1);             % length of each branch [um]
N12 = zeros(2,n_junc-1);           % starting and ending junction of each branch

% set up lengths of each branch [um]
L(1) = 50;
L(2) = 10.34;
L(3) = 3.2;
L(4) = 34.5;
L(5) = 23.5;
L(6) = 86;
L(7) = 12.3;
L(8) = 22.4;

dx_ref = .5;                       % desired dx [um]
dx = zeros(1, n_junc-1);           % actual dx of each branch [um]
nseg = zeros(1, n_junc-1);         % # of nodes in each branch [um]

for ii = 1:length(L)
    % get # of nodes in branch
    nseg(ii) = ceil(L(ii)/dx_ref); 
    % get dx of branch
    dx(ii) = L(ii)/nseg(ii);
end

% set up start/end nodes of each branch
for ii = 2:n_junc
    N12(:,ii-1) = [R(ii); N(ii)];
end

% number all nodes and mark with respective branches
n = 0;
n_index = zeros(1,sum(nseg(:)));
branch_start = zeros(size(L));


for ii = 1:length(L)    % over all branches
    for jj = 1:nseg(ii) % over all segments in branch
        if jj == 1
            if ii == 1
                n = 1;
                n_index(1) = 1;
                branch_start(1) = 1;
            elseif find(NS==ii) 
                branch_start(ii) = branch_start(find(NS==ii));
            elseif
                
            else
                branch_start(ii) = n+1;
            end
        else
             n = n +1;
        end
           
    end
end

% mark node index of each junction

