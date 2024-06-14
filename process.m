%% Load scene
% option to visualize results
openfigures = true;

% scene to load
example = 'soap';

% load scene frame (for visualization and bilateral filtering)
load(['dataset/', example, '/scene.mat']);
scene = double(scene) / (2 ^ 14 - 1);
if (openfigures)
	figure; imshow(scene); title('Scene frame'); 
end

%% OCT "groundtruth" for comparison
load(['dataset/', example, '/depth.mat']);
depthOCT = depth;
depthOCT = depthOCT - max(depthOCT(:));
valsOCT = prctile(depthOCT(:), [1 99]);
if (openfigures)
	figure; imagesc(depthOCT); axis equal tight; colorbar; clim(valsOCT); title('OCT depth'); 
end

%% Swept-angle SWI preprocessing
% filter hyperparameters
spatialWindow = 21;
intensityWindow = 0.05;

% Load SWI frames
load(['dataset/', example, '/frames.mat']);
frames = double(frames) / (2 ^ 14 - 1);

% pixels that are zero across the entire phase shifting stack are invalid
validMask = ~(var(reshape(double(frames),...
			[size(frames, 1), size(frames, 2), 16]), [], 3) == 0);
validMask = validMask(:);

%% Swept-angle SWI with Poisson model and bilateral filtering
% compute phase
phase = reconstructSWI(frames, 'gaussian', 'bilateral', 'before',...
			spatialWindow, intensityWindow, scene);

% align phase with OCT depth
vals = prctile(phase(validMask), [1 99]);
depthSWI = valsOCT(2) + (phase - vals(1)) / (vals(2) - vals(1)) *...
					(valsOCT(1) - valsOCT(2));

if (openfigures)
	figure; imagesc(depthSWI); axis equal tight; colorbar; clim(valsOCT);
	title('Swept-angle SWI depth with bilateral filtering');
end
