function [phase, amplitude, dc] = reconstructSWI(...
						frames, noiseModel,...
						filterType, filterOrder, spatialWindow, intensityWindow,...
						scene)
	% Reconstruct depth from synthetic wavelength interferometry
	% frames:			HxWx4x4 array of measurements, where the third dimension
	%					varies over subwavelength shifts and fourth over
	%					four-bucket positions (assumed in [0 1] range)
	% noiseModel:		'mc', or 'gaussian' (use amplitude estimation algorithm
	%					in the paper or based on Gaussian MLE, respectively
	% filterType:		'none', 'gaussian', or 'bilateral' (type of
	%					filtering to use)
	% filterOrder:		'before', or 'after' (filter envelopes or phase,
	%					respectively)
	% spatialWindow:	spatial window (for Gaussian and bilateral filtering)
	% intensityWindow:	intensity window (for bilateral filtering)
	% scene:			intensity image (for bilateral filtering, assumed in 
	%					[0 1] range) 

	if (nargin < 2)
		noiseModel = 'gaussian';
	end
	
	if (nargin < 3)
		filterType = 'none';
	end
	
	if (nargin < 4)
		filterOrder = 'before';
	end

	Im0Carrier = frames(:, :, 1, :);
	Im1Carrier = frames(:, :, 2, :);
	Im2Carrier = frames(:, :, 3, :);
	Im3Carrier = frames(:, :, 4, :);
	if (strcmp(noiseModel, 'mc'))
		[~, amplitudeCarrier] = mle_mc(Im0Carrier, Im1Carrier, Im2Carrier, Im3Carrier);
	elseif (strcmp(noiseModel, 'gaussian'))
		[~, amplitudeCarrier] = mle_gaussian(Im0Carrier, Im1Carrier, Im2Carrier, Im3Carrier);
	end
	
	if (max(imag(amplitudeCarrier(:))) < 10 ^ -7)
		amplitudeCarrier = real(amplitudeCarrier);
	end
	envelope = amplitudeCarrier .^ 2;
	envelope = squeeze(envelope);
  
	if (strcmp(filterOrder, 'before'))
		if (strcmp(filterType, 'gaussian'))
			% Filter the measured envelope with bilateral filtering using an ambient
			% light image of the scene as the guide image
			envelope = imgaussfilt(envelope, spatialWindow);
		elseif (strcmp(filterType, 'bilateral'))
			% Filter the measured envelope with bilateral filtering using an ambient
			% light image of the scene as the guide image
			for position = 1:4
				envelope(:, :, position) = bilateralFilter(envelope(:, :, position),...
															scene, 0, 1,...
															spatialWindow,...
															intensityWindow);
			end
		end
	end
  
	Im0Synthetic = envelope(:, :, 1);
	Im1Synthetic = envelope(:, :, 2);
	Im2Synthetic = envelope(:, :, 3);
	Im3Synthetic = envelope(:, :, 4);
	if (strcmp(noiseModel, 'mc'))
		[dc, amplitude, sina, cosa] = mle_mc(...
					Im0Synthetic, Im1Synthetic, Im2Synthetic, Im3Synthetic);
	elseif (strcmp(noiseModel, 'gaussian'))
		[dc, amplitude, sina, cosa] = mle_gaussian(...
					Im0Synthetic, Im1Synthetic, Im2Synthetic, Im3Synthetic);
	end
	
	if (strcmp(filterOrder, 'after'))
		if (strcmp(filterType, 'gaussian'))
			% Filter the recovered phase with bilateral filtering using an ambient
			% light image of the scene as the guide image
			sina = imgaussfilt(sina, spatialWindow);
			cosa = imgaussfilt(cosa, spatialWindow);
		elseif (strcmp(filterType, 'bilateral'))
			% Filter the recovered phase with bilateral filtering using an ambient
			% light image of the scene as the guide image
			sina = bilateralFilter(sina, scene, 0, 1, spatialWindow,...
									intensityWindow);
			cosa = bilateralFilter(cosa, scene, 0, 1, spatialWindow,...
									intensityWindow);
		end
	end
	phase = atan2(sina, cosa);  
end

function [dc, amplitude, sina, cosa] = mle_mc(Im0, Im1, Im2, Im3, exact)

	if (nargin < 5)
		exact = false;
	end

	dc = (Im0 + Im1 + Im2 + Im3) / 4;
	amplitude = sqrt((Im0 - dc) .^ 2 +...
					(Im1 - dc) .^ 2 +...
					(Im2 - dc) .^ 2 +...
					(Im3 - dc) .^ 2) / sqrt(8);

	if (nargout >= 3)
		if (exact)
			sina = (Im3 - Im1) ./ amplitude / 4;
			cosa = (Im0 - Im2) ./ amplitude / 4;
		else
			sina = (Im3 - Im1);
			cosa = (Im0 - Im2);
		end
	end
end

function [dc, amplitude, sina, cosa] = mle_gaussian(Im0, Im1, Im2, Im3, exact)

	if (nargin < 5)
		exact = false;
	end

	dc = (Im0 + Im1 + Im2 + Im3) / 4;
	factorEg = (Im0 - Im2) .^ 2 + (Im3 - Im1) .^ 2;
	amplitude = sqrt(factorEg) / 4;

	if (nargout >= 3)
		if (exact)
			sina = (Im3 - Im1) ./ amplitude / 4;
			cosa = (Im0 - Im2) ./ amplitude / 4;
		else
			sina = (Im3 - Im1);
			cosa = (Im0 - Im2);
		end
	end
end
