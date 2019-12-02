function outMap = createOutmap(bitMap,framesPerUp)
% image gets transposed and flipped vertically due to reflections and
% the fact that openGL is in a right handed coordinate system. so
% always rotrate 90 deg clockwise to counteract -MSC

% NOTE: This is a modified version of psycho5's CreateTexture.m

%make sure no bytes wrap around
bitMap(bitMap>1) = 1;
bitMap(bitMap<0) = 0;
outMap(size(bitMap,1),size(bitMap,2),4) = 0;

%%%%%% THIS ONLY WORKS IF framesPerUp IS DIVISIBLE BY 3
% bitDepth 7 (essentially 8 bit) 4bit 2bit and 1bit
% number of frames per color rgb. Its also not necessary if framesPerUp
% = 3

fpc = framesPerUp/3;
% puts the bitMap into values between 0 and the bitsize

nBitMap = round(bitMap*(2^(8/fpc)-1));
%wBitMap = zeros([size(nBitMap) 1 1 1]);
weights = repmat(2.^(0:8/fpc:7),[1 3]);
weights = permute(weights,[1 3 2]);

% Make weights the same variable type as nBitMap
weights = uint8(weights);

% weight the matrix so that the bits are properly shifted
wBitMap = bsxfun(@times,nBitMap,weights);

% now that the bits are aligned simply sum them together
for ii = 1:3
    outMap(:,:,ii) = sum(wBitMap(:,:,(ii-1)*fpc+1:ii*fpc),3);
end

% add alpha values to outMap. They won't be used but they keep
% makeTexture from converting the format, slowing down, and droping
% frames.
repX = ceil(size(outMap,2)/2);
repY = size(outMap,1);
alpha = repmat([100 200],[repY,repX]);
outMap(:,:,4) = alpha(:,1:size(outMap,2));

outMap = outMap(:,:,[2 1 3 4]);

% tex = Screen('MakeTexture', windowPtr, outMap, [], 1);

end