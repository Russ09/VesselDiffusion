  %Example,
  % compile needed mex file
  %mex eig3volume.c

  load('ExampleVolumeStent');
  
  % Frangi Filter the stent volume
  options.BlackWhite=false;
  options.FrangiScaleRange=[1 1];
  Vfiltered=FrangiFilter3D(V,options);

  % Show maximum intensity plots of input and result
  figure, 
  subplot(2,2,1), imshow(squeeze(max(V,[],2)),[])
  subplot(2,2,2), imshow(squeeze(max(Vfiltered,[],2)),[])
  subplot(2,2,3), imshow(V(:,:,100),[])
  subplot(2,2,4), imshow(Vfiltered(:,:,100),[])
