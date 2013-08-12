function rot = filterImageWithMembraneTemplateRotated(im, d)

dim = size(im);

  a = pi / 8;
  t = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
  dt = centeredRotate(d, 0);
%   rot1 = single(normxcorr2_mex(double(dt), double(im),'same'));
  rot1 = normxcorr2(double(dt), double(im));
  rot1 = removeBorder(rot1,dim);
  rot1 = single(rot1);
  rot = single(zeros(size(rot1,1), size(rot1,2), 8));
  rot(:,:,1) = single(rot1);
  clear rot1;
  
  dt = centeredRotate(d, a);
%   rot2 = normxcorr2_mex(double(dt), double(im),'same');
  rot2 = normxcorr2(double(dt), double(im));
  rot2 = removeBorder(rot2,dim);
  rot(:,:,2) = single(rot2);
  clear rot2;
  
  dt = centeredRotate(d, 2*a);
  rot3 = normxcorr2(double(dt), double(im));
  rot3 = removeBorder(rot3,dim);
  rot(:,:,3) = single(rot3);
  clear rot3;

  dt = centeredRotate(d, 3*a);
  rot4 = normxcorr2(double(dt), double(im));
  rot4 = removeBorder(rot4,dim);
  rot(:,:,4) = single(rot4);
  clear rot4;

  dt = centeredRotate(d, 4*a);
  rot5 = normxcorr2(double(dt), double(im));
  rot5 = removeBorder(rot5,dim);
  rot(:,:,5) = single(rot5);
  clear rot5;

  dt = centeredRotate(d, 5*a);
  rot6 = normxcorr2(double(dt), double(im));
  rot6 = removeBorder(rot6,dim);
  rot(:,:,6) = single(rot6);
  clear rot6;

  dt = centeredRotate(d, 6*a);
  rot7 = normxcorr2(double(dt), double(im));
  rot7 = removeBorder(rot7,dim);
  rot(:,:,7) = single(rot7);
  clear rot7;

  dt = centeredRotate(d, 7*a);
  rot8 = normxcorr2(double(dt), double(im));
  rot8 = removeBorder(rot8,dim);
  rot(:,:,8) = single(rot8);
  clear rot8;

