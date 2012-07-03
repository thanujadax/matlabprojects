Gabor Filtering Toolbox v 1.0.0 for Matlab
------------------------------------------


Changelog
---------
Version 1.0.0
  - First stable release. Includes now:
    * sg_version.m
    * sg_filterwithbank2.m
    * Contents.m
    * m2html generated documentation 


Version 0.2.2
  - Improved 'verbose' option for sg_createfilterbank,
    verbose=1 is now faster and verbose=2 displays all
    created filters in both frequency and spatial domain.


Version 0.2.1
  - Dependency to statistics toolbox was dropped by own
    implementation of norminv() in sg_createfilterf2.m


Version 0.2
  - The package now includes getargs.m.
  - sg_createfilterf2.m assigns all output arguments with
    all function parameters.


Version 0.1
  - Initial release



Short introduction to usage of Gabor toolbox
--------------------------------------------

First, an image is loaded and converted to double. The image has
likely three color channels even if it is a gray-scale image, so only
one of them must be selected. While not strictly necessary, the image
can be also scaled to [0,1] instead of normal [0,255].

  image=imread('someimage.jpg'); % load image
  image=image(:,:,1);            % select only one color channel
  image=double(image)./256;      % convert to double and scale to [0,1]

Then, a filterbank can be created. 

 bank=sg_createfilterbank(size(image), 0.2 , 5, 4,'verbose',1);

creates a filterbank with the size of the image, the frequency of the
highest frequency filter is 0.2, 5 filters at different frequencies
and 4 orientations are created. The created filter bank will be
displayed. Only half of the filterbank is created, because 
responses for the second half of the filter bank are complex
conjugates of the responses from the first half.

The filterbank can be used to filter images.

  r=sg_filterwithbank(image,bank,'method',1);

The responses will be returned in a special structure. The structure
can be converted to a 3-d matrix by using

  m=sg_resp2samplematrix(r);

Now, m will be a 512x512x20 (or whatever the image resolution was x20)
matrix, since there are 5*4=20 Gabor filters.  If you do not have a
good idea what to do with the responses, you can for example view 
them as an image by summing all the responses:

  imagesc(abs(sum(m,3))); colormap(gray);



For object detection and localization functionality other functions
are also available.

Sample matrix can be normalized for illumination invariance:

  m_norm=sg_normalizesamplematrix(m);


For scale invariance, extra frequencies can be first included in 
the filter bank:
  
  bank=sg_createfilterbank(size(image), 0.2 , 5, 4,'extra_freq',1);
  r=sg_filterwithbank(image,bank);
  m=sg_resp2samplematrix(r);

Then, different scales of features can be selected from sample
matrix:

  m2=sg_scalesamples(m,0,5,4);

will select the 5 highest frequencies and

  m2=sg_scalesamples(m,1,5,4);
  
the lower available 5 frequencies. Note that if normalization is
required, it must be done after this sg_scalesamples step.


For rotation invariance, the sample matrix can be rotated:

  m2=sg_rotatesamples(m,1,4);

rotates the responses by 1 orientation.
