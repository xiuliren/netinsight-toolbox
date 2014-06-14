---------------------------------------------------------------------------
This is the matlab implementation of following noise level estimation:

Xinhao Liu, Masayuki Tanaka and Masatoshi Okutomi
Noise Level Estimation Using Weak Textured Patches of a Single Noisy Image
IEEE International Conference on Image Processing (ICIP), 2012.

---------------------------------------------------------------------------

Copyright (C) 2012 Masayuki Tanaka. All rights reserved.
mtanaka@ctrl.titech.ac.jp 

---------------------------------------------------------------------------
 Contents
---------------------------------------------------------------------------
* NoiseLevel.m
 The main code of the noise level estimation.
 You can show the description by
 > help NoiseLevel
 
 demo.m also includes simple usage.
 This algorithm is implemented with only single m file. 
 

* demo.m
 Demonstration example.

* lena.png
 Sample image.

* README.txt
 This file.

---------------------------------------------------------------------------
 Note
---------------------------------------------------------------------------
We used the maximum eigenvalue of the gradient covariance matrix for the 
ICIP paper above. But, in this matlab code, we use the trace of the gradient 
covariance matrix instead of the maximum eigenvalue of that.

---------------------------------------------------------------------------
 Web page
---------------------------------------------------------------------------
http://bit.ly/NLest (http://www.ok.ctrl.titech.ac.jp/res/NLE/noise_level.html)


---------------------------------------------------------------------------
 License
---------------------------------------------------------------------------
This program is only for non-commercial use.
If you want any commercial use, please contact me.
The patent is applying.


---------------------------------------------------------------------------
 Change log
---------------------------------------------------------------------------
* 20120529 Released
