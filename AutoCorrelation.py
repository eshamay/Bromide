Ñò
ü¤Mc           @   s:
  d  d k  Z d  d k Z d   Z d   Z d   Z d   Z d   Z d   Z	 d   Z
 d	   Z d
   Z d   Z e d j o·	d Z e e  Z g  Z e D] Z e e d q [ Z e i   e i e d e d d d d e i d d  e i d d d g  Z xC e D]; Z e d Z e d j  o e i d  n e i d  qWe i   e i e d e d d d d e i d d  e i d d d e e  Z e e  Z e i   e i e e d d d d e i e e d d d d  e i d! d"  e i d e d   e i d# d d e e i  e  Z! e i" e!  Z# e i$ e! e# f  Z% e i& i& e%  Z' e i( e'  d  Z) e i& i* e)  e+ e  Z, e i   e i e e, d e+ e  ! e i d e d   e i d$ d d e i$ e i- e!  e# f  Z. e i& i& e.  Z/ e i& i* e i( e/  d   Z0 e i   e i e e0 d e+ e  ! e i d e d   e i d% d d e i   e i e e0 d e+ e  !e0 d d d& d d  e i e e, d e+ e  !e, d d d e i d e d   e i d' d d e e  Z1 e e  Z2 e i   e i e e1 d d( d d e i e e2 d d& d d  e i d! d"  e i d e d   e i d) d d g  Z3 xC e D]; Z e d Z e d j  o e3 i d  n e3 i d  q¹WxC e D]; Z e d Z e d j  o e3 i d  n e3 i d  qÿWe d d  e  Z4 g  Z5 xC e4 D]; Z e d Z e d j  o e5 i d  n e5 i d  q^We i   e i6 d  d d  e i e4 d e3 d d d d e i d d"  e i6 d  d d   e i e4 d e5 d d d d e i d d"  e i d* d d e e  Z7 e i   e i e e i  e  g e d+  e i e7 d d& e i d d  e i d e d   e i d, d d d- Z8 d- Z9 e i: e8 e9 f  Z; e i< i= d. d  Z> d/ Z? d0 Z@ e i< iA e? d  e8 e? d  e>  ZB e i< iA e@ d  e9 e@ d  e>  ZC x eD eB eC  D]~ \ ZE ZF eG eH d eE e? d   eI e8 eE e? d    ZJ eG eH d eF e@ d   eI e9 eF e@ d    ZK d e; eJ eK f <q]We i   e iL d1 e; iM d2 e iN iO e i d3 d d e
 e;  ZP e i   e iL eP d d4  d d4  f iM d2 e iN iQ e i d5 d d e e;  ZR e i   e iL eR d d4  d d4  f iM d2 e iN iQ e i d6 d d eP d d4  d f ZS eR d d4  d f ZT e i   e i e d d4  eS d d( d d e i e d d4  eT d d& d d  e i d7 d d e e;  ZU e i   e iL eU d d4  d d4  f iM d2 e iN iQ e i d8 d d e i   e i e d d4  eU d d4  d f  e i e d d4  e i  e;  g d4 d+  e i d9 d d n d S(:   iÿÿÿÿNc         C   sa   |  t  i |   } t  i i |  } t  i |  d } t  i i |  t |  } t  i |  S(   sJ  
    Given a 1D signal, return an estimation of its autocovariance
    function.

    The estimation is made by considering that the input signal
    actually desribes a full period of a wider, cyclic signal. The
    estimation is then the autocovariance of this wider signal.

    Uses the Fast Fourier Transform internally.
    i   (   t   npt   meant   fftt   abst   ifftt   lent   real(   t   signalt   centered_signalt	   ft_signalt   powerSpectralDensityt   autocovariance(    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftCyclicAutocovariance1DÆ  s
    c         C   sC   t  |   } | i d } | d j o t i | i  S| | Sd S(   sÑ   
    Given a 1D signal, return an estimation of its autocorrelation
    function.

    The autocorrelation is obtained by normalizing the autocovariance
    function computed by fftCyclicAutocovariance1D.
    i    g        N(   R   t   flatR    t   zerost   shape(   R   R   t   variance(    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftCyclicAutocorrelation1Dó  s
    c         C   sç   |  t  i |   } t  i |  } t  i | | f  } t  i i |  } t  i |  d } t  i i |  } t  i |  } t  i | | f  } t  i i |  }	 t  i i t  i |	  d  }
 | |
 } t  i | d t	 |   ! S(   sO  
    Compute the autocovariance of the input 1D signal.

    Consider the input signal to be a representative sample of a wider
    signal that has no other pattern that those present on the sample
    (this is what "representative" stands for) and especially no
    pattern whose scale is higer or equal to the input signal's size
    (this is for the difference with fftCyclicAutocovariance1D).

    The autocovariance is computed by a FFT and with a zero padding
    made in order to double the size of the signal. However the
    returned function is of the same size as the signal.
    i   i    (
   R    R   t
   zeros_liket   concatenateR   R   R   t	   ones_likeR   R   (   R   R   t   zero_paddingt   padded_signalR	   t   pseudo_powerSpectralDensityt   pseudo_autocovariancet   input_domaint   maskt   ft_maskt   mask_correction_factorsR   (    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftAutocovariance1Dt  s    	
	c         C   s@   t  |   } | d } | d j o t i | i  S| | Sd S(   sË   
    Given a 1D signal, return an estimation of its autocorrelation
    function.

    The autocorrelation is obtained by normalizing the autocovariance
    function computed by fftAutocovariance1D.
    i    g        N(   R   R    R   R   (   R   R   R   (    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftAutocorrelation1DÃ  s
    
c         C   sÔ   t  i |   } t  i |  | f  } t  i i |  } t  i |  d } t  i i |  } t  i t  } t  i | | f  } t  i i |  } t  i i t  i |  d  }	 | |	 }
 t  i |
 d t	 |   ! S(   s  
    Compute an autocovariance function without any centering, on a
    given binary signal considered as a discrete description of a set.

    The computed autocovariance will represent, for each translation
    vector h, the probability that a point belongs both to the initial
    set and its translated image.
    
    The measure computed by FFT, using a zero padding as with
    fftAutocovariance1D.
    i   i    (
   R    R   R   R   R   R   R   R   R   R   (   t   binarySignalR   R   R	   R   R   R   R   R   R   R   (    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftSetAutocovariance1D  s    	
	c         C   sa   |  t  i |   } t  i i |  } t  i |  d } t  i i |  t |  } t  i |  S(   sU  
    Given a n-dimensional signal, return an estimation of its
    autocovariance function.

    The estimation is made by considering that the input signal
    actually desribes a full period of a wider, cyclic signal. The
    estimation is then the autocovariance of this wider signal.

    Uses the Fast Fourier Transform internally.
    i   (   R    R   R   t   fftnR   t   ifftnR   R   (   R   R   R	   R
   R   (    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftCyclicAutocovariance]  s
    c         C   sC   t  |   } | i d } | d j o t i | i  S| | Sd S(   sÚ   
    Given a n-dimensional signal, return an estimation of its
    autocorrelation function.

    The autocorrelation is obtained by normalizing the autocovariance
    function computed by fftCyclicAutocovariance.
    i    g        N(   R#   R   R    R   R   (   R   R   R   (    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftCyclicAutocorrelation  s
    c         C   s  |  t  i |   } g  } | i D] } | d | d q! ~ } t  i i | |  } t  i |  d } t  i i |  } t  i |  } t  i i | |  }	 t  i i t  i |	  d  }
 | |
 } g  } |  i D] } | t |  qÔ ~ } t  i	 | |  S(   sÇ  
    Compute the autocovariance of the input n-dimensional signal.

    Consider the input signal to be a representative sample of a wider
    signal that has no other pattern that those present on the sample
    (this is what "representative" stands for) and especially no
    pattern whose scale is higer or equal to the input signal's size
    on each of its dimensions (this is for the difference with
    fftCyclicAutocovariance).

    The autocovariance is computed by a FFT and with a zero padding
    made in such a way that the padded signal is `2**n` bigger than
    the input one (where n is the dimension). However the returned
    function is of the same size as the signal on every dimension.
    i   i   (
   R    R   R   R   R!   R   R"   R   t   sliceR   (   R   R   t   _[1]t   st   padded_shapeR	   R   R   R   R   R   R   t   _[2]t   it   crop_slices(    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftAutocovariance°  s    	,
	*c         C   sC   t  |   } | i d } | d j o t i | i  S| | Sd S(   sÔ   
    Given a n-dimensional signal, return an estimation of its
    autocorrelation function.

    The autocorrelation is obtained by normalizing the autocovariance
    function computed by fftAutocovariance.
    i    g        N(   R,   R   R    R   R   (   R   R   R   (    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftAutocorrelationû  s
    c         C   sî   g  } |  i  D] } | d | d q ~ } t i i |  |  } t i |  d } t i i |  } t i |   } t i i | |  } t i i t i |  d  }	 | |	 }
 g  } |  i  D] } | t |  qÁ ~ } t i |
 |  S(   s  
    Compute an autocovariance function without any centering, on a
    given binary signal considered as a discrete description of a set.

    The computed autocovariance will represent, for each translation
    vector h, the probability that a point belongs both to the initial
    set and its translated image.
    
    The measure computed by FFT, using a zero padding as with
    fftAutocovariance.
    i   i   (	   R   R    R   R!   R   R"   R   R%   R   (   t   binary_signalR&   R'   R(   R	   R   R   R   R   R   R   R)   R*   R+   (    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   fftSetAutocovariance"  s    ,
	*t   __main__iè  i(   i    t   colort   kt	   linewidthg      ð?i-   s
   signal.pngt   dpiiÈ   ic   i   i   g      ø?s   signal_bin.pngt   redi   t   bluei   gÍÌÌÌÌÌð¿gÍÌÌÌÌÌð?s   autocorr_comparison.pngs   erroneous_autocov_signal.pngs   autocov_mask.pngt   bs   autocov_mask_hull.pngt   rs   bias_of_cyclic_estimation.pngs   signal_bin_cyclic.pngs   k:s   set_autocov.pngi   i   id   i<   g     ào@t   cmaps   boolean_model_of_rectangles.pngiÿ   s   cyclic_autocorr_img.pngs   truncated_autocorr_img.pngs   autocorr_img_profile_x.pngs   set_autocov_img.pngs   set_autocov_img_profile_x.png(V   t   pylabt   plt   numpyR    R   R   R   R   R    R#   R$   R,   R-   R/   t   __name__t   domain_lengtht   ranget   domainR&   t   stepR   t   figuret   vlinest   ylimt   savefigt
   signal_bint
   step_mod99t   appendt   cyclic_autocorrt   truncated_autocorrt   plott   xlimR   R   R   R   R   R   R   R	   R   R   R   R   R   R   R   R   R   t   cyclic_autocorr_signal_bint   truncated_autocorr_signal_bint   signal_bin_cyclict
   domain_padt   signal_bin_stationaryt   subplott   set_autocovariancet	   img_widtht
   img_heightR   t   imgt   randomt   poissont   num_rectanglest
   rect_widtht   rect_heightt   randintt   x_coordst   y_coordst   zipt   xt   yR%   t   maxt   mint   rect_drawing_area_xt   rect_drawing_area_yt   imshowt   Tt   cmt   grayt   cyclic_autocorr_imgt   jett   truncated_autocorr_imgt   cyclic_autocorr_profilet   truncated_autocorr_profilet   set_autocovariance_img(    (    (    s&   /home/shachar/work/Bromide/Autocorr.pyt   <module>T   s2  ÿ r	-		O	L	N	-	&	K	'	G%
 


$	


$$
 
 

 
	

#(		 
 
2
2
%%'
2
),