The CW-WTC method requires the support of the Cross Wavelet and Wavelet Coherence Toolbox.http://grinsted.github.io/wavelet-coherence/!

It should be noted that this method is only an optimisation for the WTC method, which focuses on setting the period range, extracting the maximum significant domain therein,
and calculating the period and phase angle based on a time-by-time weighted average of the wavelet coherence coefficients. 
The phase angle is then converted to time lag. The maximum significant domain is judged by the set cycle range and the range within the non-COI.
However, the wavelet coherence coefficients, period and phase angles are extracted later and the significant domain is not cut by the period range or COI.
