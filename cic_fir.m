function Hd = cic_fir(decf,diffd,numsecs)

Hd = dsp.CICDecimator('DecimationFactor',decf,...
                      'DifferentialDelay',diffd,...
                      'NumSections',numsecs);

end
