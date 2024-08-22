require(dplyr)
require(tidyr)
require(ggplot2)
require(PaleoSpec)

############add CI intervals to DML series
source("./FilterSpec.R")

ice<-read.csv("ice_cores_dml_snr[1].csv")

#calculate alpha and beta values and plot signal curve 
ice_sig_betas <- lm(log(signal)~log(freq), data = ice)
signalBeta = summary(ice_sig_betas)$coefficients[2]
signalAlpha = summary(ice_sig_betas)$coefficients[1]
signalAlphaExp<-exp(signalAlpha)

plot(ice$freq, ice$signal, type = "l", log = "xy")
lines(ice$freq, signalAlphaExp * ice$freq^(signalBeta), col = "red")

#calculate noise alpha and beta values
ice_noise_betas <- lm(log(noise)~log(freq), data = ice)
noiseBeta = summary(ice_noise_betas)$coefficients[2]
noiseAlpha = summary(ice_noise_betas)$coefficients[1]

#empty dataframe
sims<-data.frame(matrix(ncol = 5, nrow  = 0))
colnames(sims)<-c("freq","noise","signal","ratio","run")

#take real alpha and beta signal and noise estimates and create synthetic data
for(x in 1:100){
  run<-as.character(x)
  length<-995 #995 for DML, 475 for NGT -> length of synthetic time series
  sig_beta<-signalBeta*-1 
  noise_beta<-noiseBeta*-1
  sig_alpha<-exp(signalAlpha)
  noise_alpha<-exp(noiseAlpha)
  n_sites<-3 #3 for DML, 14 for NGT -> simluate number of records in cluster
  simclim<-ts(SimPLS(N = length, beta = sig_beta, alpha = sig_alpha)) #simulate signal spectrum
  noise_list<-list()
  noise_df<-data.frame(matrix(ncol = length, nrow = 0))
    for (j in seq_along(1:n_sites)){
      #for n sites in each cluster, simulate a series with real noise parameters and add to the signal spec 
      sim<-ts(SimPLS(N = length, beta = noise_beta, alpha = noise_alpha))
      simplus<-sim+simclim
      simspec<-SpecMTM(simplus)
      simspec<-list(simspec)
      #combine all site spectra into a stack
      noise_list<-c(noise_list, simspec)
      noise_df[nrow(noise_df) + 1,] = simplus
    }
    #calculate the SNRs with synthetic data
    colnames(noise_df)<-seq_along(1:length)
    noise_df<-ts(colMeans(noise_df))
    spectrum_of_stack<-SpecMTM(noise_df)
    freq<-spectrum_of_stack$freq
    noise_list<-lapply(noise_list, function(x) rbind(x$spec))
    noise_list<-do.call(rbind, noise_list)
    mean_of_spectra<-colMeans(noise_list)
    spec_ests<-cbind.data.frame(freq, spectrum_of_stack$spec, mean_of_spectra)
    colnames(spec_ests)<-c("freq","spectrum_of_stack","mean_of_spectra")
    #############################
    spec_ests$noise<-(spec_ests$mean_of_spectra - spec_ests$spectrum_of_stack)/(1-1/n_sites)
    spec_ests$signal<-spec_ests$mean_of_spectra - spec_ests$noise
    ###############################optional smoothing of signal and noise
    # spec_ests<-spec_ests%>%
    #   filter(!is.na(signal))
    # signal<-list(spec_ests$freq, spec_ests$signal)
    # names(signal)<-c("freq","spec")
    # class(signal)<-'spectrum'
    # signal<-FilterSpec(signal, spans = c(3,5))
    # signal<-FilterSpecLog(signal, df = 0.01)
    # spec_ests$signal<-signal$spec
    # noise<-list(spec_ests$freq, spec_ests$noise)
    # names(noise)<-c("freq","spec")
    # class(noise)<-'spectrum'
    #noise<-FilterSpec(noise, spans = c(3,5))
    #noise<-FilterSpecLog(noise, df = 0.01)
    #spec_ests$noise<-noise$spec
    #################################calculate ratio
    spec_ests$ratio<-spec_ests$signal/spec_ests$noise
    spec_ests$spectrum_of_stack<-NULL
    spec_ests$mean_of_spectra<-NULL
    spec_ests$run<-run
    sims<-list(sims, spec_ests)
    sims<-do.call(rbind, sims)
}


sims_mean_df<-sims %>%
  group_by(freq)%>%
  summarise(signal_sim = mean(signal, na.rm = TRUE),
            noise_sim = mean(noise, na.rm = TRUE),
            snr_sim = mean(ratio, na.rm = TRUE),
            signal_0.9 = quantile(signal, 0.9, na.rm = TRUE),
            signal_0.1 = quantile(signal, 0.1, na.rm = TRUE),
            noise_0.9 = quantile(noise, 0.9, na.rm = TRUE),
            noise_0.1 = quantile(noise, 0.1, na.rm = TRUE),
            snr_0.9 = quantile(ratio, 0.9, na.rm = TRUE),
            snr_0.1 = quantile(ratio, 0.1, na.rm = TRUE))


#calculate ratio of CI estimates to simulated curves to apply multiplicatively to the 'real' estimates
sims_mean_df$snr_0.1<-sims_mean_df$snr_0.1/sims_mean_df$snr_sim
sims_mean_df$snr_0.9<-sims_mean_df$snr_0.9/sims_mean_df$snr_sim

sims_mean_df$signal_0.1<-sims_mean_df$signal_0.1/sims_mean_df$signal_sim
sims_mean_df$signal_0.9<-sims_mean_df$signal_0.9/sims_mean_df$signal_sim

sims_mean_df$noise_0.1<-sims_mean_df$noise_0.1/sims_mean_df$noise_sim
sims_mean_df$noise_0.9<-sims_mean_df$noise_0.9/sims_mean_df$noise_sim


#smooth to create nice confidence intervals with no noise
freq<-sims_mean_df$freq
sims_mean_df<-sims_mean_df[2:10]
sims_smooth<-lapply(sims_mean_df, function(x) loess(x ~ freq, span = 0.2))
sims_smooth<-lapply(sims_smooth, function(x) cbind.data.frame(x$fitted))
sims_smooth<-do.call(cbind, sims_smooth)
colnames(sims_smooth)<-colnames(sims_mean_df)
sims_smooth$freq<-freq


#combine 'real' estimates with simulated estimates and CI intervals. merge by freq
#round freq and convert to character vector (and then back to numeric) to merge.
ice$freq<-round(ice$freq, digits = 3)
sims_smooth$freq<-round(sims_smooth$freq, digits = 3)

ice$freq<-as.character(ice$freq)
sims_smooth$freq<-as.character(sims_smooth$freq)

ice_join<-merge(ice, sims_smooth, by = "freq")
ice_join$freq<-as.numeric(as.character(ice_join$freq))

##plot estimates
ggplot(ice_join, aes(x=freq, y=signal*signal_0.9))+
  geom_line()+
  geom_line(data = ice_join, aes(x=freq, y=signal*signal_0.1))+
  geom_line(data = ice_join, aes(x=freq, y = signal), color = "red",linetype = "longdash")+
  scale_x_continuous(trans = c("log10"))+
  scale_y_continuous(trans = c("log10"))

ggplot(ice_join, aes(x=freq, y=noise*noise_0.9))+
  geom_line()+
  geom_line(data = ice_join, aes(x=freq, y=noise*noise_0.1))+
  geom_line(data = ice_join, aes(x=freq, y = noise), color = "red",linetype = "longdash")+
  scale_x_continuous(trans = c("log10"))+
  scale_y_continuous(trans = c("log10"))

ggplot(ice_join, aes(x=freq, y=snr*snr_0.9))+
  geom_line()+
  geom_line(data = ice_join, aes(x=freq, y=snr*snr_0.1))+
  geom_line(data = ice_join, aes(x=freq, y = snr), color = "red",linetype = "longdash")+
  scale_x_continuous(trans = c("log10"))

#write.csv(ice_join, file = "dml_withCI.csv")

