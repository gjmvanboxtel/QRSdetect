# afonso.R - QRS detection algorithm proposed by Afonso et al. (1999)
# Copyright (C) 2019  Geert van Boxtel
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Version history
# 20190117  GvB       Initial setup for package QRSdetect
#
#---------------------------------------------------------------------------------------------------------------------

#' Afonso et al. QRS detection
#'
#' Detect R peaks of the QRS complex in a raw ECG record, based on the filter-bank
#' algorithm proposed by Afonso et al. (1999)
#'
#' The Afonso et al. algorithm uses filter banks (polyphase implementation) and
#' determines candidate R-peaks on downsampled signals in different frequency
#' bands. Initially, a large number of false positives are generated. Then, logic
#' is added to decrease the number of false positives while maintaining a low
#' number of false negatives. This method is currently one of the most accurate
#' available (> 99.5% accuracy in various databases), and also one of the fastest
#' (because the logic is applied on downsampled signals).
#'
#' The present R implementation is based on the Matlab/Octave version named
#' nqrsdetect.m, Copyright (C) 2006 by Rupert Ortner, retrieved from the internet
#' pages maintained by Alois Schloegl (http://pub.ist.ac.at/~schloegl/). A few
#' improvements and minor bug fixes were made, as well as comments added.
#'
#' @param ecg The input single-channel input vector (raw ECG)
#' @param fs The frequency in Hz with which the ecg was sampled
#'
#' @return Numeric array containing the indices (sample numbers) at which the fiducial R-peaks were found
#'
#' @references Afonso, V.X., Tompkins, W.J., Nguyen, T.Q., & Luo, S.
#' (1999). ECG beat detection using filter banks. IEEE Transactions on Biomedical
#' Engineering, 46(2), 192-202. DOI: \url{https://dx.doi.org/10.1109/10.740882}{10.1109/10.740882}
#"
#' @author Geert van Boxtel
#'
#' @examples
#' data(ecg)
#' fs <- 360
#' pks <- afonso(ecg$MLII, fs)
#'
#' \dontrun{
#' # plot first 5 seconds of data
#' N <- 5 * fs
#' plot (ecg$time[1:N], ecg$MLII[1:N], type = "l", main = "MIT-BIH database, record 100",
#'       xlab = "Time (s)", ylab = "Amplitude (mV)")
#' points (pks[which(pks<=N)]/fs, ecg$MLII[pks[which(pks<=N)]], col="red")
#' }
#' @export

afonso <- function (ecg, fs) {

  # Calculate analysis filter order, bandwidth, and downsampling rate
  N <- round(fs)         #filter order
  Bw <- 5.6              #filter bandwidth
  Bwn <- 1/(fs/2)*Bw     #normalised filter bandwidth
  M <- round((fs/2)/Bw)  #downsampling rate

  # Calculate analysis filter bandwidths
  Wn0 <- Bwn #In nqrsdetect, but not used
  Wn1 <- c(Bwn, 2*Bwn)
  Wn2 <- c(2*Bwn, 3*Bwn)
  Wn3 <- c(3*Bwn, 4*Bwn)
  Wn4 <- c(4*Bwn, 5*Bwn)

  #impulse response of the analysis filters
  #h0 <- signal::fir1(N, Wn0) #In nqrsdetect, but not used
  h1 <- signal::fir1(N, Wn1, 'pass')
  h2 <- signal::fir1(N, Wn2, 'pass')
  h3 <- signal::fir1(N, Wn3, 'pass')
  h4 <- signal::fir1(N, Wn4, 'pass')

  # Downsample and filter with polyphase implementation
  downdim <- ceiling(length(ecg) / M)
  y <- matrix(0, downdim, 5)
  # y[,1] <- polyphase(ecg, h0, M) In nqrsdetect, but not used
  y[,2] <- polyphase(ecg, h1, M)
  y[,3] <- polyphase(ecg, h2, M)
  y[,4] <- polyphase(ecg, h3, M)
  y[,5] <- polyphase(ecg, h4, M)

  # Cut off initial transient because of filtering
  cut <- ceiling(N/M)
  # y1 <- c(rep(0, cut), y[cut:length(y[,1]),1]) In nqrsdetect, but not used
  y2 <- c(rep(0, cut), y[cut:length(y[,2]), 2])
  y3 <- c(rep(0, cut), y[cut:length(y[,3]), 3])
  y4 <- c(rep(0, cut), y[cut:length(y[,4]), 4])
  y5 <- c(rep(0, cut), y[cut:length(y[,5]), 5])

  # Compute composite subbands
  # Afonso et al. (1999), p. 195, eq. (13)
  P1 <- abs(y2) + abs(y3) + abs(y4)
  P2 <- abs(y2) + abs(y3) + abs(y4) + abs(y5)
  P3 <- abs(y3) + abs(y4) + abs(y5)

  # Compute features according to levers in Afonso et al. (1999)
  FL1 <- MWI(P1)
  FL2 <- MWI(P2)
  FL4 <- MWI(P3)

  #--------------------------------------
  # Level 1 of Afonso et al. (1999)
  # Determine candidate R-peaks (many false positives)
  # Method: search for inflection points on P1
  # Make 2 copies of signal; shift them relative to each other
  # Determine rising and falling edges
  # Inflections are present in both copies of signal
  d <- sign(diff(FL1))
  d1 <- c(0,d)
  d2 <- c(d,0)
  f1 <- which(d1==1)
  f2 <- which(d2==-1)
  EventsL1=intersect(f1,f2) #Detected events

  #--------------------------------------
  # Level 2 of Afonso et al. (1999)
  # Eliminate some false positives resulting from level 1
  # Do this by comparing the ouptput of level 1 against 2 different
  # thresholds, a low one and a high one. Each putative beat is classified
  # as signal or noise. This procedure results in 2 'channels'; one with
  # many beats classified as signals and few beats as noise (channel 1), and
  # one with many beats classified as noise and few as signals (channel 2).
  # In level 3, some logic is then applied to the two channels

  # Check number of L1 events and compute mean level
  lEventsL1 <- length(EventsL1)
  if (lEventsL1 <= 0) return (NULL)
  meanL1 <- sum(FL2[EventsL1]) / length(EventsL1)

  # Determine noise and signal levels and set detection thresholds
  NL <- meanL1 - meanL1*0.1
  SL <- meanL1 + meanL1*0.1
  threshold1 <- 0.08                 #Threshold detection channel 1
  threshold2 <- 0.7                  #Threshold detection channel 2

  chan1 <- detectionblock(FL2, EventsL1, NL, SL, threshold1)
  chan2 <- detectionblock(FL2, EventsL1, NL, SL, threshold2)

  #--------------------------------------
  # Level 3 of Afonso et al. (1999)
  # This level fuses the beat detection status from each of the 2 one-channel
  # detection algorithms in level 2 by incorporating a set of if-then-else rules.
  # The rules incorporate the fact that the 2 one-channel detection blocks have
  # complementary detection rates. There are four possible cases to design
  # rules for. If both channels indicate a beat then the output of level 3
  # classifies the current event as a beat. Etc. See Afonso et al (1999), p 197,
  # for a complete description.

  classL3 <- 0
  for (i in 1:lEventsL1) {

    c1 <- chan1$class[i]
    c2 <- chan2$class[i]

    # See Table at left bottom of p. 195
    # (programmed a little differently than in nqrsdetect.m)

    if (c1 == 1 & c2 == 1) classL3 <- c(classL3, 1)     # Classification as signal

    if (c1 == 0 & c2 == 1) classL3 <- c(classL3, 1)     # Classification as signal, should not occur

    if (c1 == 0 & c2 == 0) classL3 <- c(classL3, 0)     # Classification as noise

    if (c1 == 1 & c2 == 0) {
      delta1 <- (chan1$ds[i] - threshold1) / (1-threshold1)
      delta2 <- (threshold2 - chan2$ds[i]) / threshold2
      if (delta1 > delta2) {
        classL3 <- c(classL3, 1)   # Classification as Signal
      } else {
        classL3 <- c(classL3, 0)   # Classification as Noise
      }
    }

  }

  l <- length(classL3)
  if (l > 1) {
    classL3 <- classL3[2:l]
  }
  signalL3 <- EventsL1[which(classL3 > 0)]  # Signal Level 3
  noiseL3 <- EventsL1[which(classL3 == 0)]  # Noise Level 3

  #--------------------------------------
  # Level 4 of Afonso et al. (1999)
  # This level incorporates another one-channel detection block and uses feature
  # input to the MWI. This level reduces FN’s (events which were inaccurately
  # missed as beats by level 3). The beat detection rates after level 3 are
  # higher than those from the detections in level 2.

  # Lower threshold; calculate initial signal and noise levels
  threshold <- 0.3
  ls <- length(signalL3)
  ln <- length(noiseL3)
  SL <- NL <- 1
  if (ls > 0) SL <- sum(FL4[signalL3]) / ls
  if (ln > 0) NL <- sum(FL4[noiseL3]) / ln

  signalL4 <- 0
  noiseL4 <- 0
  dsL4 <- 0
  for (i in 1:lEventsL1) {

    evt <- EventsL1[i]

    if (classL3[i]==1) {   #Classification after Level 3 as Signal
      signalL4 <- c(signalL4, EventsL1[i])
      SL <- history(SL, FL4[evt])
      ds <- (FL4[evt] - NL) / (SL - NL)
      ds <- max(ds, 0)
      ds <- min(ds, 1)
      dsL4 <- c(dsL4, ds)
    } else {                          #Classification after Level 3 as Noise
      ds=(FL4[evt] - NL) / (SL - NL)
      ds <- max(ds, 0)
      ds <- min(ds, 1)
      dsL4 <- c(dsL4, ds)
      if (ds > threshold) {           #new classification as Signal
        signalL4 <- c(signalL4, evt)
        SL <- history(SL, FL4[evt])
      } else {                        #new classification as Noise
        noiseL4 <- c(noiseL4, evt)
        NL <- history(NL, FL4[evt])
      }
    }
  }
  ls <- length(signalL4)
  ln <- length(noiseL4)
  lds <- length(dsL4)
  if (ls > 1) signalL4 <- signalL4[2:ls]
  if (ln > 1) noiseL4 <- noiseL4[2:ln]
  if (lds > 1) dsL4 <- dsL4[2:lds]

  #--------------------------------------
  # Level 5 of Afonso et al. (1999)
  # The previous levels do not incorporate any timing information in the decision
  # logic. Level 5, thus, includes decision logic to eliminate possible false
  # detection during the refractory period.
  #
  # NOTE GvB: I do not think that Levels 5 and 6 of nqrsdetect.m actually
  # implement Level 5 of Afonso et al. (1999). There is no Level 6 in Afsonso.
  # I *THINK* that levels 5 and 6 of nqrsdetect.m are intended to correct
  # long and short R-R intervals. I have tried to implement the logic of
  # nqrsdetect.m here, separately for long and short R-R intervals, but adapted
  # somewhat because I suspected some bugs in that code (maybe I did not completely
  # understand it)

  signalL5 <- signalL4
  noiseL5 <- noiseL4
  periods <- diff(signalL4)

  # Compute a running mean of the RR intervals
  m1 <- 50     #100 in nqrsdetect.m
  a <- 1
  b <- rep(1/(m1), m1)
  meanperiod <- as.vector(signal::filter(b, a, periods))
  # use mean of (51:end) as value for first 50
  # (not present in nqrsdetect.m)
  lmp <- length(meanperiod)
  if (lmp > m1) meanperiod[1:m1] <- rep(mean(periods[1:m1]), m1)

  # First, detect R-R intervals that are too long (> 1.5 * running mean)
  # If such an interval is found, then it is possible that a beat was missed.
  # Check the noiseL4 array in the interval between the two consecutive signal
  # beats and check the Detection Strength value against a lower threshold.

  # Lower threshold; calculate initial signal and noise levels
  threshold <- 0.2
  ls <- length(signalL4)
  ln <- length(noiseL4)
  SL <- NL <- 1
  if (ls > 0) SL <- sum(FL4[signalL4]) / ls
  if (ln > 0) NL <- sum(FL4[noiseL4]) / ln

  lp <- length(periods)
  sq <- NULL
  if(lp > 0) {sq <- (1:lp)}
  for (i in sq) {
    # nqrsdetect.m did not involve an [i] index on meanperiod (bug??)
    if (periods[i] > meanperiod[i]*1.5) {     # if RR-interval is to long
      # check signal array
      interval <- seq(signalL4[i], signalL4[i+1])
      # is a noise beat present in that interval?
      critical <- intersect(interval, noiseL4)
      # if so, then go back and test that noise beat against a lower threshold
      lc <- length(critical)
      sq2 <- NULL
      if (lc > 0) {sq2 <- (1:lc)}
      for (j in sq2) {
        ds <- (FL4[critical[j]] - NL) / (SL - NL)
        if (ds > threshold) {         #Classification as Signal
          signalL5 <- union(signalL5, critical[j])
          noiseL5 <- noiseL5[-critical[j]]
        }
      }
    }
  }
  # union in R does not sort vector as Matlab/Octave does
  signalL5 <- sort(signalL5)

  # Next, detect intervals that are too short (< 0.5 * running mean)
  # compute new periods because they may have been adapted in the
  # preceding step (long R-R intervals deleted). If a short interval
  # is detected, then the amplitude of that peak and of its neighbour at FL2
  # are compared to the mean FL2 peak amplitude. The peak that is farthest away
  # from the mean peak level is deleted.
  # GvB 25-apr-2012: changed algorithm. Removing points within the loop leads
  # to problems with NA's. Make vector with values to be removed but remove them
  # only after the loop.
  periods <- diff(signalL5)
  lp <- length(periods)
  meanperiod <- as.vector(signal::filter(b, a, periods))
  lmp <- length(meanperiod)
  if (lmp > m1) meanperiod[1:m1] <- rep(mean(periods[1:m1]), m1)
  level <- mean(FL2[signalL5])
  sq <- NULL
  if (lp > 0) sq <- (1:lp)
  rm <- 0    #vector for points to be removed
  for (i in sq) {
    if (periods[i] < meanperiod[i]*0.5) {     #if RR-interval is to short
      d1 <- abs(FL2[signalL5[i]] - level)
      d2 <- abs(FL2[signalL5[(i+1)]] - level)
      if (d1 > d2) {
        rm <- c(rm, (i+1))
      } else {
        rm <- c(rm, i)
      }
    }
  }
  # now remove beats in the rm vector
  if (length(rm) > 1) {
    signalL5 <- signalL5[-rm[2:length(rm)]]
    noiseL5 <- sort(union(noiseL5,rm))
  }

  #-------------------------------------
  # Now, finally, convert the signal peaks at level 5
  # to the original ECG signal (not downsampled)
  peaks <- conversion(ecg, FL2, signalL5, M, N, fs)

  # # Added GvB: find local maxima in unfiltered ECG
  # # use an interval of 50 ms around the candidate peak, half at each side
  # half.int <- (round(0.05*fs/2))
  # cand.pks <- sort(peaks)
  # peaks <- NULL
  # for (i in 1:length(cand.pks)) {
  #   peaks <- c(peaks, which.max(S[(cand.pks[i]-half.int+1):(cand.pks[i]+half.int)]) + cand.pks[i]-half.int)
  # }
  return(peaks)
}

#------------------------------------------------------------------------------
# downsample - downsample a signal without interpolation
#              similar to MatlaB/Octave function
#
# parameters:
# S:          signal to downsample
# rate:       downsample factor
# phase:      starting sample number
#
# returns:    downsampled signal
#
# Note: it is usually better do use the functions 'resample' or 'decimate'
# from library(signal), as these functions will handle correct anti-aliasing
# downsampling and interpolation.
#
# 20120411 Geert van Boxtel
#------------------------------------------------------------------------------
downsample <- function (S, rate, phase=1) {
  down <- seq(phase, length(S), rate);
  return(S[down]);
}

#------------------------------------------------------------------------------
# polyphase - polyphase implementation of decimation filters
#
# y=polyphase(S,h,M)
#
# INPUT
#   S       ecg signal data
#   h       filter coefficients
#   M       downsampling rate
#
# OUTPUT
#   filtered signal
#
# DEPENDS
#   library(signal), function downsample
#
# ported to R from Matlab code in nqrsdetect by Geert van Boxtel 4-apr-2012;
# also included some comments
#
# test polyphase:
# S <- seq(1:1230305)
# h <- c(1,2,3,4,5,6,7,8,9,10,0,0,3,2,5,4,3,2,1,0)
# M<-10
# non_polyphase<- function (S,h,M) {
#   Sf <- signal::filter(h,1,S)
#   y <- signal::resample(Sf, M)
#   return (y)
# }
# system.time(y_nonp <- non_polyphase(S,h,M))
# system.time(y_poly <- polyphase(S,h,M))
# which(y_nonp!=y_poly)
#------------------------------------------------------------------------------

polyphase <- function (S, h, M) {

  # Determining polyphase components ek
  # Select from filter coefficients h
  ncols <- ceiling(length(h) / M)
  e <- matrix(0, M, ncols)
  l <- 1
  m <- length(h) %% M
  while (m > 0) {
    el <- array(0, ncols)
    for (n in 1:ncols) {
      el[n] <- h[M*(n-1)+l]
    }
    e[l,1:ncols] <- el[1:ncols]
    l <- l + 1
    m <- m - 1
  }
  for (i in l:M) {
    ncols2 <- floor(length(h)/M)
    el <- array(0,ncols2)
    for (n in 1:ncols2) {
      el[n] <- h[M*(n-1)+i]
    }
    e[i,1:ncols2] <- el[1:ncols2]
  }

  # Downsampling and filtering
  mx <- ceiling((length(S) + M) / M)
  Sdelay <- S
  w <- matrix (0, M, mx)
  for (i in 1:M) {
    Sd <- downsample(Sdelay, M)
    a <- signal::filter(e[i,], 1, Sd)
    if (length(a) < mx) a <- c(a, rep(0, (mx-length(a))))
    w[i, 1:mx] <- a
    Sdelay <- c(rep(0,i), S)
  }
  y <- colSums(w)
  return(y[1:(length(y)-1)])
}



#------------------------------------------------------------------------------
# MWI - Moving window integrator, computes the mean of two samples
#   y=MWI(S)
#
# INPUT
#   S       Signal
#
# OUTPUT
#   y       output signal
#------------------------------------------------------------------------------
MWI <- function (S) {

  a <- c(0, S)
  b <- c(S, 0)
  y <- (a+b)/2
  return(y[1:(length(y)-1)])

}

#------------------------------------------------------------------------------
# detectionblock - computation of one detection block
#
#   ret <- detectionblock(mwi,events, NL, SL, threshold)
#
# INPUT
#   mwi         Output of the MWI
#   events      Events of Level 1 (see Afonso et al)
#   NL          Initial Noise Level
#   SL          Initial Signal Level
#   threshold   Detection threshold (between [0,1])
#
# OUTPUT (as list)
#   signal      Events which are computed as Signal
#   noise       Events which are computed as Noise
#   ds          Detection strength of the Events
#   class       Classification: 0=noise, 1=signal
#------------------------------------------------------------------------------

detectionblock <- function (mwi, events, NL, SL, threshold) {

  # Initialize variables (initial 0 is stripped off later)
  signal <- 0
  noise <- 0
  ds <- 0
  class <- 0

  # Initial signal and noise levels
  sum.signal <- SL
  sum.noise <- NL

  # Start loop over events
  for (i in 1:length(events)) {
    evt <- events[i]

    # Compute detection strength and limit it to range (0,1)
    detstr <- (mwi[evt] - NL) / (SL - NL)
    detstr <- max(detstr, 0)
    detstr <- min(detstr, 1)
    ds <- c(ds, detstr)

    if (detstr > threshold) {
      # Classification as signal
      signal <- c(signal, evt)
      class <- c(class, 1)
      sum.signal <- sum.signal + mwi[evt]
      # Update the signal level
      # Note: difference with Matlab implementation:
      # length(signal) is taken as the divisor here - not length(signal)+1
      # This is because of the initial 0 in the R implementation
      SL <- sum.signal / (length(signal) + 0)
    } else {
      # Classification as noise
      noise <- c(noise, evt)
      class <- c(class, 0)
      sum.noise <- sum.noise + mwi[evt]
      # Update the noise level
      # Note: difference with Matlab implementation
      # length(noise) is taken as the divisor here - not length(noise)+1
      # This is because of the initial 0 in the R implementation
      NL <- sum.noise / (length(noise) + 0)
    }
  }

  # Return variables as a list
  return(list (signal= signal[2:length(signal)],
               noise = noise[2:length(noise)],
               ds    = ds[2:length(ds)],
               class = class[2:length(class)]))
}

#------------------------------------------------------------------------------
#  history - computes y[n]=(1-lambda)*x[n]+lambda*y[n-1]
#
#   yn=history(ynm1,xn)
#------------------------------------------------------------------------------
history <- function (ynm1, xn) {
  lambda <- 0.8 #forgetting factor
  return((1 - lambda)*xn + lambda*ynm1)
}


#------------------------------------------------------------------------------
#
# conversion - sets the fiducial points of the downsampled Signal on the
# samplepoints of the original Signal
#
#   [pnew]=conversion(S,FL2,pold,M,N,fs)
#
# INPUT
#   S           Original ECG Signal
#   FL2         Feature of Level 2 [1]
#   pold        old fiducial points
#   M           M downsampling rate
#   N           filter order
#   fs          sample rate
#
# OUTPUT
#   pnew        new fiducial points
#
#------------------------------------------------------------------------------

conversion <- function (S, FL2, pold, M, N, fs) {

  signaln <- pold
  p <- M
  q <- 1
  FL2res <- signal::resample(FL2,p,q)    #Upsampling of FL2
  nans <- which(is.nan(S))
  S[nans] <- mean(S)             #Replaces NaNs in Signal (necessary for filtering)

  # Convert beats to original sampling rate
  signaln1 <- array(0, length(signaln))
  for (i in 1:length(signaln)) {
    signaln1[i] <- signaln[i] + (M-1) * (signaln[i] - 1)
  }

  # Set the fiducial points on the maximum of (upsampled) FL2
  signaln2 <- signaln1
  range <- 1:length(FL2res)
  signaln3 <- array(0, length(signaln))
  for (i in 1:length(signaln2)) {
    start <- max((signaln2[i] - M), 1)
    stp <- min((signaln2[i] + M), length(FL2res))
    interv <- start:stp
    FL2int <- FL2res[interv]
    pk <- which.max(FL2int)
    if (length(pk) == 0) pk <- signaln2[i] - start
    else pk <- pk[1]
    delay <- N/2 + M
    signaln3[i] <- pk + signaln2[i] - M - delay
  }

  #Set the fiducial points on the maximum of the original signal
  Bw <- 5.6
  Bwn <- 1 / (fs/2) * Bw
  Wn <- c(Bwn, 5*Bwn)
  N1 <- 32
  b <- signal::fir1(N1, Wn,'pass')
  Sf <- signal::filtfilt(b, S)     #Filtered Signal with bandwidth 5.6-28 Hz
  beg <- round(1.5*M)
  fin <- 1*M
  signaln4 <- array(0, length(signaln))
  for (i in 1:length(signaln3)) {
    start <- max(signaln3[i]-beg, 1)
    stp <- min(signaln3[i]+fin, length(Sf))
    interv <- start:stp
    Sfint <- Sf[interv] - mean(Sf[interv])
    pk <- which.max(Sfint)
    if (length(pk) == 0) pk <- signaln3[i]-start
    else pk <- pk[1]
    signaln4[i] <- pk + signaln3[i] - beg - 1
  }
  signal <- signaln4

  # Delete first and last points because of initial transient of the filter
  # in polyphase_imp
  cutbeginning <- which(signal < N)
  fpointsb <- signal[cutbeginning]
  cutend <- which(signal > length(S)-N)
  fpointse <- signal[cutend]
  # GvB 14-may-2012
  if (length(fpointsb) > 0 & length(fpointse) > 0) rpeaks <- signal [-c(fpointsb, fpointse)]
  else rpeaks <- signal
  return(rpeaks)
}
