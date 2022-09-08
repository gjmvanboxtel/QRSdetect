# bahoura.R - QRS detection by wavelet transform as proposed by Bahoura, Hassani & Hubin (1997)
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
# 20190211  GvB       Initial setup for package QRSdetect
# 20220907  GvB       Use package gsignal instead of signal
#
#---------------------------------------------------------------------------------------------------------------------

#' Bahoura et al. QRS detection
#'
#' Detect R peaks of the QRS complex in a raw ECG record, based on the
#' algorithm proposed by Bahoura, Hassani & Hubin (1997)
#'
#' This function attempts to detect ECG fiducial points in a single-channel electrocardiogram signal using the
#' algorithm proposed by Bahoura, Hassani & Hubin (1997), which is a simplified version of the
#' method proposed by Li, Zheng & Tai (1995).
#'
#' @param ecg The input single-channel input vector (raw ECG)
#' @param fs The frequency in Hz with which the ecg was sampled
#'
#' @return Numeric array containing the indices (sample numbers) at which the fiducial R-peaks were found
#'
#' @references Bahoura, M., Hassani, M., & Hubin, M. (1997). DSP implementation of wavelet transform for real time ECG
#'             wave forms detection and heart rate analysis. Computer Methods and Programs in Biomedicine, 52, 35-44.
#'             DOI: \href{https://doi.org/10.1016/S0169-2607(97)01780-X}{10.1016/S0169-2607(97)01780-X}.\cr
#'             Li, C., Zheng, C., & Tai, C. (1995). Detection of ECG Characteristic Points Using Wavelet Transforms.
#'             IEEE Transactions on Biomedical Engineering, 42(1), 21-28.
#'             DOI: \href{https://doi.org/10.1109/10.362922}{10.1109/10.362922}.
#'
#' @author Geert van Boxtel
#'
#' @examples
#' data(rec100)
#' fs <- 360
#' pks <- bahoura(rec100$MLII, fs)
#'
#' \dontrun{
#' # plot first 5 seconds of data
#' N <- 5 * fs
#' plot (rec100$time[1:N], rec100$MLII[1:N], type = "l", main = "MIT-BIH database, record 100",
#'       xlab = "Time (s)", ylab = "Amplitude (mV)")
#' points (pks[which(pks<=N)]/fs, rec100$MLII[pks[which(pks<=N)]], col="red")
#' }
#' @export

bahoura <- function (ecg, fs) {

  # pad signal with one second to beginning and end of data to ramp up and down the filters
  pecg <- c(rev(ecg[1:fs]), ecg, rev(ecg[(length(ecg) - fs + 1):length(ecg)]))

  # impulse responses (Bahoura et al p 44)
  h <- c(1/8, 3/8, 3/8, 1/8)
  g <- c(2, -2)

  # use signal resampled to 250 Hz, which was also used by Bahoura et al.
  R <- gsignal::resample(pecg, 25000, round(fs * 100))

  # compute wavelet and scaling coefficients for first 4 levels
  S1 <- gsignal::conv(R, h)
  W1 <- gsignal::conv(R, g)
  R <- gsignal::resample(R, 12500, 25000)

  S2 <- gsignal::conv(R, h)
  W2 <- gsignal::conv(R, g)
  R <- gsignal::resample(R, 6250, 12500)

  S3 <- gsignal::conv(R, h)
  W3 <- gsignal::conv(R, g)
  R <- gsignal::resample(R, 3125, 6250)

  S4 <- gsignal::conv(R, h)
  W4 <- gsignal::conv(R, g)

  # Compute modulus maxima lines
  getMML <- function (W, j) {

    # find minima and maxima
    thr <- RMS(W)
    max.idx <- which(diff(sign(diff(W)))==-2)+1
    max.idx <- max.idx[which(W[max.idx] > thr)]
    min.idx <- which(diff(sign(diff(W)))==+2)+1
    min.idx <- min.idx[which(W[min.idx] < -thr)]
    idx <- sort(c(max.idx,min.idx))
    pks <- W[idx]

    # max QRS peak width = 120 ms; fs depends on scale
    mw <- 0.120 * (250 / (2 ^ (j - 1)))
    mml <- NULL
    l <- length(idx)
    i <- 2
    while (i <= l) {
      # pos and neg peaks should be closer together than 120 ms, and peaks should be positive - negative
      if (abs(idx[i] - idx[(i-1)]) <= mw && sign(pks[i]) == -1 && sign(pks[(i-1)]) == 1) {
        mml <- c(mml, c(idx[(i-1)], idx[i]))
        i <- i + 2
      } else {
        i <- i + 1
      }
    }
    # correct for filtering delay
    return(mml - floor((2 ^ j - 1) / 2))
  }

  MML1 <- getMML(W1, 1)
  MML2 <- getMML(W2, 2)
  MML3 <- getMML(W3, 3)
  MML4 <- getMML(W4, 4)

  # This function finds the concordance between modulus maxima lines (with a tolerance)
  intersect.tol <- function (x, y, tol) {
    d <- outer(x, y, FUN=function(x, y) abs(x - y))
    m <- which(d <= tol, arr.ind = TRUE)
    # rows are x, columns are y
    rc <- c(x[m[,1]], y[m[,2]])
    return(unique(rc[1:(length(rc) / 2)]))
  }

  # upsampling means multiplying by two
  peaks <- intersect.tol(MML3, 2 * MML4, tol=10)
  peaks <- intersect.tol(MML2, 2 * peaks, tol=10)
  peaks <- intersect.tol(MML1, 2 * peaks, tol=10)

  # The R-wave is the zero-crossing of the final modulus maxima lines on wavelet scale 1
  peaks <- matrix(peaks, ncol=2, byrow=TRUE)
  idx <- apply(peaks, 1, function(x) x[1] + which.max(abs(W1[x[1]:x[2]])) - 1)

  # upsample to original fs
  idx <- round(idx * (fs / 250))

  # compute local maxima in the unfiltered ECG
  # use an interval of 100 ms around the candidate peak, half at each side
  half.int <- (round(0.1 * fs / 2))
  cand.idx <- sort(idx)
  idx <- NULL
  for (i in 1:length(cand.idx)) {
    idx <- c(idx, which.max(pecg[(cand.idx[i]-half.int+1):(cand.idx[i]+half.int)]) + cand.idx[i]-half.int)
  }

  return(idx[idx > fs & idx < length(pecg) - fs + 1] - fs)

}

