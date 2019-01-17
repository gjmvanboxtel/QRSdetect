# pantompkins.R - QRS detection algorithm proposed by Pan & Tompkins (1985)
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
# 20190116  GvB       Initial setup for package QRSdetect
#
#---------------------------------------------------------------------------------------------------------------------

#' Pan & Tompkins QRS detection
#'
#' Detect R peaks of the QRS complex in a raw ECG record, based on the
#' algorithm proposed by Pan & Tompkins (1985)
#'
#' This function attempts to detect ECG fiducial points in a single-channel electrocardiogram signal using the
#' algorithm proposed by Pan & Tompkins (1985). The algorithm was originally developed to be implemented in hardware.
#' Here a software implementation with digital filters, and some improvements:
#' \enumerate{
#'   \item Filter the signal 5-15 Hz (Butterworth); pre- and postpad the signal with 1 s of reversed data
#'         to ramp up and down the filter
#'   \item Compute the derivative
#'   \item Square the derivative
#'   \item Apply a moving average (0.15 s filter length)
#'   \item Search for the peak using findpeaks (mind = 200, minh = 2*sd of moving average), lookup local maximum in signal
#' }
#'
#' @param ecg The input single-channel input vector (raw ECG)
#' @param fs The frequency in Hz with which the ecg was sampled
#'
#' @return Numeric array containing the indices (sample numbers) at which the fiducial R-peaks were found
#'
#' @references Pan, J., & Tompkins, W.J. (1985). A real-time QRS detection algorithm. IEEE Transactions on
#'             Biomedical Engineering, Vol. BME-32, Issue 3, 230-236,
#'             DOI: \href{https://dx.doi.org/10.1109/TBME.1985.325532}{10.1109/TBME.1985.325532}
#'
#' @author Geert van Boxtel
#'
#' @examples
#' data(ecg)
#' fs <- 360
#' pks <- pantompkins(ecg$MLII, fs)
#'
#' \dontrun{
#' # plot first 5 seconds of data
#' N <- 5 * fs
#' plot (ecg$time[1:N], ecg$MLII[1:N], type = "l", main = "MIT-BIH database, record 100",
#'       xlab = "Time (s)", ylab = "Amplitude (mV)")
#' points (pks[which(pks<=N)]/fs, ecg$MLII[pks[which(pks<=N)]], col="red")
#' }
#' @export

pantompkins <- function (ecg, fs) {

  l <- length(ecg)

  ######################################################################
  # preprocessing: filter, derivative, squaring, moving average

  # Step 1: filter the input signal 5-15 Hz (use Butterworth, which has been shown to work best for ECG).
  # Pad the signal with 1 second of reversed data at the beginning and end to ramp up and down the filter,
  # then chop it off again after the filtering
  fecg <- signal::filtfilt(signal::butter(5,15*2/fs,"low"), signal::filtfilt(signal::butter(5,5*2/fs,"high"),
                           c(rev(utils::head(ecg,fs)),ecg,rev(utils::tail(ecg,fs)))))[(fs+1):(l+fs)]

  # Step 2: Compute the derivative
  decg <- array(0, l);
  for (i in 5:l) {
    decg[i] <- (2*fecg[i] + fecg[(i-1)] - fecg[(i-3)] - 2*fecg[(i-4)]) / 8;
  }

  # Step 3: Square the derivative
  secg <- decg^2;

  # Step 4: Apply a moving average
  # Use the filter function from the stats package for this purpose
  mecg <- stats::filter (secg, filter=rep(1/(0.15*fs),0.15*fs), method="convolution", sides=2);

  ######################################################################
  # compute fiducial marks

  # find peaks. A minimum dstance of 200 ms used because physiologically no RR interval can be shorter than 200 ms
  # 2* sd (mecg) seems a good value for minh
  peaks <- findpeaks (mecg, mind = round(0.2*fs), minh = 2*stats::sd(mecg, na.rm = TRUE))
  if (is.null(peaks) || length(peaks) <= 0) return (NULL)

  # compute local maxima in the unfiltered ECG
  # use an interval of 50 ms around the candidate peak, half at each side
  half.int <- (round(0.05*fs/2))
  cand.idx <- sort(peaks$idx)
  idx <- NULL
  for (i in 1:length(cand.idx)) {
    idx <- c(idx, which.max(ecg[(cand.idx[i]-half.int+1):(cand.idx[i]+half.int)]) + cand.idx[i]-half.int)
  }

  return (idx)
}

