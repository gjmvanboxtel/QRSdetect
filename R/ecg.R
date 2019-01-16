#' Electrocardiogram data
#'
#' A data set containing the raw electrocardiogram data from MIT-BIH database
#' record 100, The data were collected with a sampling rate of 360 Hz and
#' digitized with 11 bit resolution over a \eqn{\pm } 5 mV range.
#'
#' @format A dataset with 3 variables and 108000 observations (the first 5 minutes out of the
#'         recording of about 30 minutes):
#' \describe{
#'   \item{time}{the time base in seconds (0 - 299.997)}
#'   \item{MLII}{ECG from a Modified Leg lead, in millivolts}
#'   \item{V5}{ECG from a V5 lead, in millivolts}
#' }
#'
#' @source \url{https://physionet.org/physiobank/database/mitdb/}, (doi:10.13026/C2F305).
#'
#' @references Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database. IEEE Eng in Med and Biol
#'             20(3):45-50 (May-June 2001). (PMID: 11446209)\cr
#'             Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG, Mietus JE, Moody GB, Peng C-K,
#'             Stanley HE. PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex
#'             Physiologic Signals. Circulation 101(23):e215-e220, 2000 (June 13).
#'
"ecg"

# Code used to export dataset:
# dfn <- paste('/media/geert/DATA/Data/MIT-BIH Arrhythmia Database/100-dat.txt',sep='')
# header <- try(scan(dfn, nlines = 1, what = character(), quiet=T), silent=T)
# if (class(header) == "try-error") cat(paste("\nError scanning header of data file, rec=",rec))
# dat <- try(read.table(dfn, skip = 2, header=F), silent = T)
# if (class(dat) == "try-error") cat(paste("\nError reading data file, rec=",rec))
# names(dat) <- header[2:4]
# ecg <- dat[1:108000,]    # 5 min * 60 s * 360 Hz
# devtools::use_data(ecg)
