#' Canadian Incident COVID-19 Cases
#'
#' This dataset contains 3+ years of incident COVID-19 case counts
#' as reported by [opencovid.ca](https://opencovid.ca/) as of
#' July 4, 2023.
#'
#' @format
#' A data frame with 1,253 rows and 2 columns:
#' \describe{
#'   \item{date}{The observed date. A `Date` object.}
#'   \item{incident_cases}{The number of new recorded cases for this date.}
#' }
#'
#' @source This data is available under the CC-BY-4.0 License.
#'
#' See:
#'
#' Berry, I., O’Neill, M., Sturrock, S. L., Wright, J. E., Acharya, K.,
#'   Brankston, G., Harish, V., Kornas, K., Maani, N., Naganathan, T., Obress, L.,
#'   Rossi, T., Simmons, A. E., Van Camp, M., Xie, X., Tuite, A. R., Greer, A. L.,
#'   Fisman, D. N., & Soucy, J.-P. R. (2021). A sub-national real-time
#'   epidemiological and vaccination database for the COVID-19 pandemic in
#'   Canada. Scientific Data, 8(1).
#'   \doi{10.1038/s41597-021-00955-2}
#'
"cancovid"
