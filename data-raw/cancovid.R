## code to prepare `cancovid` dataset goes here
## run on 4 July 2023
res <- httr::GET(
  url = "https://api.opencovid.ca/timeseries?",
  query = list('fmt' = 'csv', 'geo' = 'can', 'stat' = 'cases')
)
cancovid <- httr::content(res)
first_case <- min(which(cancovid$value > 0))
if (first_case > 1) {
  cancovid <- cancovid[-seq(first_case - 1L), ]
}
library(dplyr)
cancovid <- cancovid %>%
  select(date, incident_cases = value_daily)

usethis::use_data(cancovid, overwrite = TRUE)
