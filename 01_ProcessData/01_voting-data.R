library(jsonlite)
url <- "https://static01.nyt.com/elections-assets/2020/data/api/2020-11-03/national-map-page/national/president.json"
jdata <- read_json(url, simplifyDataFrame = TRUE)

voting <- jdata$data$races$counties %>%
  map(~select(.x, fips, votes, margin2020, last_updated, results)) %>%
  map(~mutate(.x, trump = results$trumpd, biden = results$biden)) %>%
  map(~select(.x, -results)) %>%
  invoke(rbind, .) %>%
  mutate(pct_rep_lead = 100*(trump - biden)/votes,
         rep_lead_ntile = ntile(pct_rep_lead, 4))

saveRDS(voting, "../00_RawData/other/voting-agg.rds")
