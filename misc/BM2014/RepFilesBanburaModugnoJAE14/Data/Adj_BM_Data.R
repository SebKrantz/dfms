library(fastverse)

Data <- readxl::read_xlsx("~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/Data2009-10-15.xlsx",
                         skip = 1, .name_repair = "none", sheet = 'MonthlyLong')
nam <- names(Data)
DataQ <- readxl::read_xlsx("~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/Data2009-10-15.xlsx",
                          skip = 1, .name_repair = "none", sheet = 'QuarterlyLong')
namQ <- names(DataQ)
namall <- c(nam, namQ)

TransfM <- readxl::read_xlsx("~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/Legend.xlsx",
                           skip = 0, .name_repair = "none", sheet = 'TransfM')
all(TransfM$`Name short` %in% nam)
TransfM %<>% fsubset(`Name short` %in% nam)

TransfQ <- readxl::read_xlsx("~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/Legend.xlsx",
                             skip = 0, .name_repair = "none", sheet = 'TransfQ')
all(TransfQ$`Name short` %in% namQ)
TransfQ %<>% fsubset(`Name short` %in% namQ)

Models <- readxl::read_xlsx("~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/Legend.xlsx",
                             skip = 0, .name_repair = "none", sheet = 'Models')
all(Models$Series %in% namall)
Models %<>% fsubset(Series %in% namall)

writexl::write_xlsx(list(TransfM = TransfM, TransfQ = TransfQ, Models = Models),
                    path = "~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/Legend.xlsx")


writexl::write_xlsx(X, path = "~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/X.xlsx")
X_diff <- fscale(fdiff(X)) %>% ss(-1L)

writexl::write_xlsx(X_diff, path = "~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/X_diff.xlsx")


X_diff <- readxl::read_xlsx("~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/X_diff.xlsx") # %>% make.names()

DFM(X_diff, 2)

# Generate Complete Datasets ------------------------------------------------------------------------------

library(janitor)
### Monthly Dataset
prep_data <- function(Data) {
  BM14_M <- Data %>% clean_names()
  anyDuplicated(names(BM14_M))
  names(BM14_M)[1L] <- "date"
  vlabels(BM14_M) <- BM14_M[1L, ] %>% as.character() %>% iconv("", "UTF-8") %>% # tolower() %>%
     trimws() %>% gsub("  ", " ", .) %>% tools::toTitleCase()
  vlabels(BM14_M, "code") <- BM14_M[2L, ] %>% as.character() %>% # iconv("", "UTF-8") %>% # tolower() %>%
     trimws() %>% gsub("  ", " ", .)
  BM14_M %<>% ss(-(1:2))
  settfmv(BM14_M, 1L, as.Date)
  settfmv(BM14_M, -1L, function(x) copyMostAttrib(as.numeric(x), x))
  attr(BM14_M$date, "label") <- "Date"
  BM14_M %<>% ss(pcountNA(.) < fncol(.)-1L)
  BM14_M
}
BM14_M <- prep_data(Data)
### Quarterly Dataset
BM14_Q <- prep_data(DataQ)
### Models
BM14_Models <- Models %>% gv(2:5) %>% clean_names()
settfm(BM14_Models, series = make_clean_names(series))
BM14_Models %<>% ss(!pallNA(gv(.,-1)))
BM14_Models %<>% tfmv(-1, function(x) as.logical(replace_NA(x, 0)))
settfm(BM14_Models, freq = iif(series %in% names(BM14_M), "M", "Q"))
setcolorder(BM14_Models, c("series", "freq"))
BM14_Models %>% sbt(freq == "Q", series) %>% unlist() %in% names(BM14_Q)
# Subsetting Monthly and Quarterly Data: Only Series used in Models
BM14_M %<>% gv(c("date", sbt(BM14_Models, freq == "M", series)[[1]]))
BM14_Q %<>% gv(c("date", sbt(BM14_Models, freq == "Q", series)[[1]]))
# Adding Transformation Info and series labels + codes
TransfM %<>% clean_names()
settfm(TransfM, name_short = make_clean_names(name_short))
TransfQ %<>% clean_names()
settfm(TransfQ, name_short = make_clean_names(name_short))
TransfM %>% merge(namlab(BM14_M) %>% rnm(Variable = name_short), by = "name_short") # %>% View() # Same
BM14_Models$label <- c(vlabels(BM14_M), vlabels(BM14_Q))[BM14_Models$series] %>% unname()
BM14_Models$code <- c(vlabels(BM14_M, "code"), vlabels(BM14_Q, "code"))[BM14_Models$series] %>% unname() # No clean names everywhere
BM14_Models$log_trans <- c(with(TransfM, setNames(transf_log, name_short)),
                           with(TransfQ, setNames(transf_log, name_short)))[BM14_Models$series] %>% as.logical()
setcolorder(BM14_Models, c("series", "label", "code", "freq", "log_trans"))

# Saving all
save(BM14_Models, BM14_M, BM14_Q, file = "~/Documents/R/DFM/misc/BM2014/RepFilesBanburaModugnoJAE14/RepFilesBanburaModugnoJAE14/Data/BM14_prepared.RData")
.c(BM14_Models, BM14_M, BM14_Q) %=% list(qDF(BM14_Models), qDF(BM14_M), qDF(BM14_Q))
usethis::use_data(BM14_Models, BM14_M, BM14_Q, overwrite = TRUE)

