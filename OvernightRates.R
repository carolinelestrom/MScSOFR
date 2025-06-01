####################################################################################################################
####################################################################################################################
#----------------------------------------- Data --------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
SOFR <- read.csv("~/Documents/KU/MSc/Data/SOFR.csv")
IORB <- read.csv("~/Documents/KU/MSc/Data/IORB.csv")
TargetL <- read.csv("~/Documents/KU/MSc/Data/DFEDTARL.csv")
TargetU <- read.csv("~/Documents/KU/MSc/Data/DFEDTARU.csv")
IORBSample <- IORB %>%
  filter(row_number() %in% c(1:1401))
SOFRSample <- SOFR %>%
  filter(row_number() %in% c(305:1305))
TargetLSample <- TargetL %>%
  filter(row_number() %in% c(426:1826))
TargetUSample <- TargetU %>%
  filter(row_number() %in% c(426:1826))
SOFRSample$observation_date <- as.Date(SOFRSample$observation_date)
TargetLSample$observation_date <- as.Date(TargetLSample$observation_date)
TargetUSample$observation_date <- as.Date(TargetUSample$observation_date)
#SOFRSample$SOFR <- SOFRSample$SOFR/100


### Generate a complete sequence of dates from start to end
AllDates <- seq(min(SOFRSample$observation_date), max(SOFRSample$observation_date), by = "day")

### Create full data frame with all dates
SOFRSampleWeekend <- data.frame(observation_date = AllDates) %>%
  left_join(SOFRSample, by = "observation_date") %>%
  arrange(observation_date)

### Fill NAs forward to replace missing values (including weekends)
SOFRSampleWeekend <- SOFRSampleWeekend %>%
  mutate(SOFR = zoo::na.locf(SOFR, na.rm = FALSE))



nrow(IORBSample)             ### 1363
nrow(SOFRSampleWeekend)      ### 1363
nrow(TargetUSample)          ### 1363
nrow(TargetLSample)          ### 1363


### Define the start and end dates
SampleStart <- as.Date("2021-07-29")
SampleEnd <- as.Date("2025-05-29")

### Create a sequence of daily dates
SampleTime <- seq.Date(from = SampleStart, to = SampleEnd, by = "day")

Sample <- data.frame("IORB" = IORBSample$IORB,
                     "SOFR" = SOFRSampleWeekend$SOFR,
                     "TargetL" = TargetLSample$DFEDTARL,
                     "TargetU" = TargetUSample$DFEDTARU,
                     "Date" = SampleTime)

### Specified jump dates
SampleJumps <- data.frame(Jump = as.Date(c("2021-07-29", "2021-09-23", "2021-11-04", "2021-12-16",
                                           "2022-01-27", "2022-03-17", "2022-05-05", "2022-06-16", "2022-07-28", "2022-09-22", "2022-11-03", "2022-12-15",
                                           "2023-02-02", "2023-03-23", "2023-05-04", "2023-06-15", "2023-07-27", "2023-09-21", "2023-11-02", "2023-12-24",
                                           "2024-02-01", "2024-03-21", "2024-05-02", "2024-06-13", "2024-08-01", "2024-09-19", "2024-11-08", "2024-12-19",
                                           "2025-01-30", "2025-03-20", "2025-05-08"))
                          , Type = "Jumps")




### Plot it, plot it real good
ggplot(Sample, aes(x = Date)) +
  geom_line(aes(y = SOFR, color = "SOFR"), size = 1.3) +
  geom_line(aes(y = IORB, color = "IORB"), size = 1.3) +
  geom_line(aes(y = TargetL, color = "Target Range"), linetype = "dotted", size = 1.3) +
  geom_line(aes(y = TargetU, color = "Target Range"), linetype = "dotted", size = 1.3) +
  #geom_vline(xintercept = as.numeric(SampleJumps), linetype = "dashed", color = "black", alpha = 0.7) +
  geom_vline(data = SampleJumps, aes(xintercept = Jump, color = Type), linetype = "dashed", alpha = 0.7) +
  labs(title = "Overnight Rates",
       x = "Time", y = "") +
  scale_color_manual(name = "Rates", values = c("SOFR" = "#901a1E", "IORB" = "#39641c", "Target Range" = "#666666", "Jumps" = "black")) +
  scale_x_date(
    breaks = as.Date(c("2021-04-01", "2021-08-01", "2021-12-01", "2022-04-01", "2022-08-01", "2022-12-01",
                       "2023-04-01", "2023-08-01", "2023-12-01", "2024-04-01", "2024-08-01", "2024-12-01", "2025-04-01")),
    labels = c("2021-04", "2021-08", "2021-12", "2022-04", "2022-08", "2022-12",
               "2023-04", "2023-08", "2023-12", "2024-04", "2024-08", "2024-12", "2025-04"),
    limits = c(as.Date("2021-07-29"), as.Date("2025-04-21"))
  ) +
  #theme(legend.position = "none") +
  theme_minimal()


####################################################################################################################
####################################################################################################################
#--------------------------------------- Final Plot ----------------------------------------------------------------
####################################################################################################################
####################################################################################################################
ggplot(Sample, aes(x = Date)) +
  geom_line(aes(y = SOFR, color = "SOFR"), size = 1.7) +
  geom_line(aes(y = IORB, color = "IORB"), size = 1.7) +
  geom_line(aes(y = TargetL, color = "Target Range"), linetype = "dotted", size = 1.7) +
  geom_line(aes(y = TargetU, color = "Target Range"), linetype = "dotted", size = 1.7) +
  #geom_vline(xintercept = as.numeric(SampleJumps), linetype = "dashed", color = "black", alpha = 0.7) +
  geom_vline(data = SampleJumps, aes(xintercept = Jump, color = Type), linetype = "dashed", alpha = 0.7, size = 1) +
  labs(title = "Overnight Rates",
       x = "Time", y = "") +
  scale_x_date(
    breaks = as.Date(c("2021-04-01", "2021-08-01", "2021-12-01", "2022-04-01", "2022-08-01", "2022-12-01",
                       "2023-04-01", "2023-08-01", "2023-12-01", "2024-04-01", "2024-08-01", "2024-12-01", "2025-04-01")),
    labels = c("2021-04", "2021-08", "2021-12", "2022-04", "2022-08", "2022-12",
               "2023-04", "2023-08", "2023-12", "2024-04", "2024-08", "2024-12", "2025-04"),
    limits = c(as.Date("2021-07-29"), as.Date("2025-04-21"))
  ) +
  scale_color_manual(name = "Rates", values = c("Target Range" = "#4D5D53",
                                                "SOFR" = "#901a1E", 
                                                "N-ME" = "#FFBCD9", 
                                                "IORB" = "#39641c", 
                                                "S/ME" = "#FF8C00",
                                                "Jumps" = "#4D5D53")) +
  #theme(legend.position = "none") +
  theme_minimal() %+replace%
  theme(
    plot.title = element_text(size = 47, hjust = 0),
    axis.title = element_text(size = 37),
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    legend.position = "top",
    legend.text = element_text(size = 23),
    legend.title = element_text(size = 27),
    legend.background = element_rect(fill = "white", color = "black", size = 0.7),
    legend.box.background = element_rect(color = "black", size = 0.7)
  )





ggplot(Sample, aes(x = Date)) +
  geom_line(aes(y = SOFR, color = "SOFR"), size = 1.7) +
  geom_line(aes(y = IORB, color = "IORB"), size = 1.7) +
  geom_line(aes(y = TargetL, color = "Target Range"), linetype = "dashed", size = 1.7) +
  geom_line(aes(y = TargetU, color = "Target Range"), linetype = "dashed", size = 1.7) +
  #geom_vline(xintercept = as.numeric(SampleJumps), linetype = "dashed", color = "black", alpha = 0.7) +
  geom_vline(data = SampleJumps, aes(xintercept = Jump), color = "#4D5D53", linetype = "dashed", alpha = 1, size = 1) +
  labs(title = "Overnight Rates",
       x = "Time", y = "") +
  scale_x_date(
    breaks = as.Date(c("2022-01-01", "2023-01-01", "2024-01-01", "2025-01-01")),
    labels = c("2022-01-01", "2023-01-01", "2024-01-01", "2025-01-01"),
    limits = c(as.Date("2021-07-29"), as.Date("2025-05-29"))
  ) +
  scale_color_manual(name = "Rates", values = c("Target Range" = "#39641c",
                                                "SOFR" = "#901a1E", 
                                                "N-ME" = "#FFBCD9", 
                                                "IORB" = "#FF8C00", 
                                                "S-ME" = "#FF8C00",
                                                "Jumps" = "#4D5D53")) +
  #theme(legend.position = "none") +
  theme_minimal() %+replace%
  theme(
    plot.title = element_text(size = 47, hjust = 0.5, vjust = 1.5),
    axis.title = element_text(size = 37),
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    legend.position = "top",
    legend.text = element_text(size = 23),
    legend.title = element_text(size = 27),
    legend.background = element_rect(fill = "white", color = "black", size = 0.7),
    legend.box.background = element_rect(color = "black", size = 0.7)
  )


####################################################################################################################
####################################################################################################################
#----------------------------------- Color Pallette ----------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Black/White ish
"#4D5D53", "#F0EAD6"
### Bielefeld Vibes
"#FF8C00", "#404080", "#87A96B", "#FFBCD9", "#996666"
### KU
#"#901a1E", "#39641c", "#666666", "#ffbd38", "#0a5963", "#122947", "#425570"
### Nice
#"hotpink3", "orchid4", "steelblue"




