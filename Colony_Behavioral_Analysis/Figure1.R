library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(patchwork)
library(purrr)
library(ggpubr)
library(stringr)

set.seed(42)

# ==============================================================================
# Read data
# ==============================================================================

AssayLarva <- read.csv("/Users/tomaskay/Desktop/OLDNeuropeptides/AssayLarva.csv")
AssayControl <- read.csv("/Users/tomaskay/Desktop/OLDNeuropeptides/AssayFood.csv")
Colony <- read.csv("/Users/tomaskay/Desktop/OLDNeuropeptides/ColonyExperiment.csv")

Colony[is.na(Colony)] <- 0

# Collapse manual annotations per ant
collapsed_AssayLarva <- AssayLarva %>%
  group_by(Ant) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )
collapsed_AssayLarva[is.na(collapsed_AssayLarva)] <- 0

collapsed_AssayControl <- AssayControl %>%
  group_by(Ant) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )
collapsed_AssayControl[is.na(collapsed_AssayControl)] <- 0

# ==============================================================================
# Colony context data
# ==============================================================================

behaviors <- c(
  "Foraging", "Nursing", "Guarding", "Incidental",
  "Manipulating", "Carrying", "Antennating"
)

Colony_long <- Colony %>%
  pivot_longer(
    cols = all_of(behaviors),
    names_to = "Behavior",
    values_to = "Count"
  ) %>%
  filter(!is.na(Count)) %>%
  mutate(
    Count = Count / 9,
    Behavior = factor(Behavior, levels = behaviors),
    Age = factor(Age, levels = c("12d", "4m"))
  )

Colony_summary <- Colony_long %>%
  group_by(Age, Behavior) %>%
  summarise(Total = sum(Count), .groups = "drop")

# ==============================================================================
# Colony plot 1: Nursing and foraging
# ==============================================================================

Colony_subset <- Colony_long %>%
  filter(Behavior %in% c("Nursing", "Foraging")) %>%
  mutate(Behavior = factor(Behavior, levels = c("Nursing", "Foraging")))

stat.test <- Colony_subset %>%
  group_by(Behavior) %>%
  t_test(Count ~ Age, var.equal = FALSE) %>%
  add_significance() %>%
  add_xy_position(x = "Behavior", dodge = 0.8)

p <- ggplot(Colony_subset, aes(x = Behavior, y = Count, fill = Age)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    colour = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.75,
    colour = "red",
    fatten = 2,
    position = position_dodge(width = 0.8),
    show.legend = FALSE
  ) +
  geom_jitter(
    aes(colour = Age),
    shape = 21,
    fill = "black",
    alpha = 0.4,
    size = 2,
    stroke = 0.3,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  ) +
  scale_fill_manual(values = c("12d" = "#019eff", "4m" = "#3700ff"), guide = "none") +
  scale_colour_manual(values = c("12d" = "black", "4m" = "black"), guide = "none") +
  labs(x = NULL, y = "Prop. time") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  stat_pvalue_manual(
    stat.test,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 1,
    size = 4,
    hide.ns = TRUE
  )

# ==============================================================================
# Colony plot 2: Other behaviors
# ==============================================================================

Colony_other <- Colony_long %>%
  filter(Behavior %in% c("Carrying", "Guarding", "Incidental", "Manipulating", "Antennating")) %>%
  mutate(
    Behavior = factor(
      Behavior,
      levels = c("Carrying", "Guarding", "Manipulating", "Antennating", "Incidental")
    )
  )

stat.other <- Colony_other %>%
  group_by(Behavior) %>%
  t_test(Count ~ Age, var.equal = FALSE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "Behavior", dodge = 0.8)

q <- ggplot(Colony_other, aes(x = Behavior, y = Count, fill = Age)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    colour = "black"
  ) +
  stat_summary(
    data = Colony_other,
    aes(x = Behavior, y = Count, group = Age),
    fun = median,
    geom = "crossbar",
    width = 0.7,
    colour = "red",
    fatten = 2,
    position = position_dodge(width = 0.8),
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_jitter(
    aes(colour = Age),
    shape = 21,
    fill = "black",
    alpha = 0.4,
    size = 2,
    stroke = 0.3,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  ) +
  scale_fill_manual(values = c("12d" = "#019eff", "4m" = "#3700ff"), guide = "none") +
  scale_colour_manual(values = c("12d" = "black", "4m" = "black"), guide = "none") +
  labs(x = NULL, y = "Prop. time") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  stat_pvalue_manual(
    stat.other,
    label = "p.adj.signif",
    tip.length = 0.01,
    bracket.size = 1,
    size = 4,
    hide.ns = TRUE
  ) +
  expand_limits(y = max(Colony_other$Count, na.rm = TRUE) * 1.12)

# Combined colony-context plot
(p | q) +
  plot_layout(widths = c(1, 2)) &
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) &
  theme(legend.position = "none")

# ==============================================================================
# Manual behavior data: prep
# ==============================================================================

behavior_cols <- c(
  "nocontact", "inspectingGrooming", "antennating",
  "carrying", "incidental", "guarding"
)

labels_map <- c(
  nocontact = "Contact",
  inspectingGrooming = "Manipulating",
  antennating = "Antennating",
  carrying = "Carrying",
  incidental = "Incidental",
  guarding = "Guarding"
)

pretty_levels <- c(
  "Contact", "Carrying", "Guarding",
  "Manipulating", "Antennating", "Incidental"
)

safe_prop <- function(x) {
  s <- sum(x, na.rm = TRUE)
  if (s == 0) rep(NA_real_, length(x)) else x / s
}

normalize_behavior_cols <- function(dat, cols) {
  miss <- setdiff(cols, names(dat))
  
  if (length(miss)) {
    dat[miss] <- 0
  }
  
  dat[cols <- intersect(cols, names(dat))] <- lapply(dat[cols], as.numeric)
  dat
}

collapsed_AssayLarva <- normalize_behavior_cols(collapsed_AssayLarva, behavior_cols)
collapsed_AssayControl <- normalize_behavior_cols(collapsed_AssayControl, behavior_cols)

# Tag age/treatment groups by row order
collapsed_AssayLarva <- collapsed_AssayLarva %>%
  mutate(Group = if_else(row_number() <= 12, "12d:Brood", "4m:Brood"))

collapsed_AssayControl <- collapsed_AssayControl %>%
  mutate(Group = if_else(row_number() <= 12, "12d:Food", "4m:Food"))

# Totals per group
totals_brood <- collapsed_AssayLarva %>%
  group_by(Group) %>%
  summarise(across(any_of(behavior_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

totals_food <- collapsed_AssayControl %>%
  group_by(Group) %>%
  summarise(across(any_of(behavior_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# Group-level tidy data
df_plot <- bind_rows(totals_brood, totals_food) %>%
  pivot_longer(cols = any_of(behavior_cols), names_to = "Behavior", values_to = "Count") %>%
  mutate(Behavior = labels_map[Behavior]) %>%
  complete(Behavior = pretty_levels, Group, fill = list(Count = 0)) %>%
  group_by(Group) %>%
  mutate(Proportion = safe_prop(Count)) %>%
  ungroup() %>%
  mutate(
    Proportion = if_else(Behavior == "Contact", 1 - Proportion, Proportion),
    Behavior = factor(Behavior, levels = pretty_levels),
    Group = factor(Group, levels = c("12d:Brood", "4m:Brood", "12d:Food", "4m:Food"))
  )

# Individual-level tidy data
collapsed_AssayLarva <- collapsed_AssayLarva %>%
  mutate(Ant = as.character(Ant))

collapsed_AssayControl <- collapsed_AssayControl %>%
  mutate(Ant = as.character(Ant))

df_individual <- bind_rows(collapsed_AssayLarva, collapsed_AssayControl) %>%
  pivot_longer(cols = any_of(behavior_cols), names_to = "Behavior", values_to = "Count") %>%
  mutate(Behavior = labels_map[Behavior]) %>%
  complete(Behavior = pretty_levels, Group, Ant, fill = list(Count = 0)) %>%
  group_by(Group, Ant) %>%
  mutate(Proportion = safe_prop(Count)) %>%
  ungroup() %>%
  mutate(
    Proportion = if_else(Behavior == "Contact", 1 - Proportion, Proportion),
    Behavior = factor(Behavior, levels = pretty_levels),
    Group = factor(Group, levels = c("12d:Brood", "4m:Brood", "12d:Food", "4m:Food"))
  )

# ==============================================================================
# Contact-only plot
# ==============================================================================

contact_df <- df_individual %>%
  filter(Behavior == "Contact")

fill_pal <- c(
  "12d:Brood" = "#019eff",
  "4m:Brood" = "#3700ff",
  "12d:Food" = "#FF9C00",
  "4m:Food" = "#FE6100"
)

pairs <- list(
  c("12d:Brood", "4m:Brood"),
  c("12d:Food", "4m:Food")
)

stat_contact <- map_dfr(pairs, function(g) {
  sub <- contact_df %>%
    filter(Group %in% g, !is.na(Proportion)) %>%
    droplevels()
  
  if (n_distinct(sub$Group) < 2 || min(table(sub$Group)) < 2) {
    return(NULL)
  }
  
  y.pos <- min(max(sub$Proportion, na.rm = TRUE) + 0.05, 0.98)
  
  t_test(sub, Proportion ~ Group, var.equal = FALSE) %>%
    add_significance() %>%
    mutate(
      group1 = g[1],
      group2 = g[2],
      y.position = y.pos
    )
})

p_contact <- ggplot(contact_df, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "black") +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.55,
    colour = "red",
    fatten = 2,
    show.legend = FALSE
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    shape = 21,
    size = 2,
    fill = "black",
    colour = "black",
    stroke = 0.25,
    alpha = 0.5
  ) +
  scale_fill_manual(values = fill_pal, guide = "none") +
  labs(x = NULL, y = "Prop. time") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

x_mid <- (length(levels(contact_df$Group)) + 1) / 2

p_contact <- p_contact +
  scale_x_discrete(labels = rep("", length(levels(contact_df$Group)))) +
  theme(
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 10, 24, 10),
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off") +
  annotate(
    "text",
    x = x_mid,
    y = 0,
    label = "Contact",
    vjust = 2.5,
    size = 4.2
  )

if (nrow(stat_contact) > 0) {
  p_contact <- p_contact +
    stat_pvalue_manual(
      stat_contact,
      label = "p.signif",
      hide.ns = FALSE,
      tip.length = 0.01,
      bracket.size = 1,
      size = 4
    )
}

# ==============================================================================
# Non-contact plots
# ==============================================================================

contact_df <- df_individual %>%
  filter(Behavior == "Contact") %>%
  filter(!is.na(Proportion), !is.na(Group))

noncontact_df <- df_individual %>%
  filter(Behavior != "Contact")

noncontact_12d <- noncontact_df %>%
  filter(str_starts(as.character(Group), "12d")) %>%
  droplevels() %>%
  mutate(Group = factor(Group, levels = c("12d:Brood", "12d:Food")))

noncontact_4m <- noncontact_df %>%
  filter(str_starts(as.character(Group), "4m")) %>%
  droplevels() %>%
  mutate(Group = factor(Group, levels = c("4m:Brood", "4m:Food")))

p_12d_noncontact <- ggplot(noncontact_12d, aes(x = Behavior, y = Proportion, fill = Group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.55,
    outlier.shape = NA,
    colour = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.55,
    colour = "red",
    fatten = 2,
    position = position_dodge(0.8),
    show.legend = FALSE
  ) +
  geom_jitter(
    aes(colour = Group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    shape = 21,
    size = 2,
    fill = "black",
    stroke = 0.25,
    alpha = 0.5
  ) +
  scale_fill_manual(values = fill_pal, guide = "none") +
  scale_colour_manual(
    values = c("12d:Brood" = "black", "12d:Food" = "black"),
    guide = "none"
  ) +
  labs(x = NULL, y = "Prop. time") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

p_4m_noncontact <- ggplot(noncontact_4m, aes(x = Behavior, y = Proportion, fill = Group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.55,
    outlier.shape = NA,
    colour = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.55,
    colour = "red",
    fatten = 2,
    position = position_dodge(0.8),
    show.legend = FALSE
  ) +
  geom_jitter(
    aes(colour = Group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    shape = 21,
    size = 2,
    fill = "black",
    stroke = 0.25,
    alpha = 0.5
  ) +
  scale_fill_manual(values = fill_pal, guide = "none") +
  scale_colour_manual(
    values = c("4m:Brood" = "black", "4m:Food" = "black"),
    guide = "none"
  ) +
  labs(x = NULL, y = "Prop. time") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# ==============================================================================
# Non-contact statistics and significance bars
# ==============================================================================

stat_noncontact_12d <- noncontact_12d %>%
  filter(Group %in% c("12d:Brood", "12d:Food"), !is.na(Proportion)) %>%
  group_by(Behavior) %>%
  group_modify(~ {
    sub <- .x
    
    if (n_distinct(sub$Group) < 2 || min(table(sub$Group)) < 2) {
      return(tibble())
    }
    
    y.pos <- min(max(sub$Proportion, na.rm = TRUE) + 0.05, 0.98)
    
    t_test(sub, Proportion ~ Group, var.equal = FALSE) %>%
      add_significance() %>%
      mutate(
        group1 = "12d:Brood",
        group2 = "12d:Food",
        y.position = y.pos
      )
  }) %>%
  ungroup()

stat_noncontact_4m <- noncontact_4m %>%
  filter(Group %in% c("4m:Brood", "4m:Food"), !is.na(Proportion)) %>%
  group_by(Behavior) %>%
  group_modify(~ {
    sub <- .x
    
    if (n_distinct(sub$Group) < 2 || min(table(sub$Group)) < 2) {
      return(tibble())
    }
    
    y.pos <- min(max(sub$Proportion, na.rm = TRUE) + 0.05, 0.98)
    
    t_test(sub, Proportion ~ Group, var.equal = FALSE) %>%
      add_significance() %>%
      mutate(
        group1 = "4m:Brood",
        group2 = "4m:Food",
        y.position = y.pos
      )
  }) %>%
  ungroup()

beh_levels <- levels(noncontact_12d$Behavior)
half_span <- 0.2

stat_noncontact_12d <- stat_noncontact_12d %>%
  mutate(
    Behavior = factor(Behavior, levels = beh_levels),
    xnum = as.numeric(Behavior),
    xmin = xnum - half_span,
    xmax = xnum + half_span,
    y.position = pmin(pmax(y.position, 0.02), 0.98)
  )

stat_noncontact_4m <- stat_noncontact_4m %>%
  mutate(
    Behavior = factor(Behavior, levels = beh_levels),
    xnum = as.numeric(Behavior),
    xmin = xnum - half_span,
    xmax = xnum + half_span,
    y.position = pmin(pmax(y.position, 0.02), 0.98)
  )

p_12d_noncontact <- p_12d_noncontact +
  stat_pvalue_manual(
    stat_noncontact_12d,
    label = "p.signif",
    hide.ns = FALSE,
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 1,
    size = 4
  ) +
  expand_limits(y = max(noncontact_12d$Proportion, na.rm = TRUE) + 0.08) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

p_4m_noncontact <- p_4m_noncontact +
  stat_pvalue_manual(
    stat_noncontact_4m,
    label = "p.signif",
    hide.ns = FALSE,
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 1,
    size = 4
  ) +
  expand_limits(y = max(noncontact_4m$Proportion, na.rm = TRUE)) +
  theme(legend.position = "none")

# ==============================================================================
# Final combined plot
# ==============================================================================

right_stack <-
  (p_12d_noncontact + theme(legend.position = "none")) /
  (p_4m_noncontact + theme(legend.position = "none"))

final_plot <- (p_contact + theme(legend.position = "none") | right_stack) +
  plot_layout(widths = c(1, 2))

final_plot
