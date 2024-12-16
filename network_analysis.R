# パッケージをロードする ----

# データの読み込みとクリーニング
library(tidyverse)
library(janitor)
# 希釈
library(iNEXT)
# ネットワーク図
library(igraph)
library(ggraph)
library(tidygraph)
library(ggrepel)

# OTUデータを読み込む ----

# 生データを読み込んで、列名（サンプル名）を標準化する
otu_raw_with_taxonomy <- read_delim(
  "data/otu_table.txt",
  delim = "\t", show_col_types = FALSE
) |>
  clean_names()

# 生データにはOTUのリード数と分類データが混合しているので、分ける
otu_taxonomy <- otu_raw_with_taxonomy |>
  select(otu_id, taxonomy)

otu_raw <- otu_raw_with_taxonomy |>
  select(-taxonomy)

# 植物のデータを読み込む ----

plants <- read_csv("data/plants.csv", show_col_types = FALSE) |>
  rename(
    sample = samplename,
    plant = Host.Plant,
    location = 採取地
  ) |>
  # サンプル名をOTUデータの標準化された名前に合わせる
  mutate(
    sample = make_clean_names(sample)
  )

# データクリーニング ----

# 闘値を設定する（次のステップに使う）
threshold_fraction <- 0.001

# 低頻度OTUを除外する
otu_clean <-
  otu_raw |>
  # 以下の処理のため、データの形を縦長に変える
  pivot_longer(names_to = "host", values_to = "abun", cols = -otu_id) |>
  # 闘値をサンプルあたり総リードの0.1%と設定し、それより少ないリード数を0とみなす
  mutate(threshold = sum(abun) * threshold_fraction, .by = otu_id) |>
  mutate(abun = case_when(
    abun < threshold ~ 0,
    .default = abun
  )) |>
  select(-threshold) |>
  # 全ての植物において個体数の合計が0のOTUを除外する
  mutate(total_otu_abun = sum(abun), .by = otu_id) |>
  filter(total_otu_abun > 0) |>
  select(-total_otu_abun) |>
  # 縦長から幅広い形に戻す
  pivot_wider(names_from = host, values_from = abun)

# 分類データを整える
tax_levels <- c(
  "kingdom", "phylum", "class", "order", "family", "genus", "species"
)

taxonomy_clean <-
  otu_taxonomy |>
  mutate(taxonomy = str_remove_all(taxonomy, ".__")) |>
  separate_wider_delim(
    taxonomy,
    delim = "; ",
    names = tax_levels,
    too_few = "align_start"
  )

# それぞれの分類階級がどれくらい同定できているのか確認する
taxonomy_clean |>
  pivot_longer(
    names_to = "rank",
    values_to = "taxon",
    -otu_id
  ) |>
  filter(!is.na(taxon)) |>
  count(rank, sort = TRUE)

# 希釈 ----

# 時間かかるので、ここではサンプルの最初の６つだけを解析する
otu_small <- otu_clean |>
  select(1:6)

inext_res <- iNEXT(
  column_to_rownames(otu_small, "otu_id"),
  q = 0,
  nboot = 20
)

# 実際のデータと推測された種数を比較する
# 1 に近ければ近いほど観察された種数と推測される種数が一致している
# 十分一致していれば、ネットワークなどで希釈する必要なし
cor(inext_res$AsyEst$Observed, inext_res$AsyEst$Estimator)

# 希釈プロットを作成して、保存する
ggiNEXT(inext_res, type = 1)

ggsave(filename = "results/rarefaction_curve.png")

# ネットワーク図 ----

# 綺麗なOTUデータに植物のデータと分類データを加える
otu_long_with_plants_and_taxonomy <-
  otu_clean |>
  pivot_longer(names_to = "sample", values_to = "abun", -otu_id) |>
  left_join(taxonomy_clean, by = "otu_id", relationship = "many-to-one") |>
  left_join(plants, by = "sample", relationship = "many-to-one")

# ネットワーク図の準備：データを用意する
# ここでは、科レベルにしている
network_comm <-
  otu_long_with_plants_and_taxonomy |>
  filter(!is.na(family)) |>
  filter(abun > 0) |>
  mutate(abun = 1) |>
  select(family, plant, abun) |>
  unique() |>
  pivot_wider(names_from = "family", values_from = abun, values_fill = 0)

# ネットワーク図の準備：dataframeをgraphに変換する
graph <- network_comm |>
  column_to_rownames("plant") |>
  as.matrix() |>
  graph_from_biadjacency_matrix(multiple = TRUE) |>
  as_tbl_graph() |>
  activate(nodes) |>
  # 図で示すデータを用意する
  # - 植物かどうか
  rename(plant = type) |>
  mutate(plant = !plant) |>
  # - 媒介中心性
  mutate(importance = centrality_betweenness()) |>
  # - 植物の名前
  mutate(plant_name = if_else(plant, name, NA_character_))

# ネットワーク図を作る
ggraph(graph) +
  geom_edge_link() +
  geom_node_point(aes(color = plant, size = importance)) +
  geom_label_repel(aes(label = plant_name, x = x, y = y))

ggsave(filename = "results/network.png")
