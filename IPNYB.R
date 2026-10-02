# 1. Persiapan --------------------------------------------------------------

library(readxl)
library(dplyr)
library(sf)
library(ggplot2)
library(spdep)
library(spatialreg)

alpha <- 0.05
project_dir <- "D:/1. Skripshit/Skripsi Bismillah/Code/R"
output_dir <- file.path(project_dir, "IPM")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

judul <- function(x) cat("\n---", x, "---\n")


# 2. Data -------------------------------------------------------------------

raw <- read_excel(
  file.path(project_dir, "Data Skripsi.xlsx"),
  sheet = "Data_Noprov",
  skip = 2,
  col_names = FALSE
)

data <- raw[, c(1, 14, 11, 8, 9, 6, 12, 23)]
names(data) <- c(
  "kab_kota", "ipm", "persen_miskin", "tpt", "tpak",
  "kepadatan", "laju_pdrb_adhk", "AIR"
)

data <- data %>%
  mutate(
    kab_kota = trimws(as.character(kab_kota)),
    across(-kab_kota, as.numeric),
    nama_join = case_when(
      grepl("Kep.*Seribu", kab_kota, ignore.case = TRUE) ~
        "Administrasi Kepulauan Seribu",
      grepl("^Jakarta (Barat|Pusat|Selatan|Timur|Utara)$", kab_kota) ~
        paste("Kota Administrasi", kab_kota),
      TRUE ~ trimws(sub("^Kabupaten\\s+", "", kab_kota))
    )
  )

formula_model <- ipm ~ persen_miskin + tpt + tpak + kepadatan +
  laju_pdrb_adhk + AIR
variabel_model <- all.vars(formula_model)

baris_valid <- complete.cases(data[, variabel_model]) &
  Reduce(`&`, lapply(data[, variabel_model], is.finite))
data <- data[baris_valid, ]

provinsi_jawa <- c(
  "DKI Jakarta", "Banten", "Jawa Barat", "Jawa Tengah",
  "Daerah Istimewa Yogyakarta", "Jawa Timur"
)

shp <- st_read(
  file.path(project_dir, "Peta Jawa", "Peta_KabKota_Jawa.shp"),
  quiet = TRUE
) %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_make_valid() %>%
  filter(WADMPR %in% provinsi_jawa, !is.na(WADMKK))

gagal_join <- anti_join(
  data %>% select(kab_kota, nama_join),
  st_drop_geometry(shp) %>% distinct(WADMKK),
  by = c("nama_join" = "WADMKK")
)

if (nrow(gagal_join) > 0L) {
  print(gagal_join)
  stop("Ada wilayah yang gagal di-join.")
}

peta <- shp %>%
  left_join(data, by = c("WADMKK" = "nama_join")) %>%
  filter(!is.na(ipm))

if (nrow(peta) != nrow(data) || anyDuplicated(peta$kab_kota)) {
  stop("Join data dan shapefile tidak satu-ke-satu.")
}

peta$.region_id <- paste(peta$WADMPR, peta$WADMKK, sep = "::")
model_data <- as.data.frame(st_drop_geometry(peta))
rownames(model_data) <- peta$.region_id


# 3. Matriks bobot spasial --------------------------------------------------

tambah_link <- function(nb, nama, a, b) {
  if (sum(nama == a) != 1L || sum(nama == b) != 1L) {
    stop("Wilayah penghubung tidak ditemukan secara unik.")
  }

  i <- match(a, nama)
  j <- match(b, nama)
  nb[[i]] <- sort(unique(c(nb[[i]][nb[[i]] > 0L], as.integer(j))))
  nb[[j]] <- sort(unique(c(nb[[j]][nb[[j]] > 0L], as.integer(i))))
  nb
}

sf_use_s2(FALSE)
nb <- suppressWarnings(poly2nb(
  peta,
  row.names = peta$.region_id,
  queen = TRUE,
  snap = 1e-5
))

link_tambahan <- data.frame(
  dari = c("Bangkalan", "Administrasi Kepulauan Seribu"),
  ke = c("Kota Surabaya", "Kota Administrasi Jakarta Utara")
)

for (k in seq_len(nrow(link_tambahan))) {
  nb <- tambah_link(
    nb,
    peta$WADMKK,
    link_tambahan$dari[k],
    link_tambahan$ke[k]
  )
}

if (any(card(nb) == 0L) || n.comp.nb(nb)$nc != 1L ||
    !is.symmetric.nb(nb, force = TRUE)) {
  stop("SWM Queen belum valid atau belum terhubung.")
}

W <- nb2listw(nb, style = "W", zero.policy = TRUE)


# 4. Eksplorasi dan autokorelasi spasial -----------------------------------

deskriptif <- t(sapply(model_data[, variabel_model], function(x) {
  c(
    N = length(x),
    Mean = mean(x),
    SD = sd(x),
    Min = min(x),
    Median = median(x),
    Max = max(x)
  )
}))

judul("Statistik deskriptif")
print(round(deskriptif, 3))

ipm_spasial <- setNames(peta$ipm, peta$.region_id)

moran_global <- moran.test(
  ipm_spasial,
  W,
  alternative = "two.sided",
  zero.policy = TRUE,
  spChk = TRUE
)

judul("Global Moran's I IPM")
print(moran_global)

lisa <- localmoran_perm(
  ipm_spasial,
  W,
  nsim = 999,
  alternative = "two.sided",
  zero.policy = TRUE,
  iseed = 123
)

p_sim <- "Pr(z != E(Ii)) Sim"
if (!p_sim %in% colnames(lisa)) {
  stop("Kolom p-value simulasi dua sisi Local Moran tidak ditemukan.")
}

z_ipm <- setNames(as.numeric(scale(peta$ipm)), peta$.region_id)
lag_z_ipm <- lag.listw(W, z_ipm, zero.policy = TRUE)

peta$lisa_I <- lisa[, "Ii"]
peta$lisa_Z <- lisa[, "Z.Ii"]
peta$lisa_p <- lisa[, p_sim]
peta$lisa_p_fdr <- p.adjust(peta$lisa_p, method = "BH")
peta$lisa_cluster <- factor(
  case_when(
    peta$lisa_p >= alpha ~ "Tidak Signifikan",
    z_ipm >= 0 & lag_z_ipm >= 0 ~ "High-High",
    z_ipm < 0 & lag_z_ipm < 0 ~ "Low-Low",
    z_ipm >= 0 & lag_z_ipm < 0 ~ "High-Low",
    TRUE ~ "Low-High"
  ),
  levels = c(
    "High-High", "Low-Low", "High-Low", "Low-High", "Tidak Signifikan"
  )
)

judul("LISA IPM")
print(table(peta$lisa_cluster))
print(
  peta %>%
    st_drop_geometry() %>%
    filter(lisa_p < alpha) %>%
    select(WADMKK, WADMPR, ipm, lisa_I, lisa_Z, lisa_p, lisa_p_fdr, lisa_cluster) %>%
    arrange(lisa_p)
)


# 5. Pemodelan spasial ------------------------------------------------------

model_ols <- lm(formula_model, data = model_data, na.action = na.fail)

judul("OLS")
print(summary(model_ols))

judul("VIF")
print(car::vif(model_ols))

uji_bp <- lmtest::bptest(model_ols)
judul("Breusch-Pagan")
print(uji_bp)

moran_residual_ols <- lm.morantest(
  model_ols,
  W,
  alternative = "two.sided",
  zero.policy = TRUE,
  spChk = TRUE
)
judul("Moran's I residual OLS")
print(moran_residual_ols)

uji_lm <- lm.RStests(
  model_ols,
  W,
  test = "all",
  zero.policy = TRUE,
  spChk = TRUE
)
judul("LM/RS Lag, Error, dan Robust LM/RS")
print(summary(uji_lm))

model_sar <- lagsarlm(
  formula_model,
  data = model_data,
  listw = W,
  zero.policy = TRUE,
  na.action = na.fail
)

model_sem <- errorsarlm(
  formula_model,
  data = model_data,
  listw = W,
  zero.policy = TRUE,
  na.action = na.fail
)

judul("SAR")
print(summary(model_sar))

judul("SEM")
print(summary(model_sem))

uji_durbin <- SD.RStests(
  model_ols,
  W,
  test = "all",
  Durbin = TRUE,
  zero.policy = TRUE
)
judul("Uji kebutuhan SDM/SDEM")
print(summary(uji_durbin))

perlu_sdm <- uji_durbin$SDM_adjRSWX$p.value < alpha
perlu_sdem <- uji_durbin$SDEM_RSWX$p.value < alpha
model_spasial <- list(SAR = model_sar, SEM = model_sem)

if (perlu_sdm) {
  model_spasial$SDM <- lagsarlm(
    formula_model,
    data = model_data,
    listw = W,
    Durbin = TRUE,
    zero.policy = TRUE,
    na.action = na.fail
  )

  judul("SDM")
  print(summary(model_spasial$SDM))
}

if (perlu_sdem) {
  model_spasial$SDEM <- errorsarlm(
    formula_model,
    data = model_data,
    listw = W,
    Durbin = TRUE,
    zero.policy = TRUE,
    na.action = na.fail
  )

  judul("SDEM")
  print(summary(model_spasial$SDEM))
}


# 6. Evaluasi dan pemilihan model ------------------------------------------

for (nama in names(model_spasial)) {
  judul(paste("LR parameter spasial", nama))
  print(LR1.Sarlm(model_spasial[[nama]]))

  judul(paste("Wald parameter spasial", nama))
  print(Wald1.Sarlm(model_spasial[[nama]]))
}

if (perlu_sdm) {
  judul("LR: SDM vs SAR")
  print(LR.Sarlm(model_spasial$SDM, model_sar))
}

if (perlu_sdem) {
  judul("LR: SDEM vs SEM")
  print(LR.Sarlm(model_spasial$SDEM, model_sem))
}

uji_moran_residual <- function(model) {
  if (inherits(model, "lm")) {
    return(lm.morantest(
      model,
      W,
      alternative = "two.sided",
      zero.policy = TRUE,
      spChk = TRUE
    ))
  }

  moran.test(
    setNames(as.numeric(residuals(model)), peta$.region_id),
    W,
    alternative = "two.sided",
    zero.policy = TRUE,
    spChk = TRUE
  )
}

semua_model <- c(list(OLS = model_ols), model_spasial)

diagnostik_residual <- bind_rows(lapply(names(semua_model), function(nama) {
  uji <- uji_moran_residual(semua_model[[nama]])
  data.frame(
    Model = nama,
    Moran_I_Residual = unname(uji$estimate[1]),
    P_value_Moran_Residual = uji$p.value
  )
}))

perbandingan_model <- data.frame(
  Model = names(semua_model),
  LogLik = vapply(semua_model, function(x) as.numeric(logLik(x)), numeric(1)),
  AIC = vapply(semua_model, AIC, numeric(1)),
  BIC = vapply(semua_model, BIC, numeric(1))
) %>%
  left_join(diagnostik_residual, by = "Model") %>%
  arrange(AIC)

judul("Perbandingan model")
print(perbandingan_model, row.names = FALSE, digits = 6)

kandidat_final <- perbandingan_model %>%
  filter(P_value_Moran_Residual >= alpha)

if (nrow(kandidat_final) == 0L) {
  warning("Semua model masih memiliki autokorelasi residual; dipilih AIC terendah.")
  kandidat_final <- perbandingan_model
}

nama_model_final <- kandidat_final$Model[which.min(kandidat_final$AIC)]
model_final <- semua_model[[nama_model_final]]

judul(paste("Model final:", nama_model_final))
if (inherits(model_final, "Sarlm")) {
  print(summary(model_final, Nagelkerke = TRUE))
} else {
  print(summary(model_final))
}

judul("Moran's I residual model final")
print(uji_moran_residual(model_final))

residual_final <- as.numeric(residuals(model_final))
residual_final <- residual_final[is.finite(residual_final)]

if (length(residual_final) < 3L || length(residual_final) > 5000L ||
    !is.finite(sd(residual_final)) || sd(residual_final) == 0) {
  stop("Residual model final tidak memenuhi syarat diagnostik normalitas.")
}

uji_shapiro_final <- shapiro.test(residual_final)
z_residual_final <- as.numeric(scale(residual_final))
keputusan_normalitas <- if (uji_shapiro_final$p.value >= alpha) {
  "Tidak cukup bukti untuk menolak normalitas residual"
} else {
  "Normalitas residual ditolak"
}

judul("Uji normalitas residual model final")
print(uji_shapiro_final)
cat("Keputusan:", keputusan_normalitas, "\n")


# 7. Dampak dan prediksi SDM ------------------------------------------------

imp_sdm <- NULL
metrik_prediksi_sdm <- NULL

if ("SDM" %in% names(model_spasial)) {
  set.seed(123)
  imp_sdm <- impacts(model_spasial$SDM, listw = W, R = 1000)

  judul("Spatial impacts SDM")
  print(summary(imp_sdm, zstats = TRUE, short = TRUE))

  if (!identical(peta$.region_id, rownames(model_data))) {
    stop("Urutan wilayah peta dan data model tidak identik.")
  }

  id_w <- attr(W$neighbours, "region.id")
  if (is.null(id_w) || !identical(as.character(id_w), peta$.region_id)) {
    stop("Urutan wilayah pada W tidak identik dengan data model.")
  }

  # Respons IPM tidak diberikan agar prediksi dibentuk dari kovariat dan W.
  data_prediksi_sdm <- model_data[
    , setdiff(names(model_data), "ipm"), drop = FALSE
  ]
  rownames(data_prediksi_sdm) <- rownames(model_data)

  prediksi_sdm_obj <- predict(
    model_spasial$SDM,
    newdata = data_prediksi_sdm,
    listw = W,
    pred.type = "TC",
    zero.policy = TRUE,
    spChk = TRUE
  )

  id_prediksi <- attr(prediksi_sdm_obj, "region.id")
  if (is.null(id_prediksi)) id_prediksi <- names(prediksi_sdm_obj)
  if (is.null(id_prediksi)) {
    stop("ID wilayah tidak tersedia pada hasil prediksi SDM.")
  }

  indeks_prediksi <- match(peta$.region_id, id_prediksi)
  nilai_prediksi <- as.numeric(prediksi_sdm_obj)

  if (anyNA(indeks_prediksi) ||
      length(nilai_prediksi) != length(id_prediksi) ||
      any(!is.finite(nilai_prediksi[indeks_prediksi]))) {
    stop("Hasil prediksi SDM tidak lengkap atau tidak sejajar dengan peta.")
  }

  peta$prediksi_sdm <- nilai_prediksi[indeks_prediksi]
  peta$residual_prediksi_sdm <- peta$ipm - peta$prediksi_sdm

  metrik_prediksi_sdm <- data.frame(
    RMSE = sqrt(mean(peta$residual_prediksi_sdm^2)),
    MAE = mean(abs(peta$residual_prediksi_sdm)),
    R2_Korelasi = cor(peta$ipm, peta$prediksi_sdm)^2
  )

  judul("Akurasi prediksi SDM pada sampel estimasi")
  print(metrik_prediksi_sdm, row.names = FALSE, digits = 5)
} else {
  cat("\n[INFO] Spatial impacts dan prediksi SDM dilewati karena SDM tidak diestimasi.\n")
  warning(
    sprintf(
      "Visualisasi SDM tidak dibuat (p = %.4g; alpha = %.2f).",
      uji_durbin$SDM_adjRSWX$p.value,
      alpha
    ),
    call. = FALSE
  )
}


# 8. Visualisasi ------------------------------------------------------------

tema_peta <- theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "grey15"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    plot.caption = element_text(hjust = 0, color = "grey35"),
    plot.margin = margin(10, 12, 10, 12)
  )

tema_grafik <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    plot.caption = element_text(hjust = 0, color = "grey35"),
    panel.grid.minor = element_blank()
  )

simpan_plot <- function(plot, nama_file, lebar, tinggi) {
  ggsave(
    filename = file.path(output_dir, nama_file),
    plot = plot,
    width = lebar,
    height = tinggi,
    dpi = 600,
    bg = "white"
  )
}

# 8.1 Peta jaringan Queen

indeks_link <- data.frame(
  dari = match(link_tambahan$dari, peta$WADMKK),
  ke = match(link_tambahan$ke, peta$WADMKK)
)

id_link_tambahan <- paste(
  pmin(indeks_link$dari, indeks_link$ke),
  pmax(indeks_link$dari, indeks_link$ke),
  sep = "-"
)

sisi_nb <- do.call(rbind, lapply(seq_along(nb), function(i) {
  j <- nb[[i]]
  j <- j[j > i]
  if (length(j) == 0L) return(NULL)
  data.frame(dari = rep.int(i, length(j)), ke = as.integer(j))
}))

sisi_nb$id <- paste(sisi_nb$dari, sisi_nb$ke, sep = "-")
sisi_nb$jenis <- factor(
  ifelse(
    sisi_nb$id %in% id_link_tambahan,
    "Koneksi tambahan",
    "Queen contiguity"
  ),
  levels = c("Queen contiguity", "Koneksi tambahan")
)
sisi_nb <- sisi_nb[order(sisi_nb$jenis), ]

# Web Mercator dipakai hanya untuk menggambar garis antarsentroid.
peta_jaringan <- st_transform(peta, 3857)
titik_jaringan <- suppressWarnings(st_point_on_surface(peta_jaringan))
koordinat_titik <- st_coordinates(titik_jaringan)

geometri_sisi <- st_sfc(
  lapply(seq_len(nrow(sisi_nb)), function(k) {
    st_linestring(
      koordinat_titik[c(sisi_nb$dari[k], sisi_nb$ke[k]), 1:2, drop = FALSE]
    )
  }),
  crs = st_crs(peta_jaringan)
)
sisi_nb_sf <- st_sf(sisi_nb, geometry = geometri_sisi)

peta_queen_modifikasi <- ggplot() +
  geom_sf(
    data = peta_jaringan,
    fill = "#F7F7F7",
    color = "#B8B8B8",
    linewidth = 0.15
  ) +
  geom_sf(
    data = sisi_nb_sf,
    aes(color = jenis, linewidth = jenis, alpha = jenis),
    lineend = "round"
  ) +
  geom_sf(
    data = titik_jaringan,
    shape = 21,
    size = 0.45,
    stroke = 0.15,
    color = "#404040",
    fill = "white"
  ) +
  scale_color_manual(
    values = c(
      "Queen contiguity" = "#35608D",
      "Koneksi tambahan" = "#D1495B"
    ),
    name = "Jenis hubungan",
    drop = FALSE
  ) +
  scale_linewidth_manual(
    values = c("Queen contiguity" = 0.22, "Koneksi tambahan" = 1.15),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c("Queen contiguity" = 0.35, "Koneksi tambahan" = 1),
    guide = "none"
  ) +
  labs(
    title = "Jaringan Ketetanggaan Queen Setelah Modifikasi",
    subtitle = "Kabupaten/kota di Pulau Jawa",
    caption = paste0(
      "Garis merah: Bangkalan-Kota Surabaya dan Kepulauan Seribu-",
      "Kota Administrasi Jakarta Utara. Garis menunjukkan topologi biner; ",
      "W untuk analisis distandardisasi per baris."
    )
  ) +
  tema_peta

simpan_plot(
  peta_queen_modifikasi,
  "IPNYB_Peta_Queen_Modifikasi.png",
  lebar = 11,
  tinggi = 6.5
)

# 8.2 Peta IPM dan LISA

peta_ipm <- ggplot(peta) +
  geom_sf(aes(fill = ipm), color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "C", direction = -1, name = "IPM") +
  labs(title = "Sebaran IPM Kabupaten/Kota di Pulau Jawa") +
  tema_peta

simpan_plot(peta_ipm, "IPNYB_Peta_IPM.png", lebar = 10, tinggi = 6)

peta_lisa <- ggplot(peta) +
  geom_sf(aes(fill = lisa_cluster), color = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c(
      "High-High" = "#B2182B",
      "Low-Low" = "#2166AC",
      "High-Low" = "#EF8A62",
      "Low-High" = "#67A9CF",
      "Tidak Signifikan" = "grey85"
    ),
    breaks = c(
      "High-High", "Low-Low", "High-Low", "Low-High", "Tidak Signifikan"
    ),
    drop = FALSE,
    name = "Klaster LISA"
  ) +
  labs(
    title = "Klaster LISA IPM Kabupaten/Kota di Pulau Jawa",
    subtitle = "Local Moran permutation (999 simulasi); Queen yang telah dimodifikasi",
    caption = paste0(
      "Klaster memakai p-value simulasi mentah < 0,05; ",
      "lisa_p_fdr menyimpan koreksi Benjamini-Hochberg."
    )
  ) +
  tema_peta

simpan_plot(peta_lisa, "IPNYB_Peta_LISA_IPM.png", lebar = 10, tinggi = 6)

# 8.3 Diagnostik normalitas residual

simpan_diagnostik_normalitas <- function(
    residual_z, nama_model, uji, keputusan, file_png) {
  png(
    file_png,
    width = 11,
    height = 5.8,
    units = "in",
    res = 600,
    bg = "white"
  )

  par_lama <- par(no.readonly = TRUE)
  on.exit({
    par(par_lama)
    dev.off()
  }, add = TRUE)

  par(
    mfrow = c(1, 2),
    mar = c(4.5, 4.5, 3, 1.2),
    oma = c(0, 0, 3.2, 0),
    las = 1,
    family = "sans"
  )

  hist(
    residual_z,
    breaks = "FD",
    probability = TRUE,
    col = "#BDD7E7",
    border = "white",
    main = "Histogram dan Kurva Kepadatan",
    xlab = "Residual terstandarisasi",
    ylab = "Kepadatan"
  )

  garis_x <- seq(min(residual_z), max(residual_z), length.out = 300)
  lines(density(residual_z), col = "#35608D", lwd = 2)
  lines(garis_x, dnorm(garis_x), col = "#D19A27", lwd = 2, lty = 2)
  legend(
    "topright",
    legend = c("Kepadatan residual", "Normal baku"),
    col = c("#35608D", "#D19A27"),
    lwd = 2,
    lty = c(1, 2),
    bty = "n",
    cex = 0.82
  )

  qqnorm(
    residual_z,
    pch = 19,
    cex = 0.65,
    col = adjustcolor("#35608D", alpha.f = 0.65),
    main = "Normal Q-Q Plot",
    xlab = "Kuantil teoretis normal",
    ylab = "Kuantil residual terstandarisasi"
  )
  qqline(residual_z, col = "#D19A27", lwd = 2, lty = 2)

  mtext(
    sprintf("Diagnostik Normalitas Residual Model Final (%s)", nama_model),
    side = 3,
    outer = TRUE,
    line = 1.65,
    font = 2,
    cex = 1.15
  )
  mtext(
    sprintf(
      "Shapiro-Wilk: W = %.4f; p-value = %.4g | %s",
      unname(uji$statistic),
      uji$p.value,
      tolower(keputusan)
    ),
    side = 3,
    outer = TRUE,
    line = 0.25,
    cex = 0.82,
    col = "grey30"
  )
}

simpan_diagnostik_normalitas(
  z_residual_final,
  nama_model_final,
  uji_shapiro_final,
  keputusan_normalitas,
  file.path(output_dir, "IPNYB_Diagnostik_Normalitas_Residual_Model_Final.png")
)

# 8.4 Visualisasi SDM

if ("SDM" %in% names(model_spasial)) {
  level_panel_sdm <- c("Aktual", "Prediksi SDM")

  peta_aktual_sdm <- peta
  peta_aktual_sdm$panel_sdm <- factor("Aktual", levels = level_panel_sdm)
  peta_aktual_sdm$nilai_ipm_sdm <- peta_aktual_sdm$ipm
  peta_aktual_sdm <- peta_aktual_sdm[, c("panel_sdm", "nilai_ipm_sdm")]

  peta_prediksi_sdm <- peta
  peta_prediksi_sdm$panel_sdm <- factor("Prediksi SDM", levels = level_panel_sdm)
  peta_prediksi_sdm$nilai_ipm_sdm <- peta_prediksi_sdm$prediksi_sdm
  peta_prediksi_sdm <- peta_prediksi_sdm[, c("panel_sdm", "nilai_ipm_sdm")]

  peta_banding_sdm <- rbind(peta_aktual_sdm, peta_prediksi_sdm)
  batas_fill_sdm <- range(peta_banding_sdm$nilai_ipm_sdm, finite = TRUE)

  peta_aktual_prediksi_sdm <- ggplot(peta_banding_sdm) +
    geom_sf(aes(fill = nilai_ipm_sdm), color = "white", linewidth = 0.15) +
    facet_wrap(~panel_sdm, nrow = 1) +
    scale_fill_viridis_c(
      option = "C",
      direction = -1,
      limits = batas_fill_sdm,
      name = "Nilai IPM"
    ) +
    labs(
      title = "IPM Aktual dan Prediksi Spatial Durbin Model (SDM)",
      subtitle = sprintf(
        "Prediksi TC | RMSE = %.3f | MAE = %.3f | R2 = %.3f",
        metrik_prediksi_sdm$RMSE,
        metrik_prediksi_sdm$MAE,
        metrik_prediksi_sdm$R2_Korelasi
      ),
      caption = paste0(
        "Skala warna kedua panel dibuat sama. Hasil ini adalah kalibrasi pada ",
        "sampel estimasi, bukan validasi out-of-sample."
      )
    ) +
    tema_peta +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", color = "grey15"),
      strip.background = element_rect(fill = "grey95", color = NA)
    )

  simpan_plot(
    peta_aktual_prediksi_sdm,
    "IPNYB_Peta_Aktual_vs_Prediksi_SDM.png",
    lebar = 14,
    tinggi = 5.8
  )

  data_actual_fitted <- peta %>%
    st_drop_geometry() %>%
    transmute(
      observasi = seq_len(n()),
      aktual = ipm,
      prediksi = prediksi_sdm
    )

  subtitle_sdm <- sprintf(
    "RMSE = %.3f | MAE = %.3f | R2 = %.3f",
    metrik_prediksi_sdm$RMSE,
    metrik_prediksi_sdm$MAE,
    metrik_prediksi_sdm$R2_Korelasi
  )

  scatter_aktual_prediksi <- ggplot(
    data_actual_fitted,
    aes(x = aktual, y = prediksi)
  ) +
    geom_point(size = 2.5, alpha = 0.75, color = "#35608D") +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.8,
      color = "#D1495B"
    ) +
    coord_equal() +
    labs(
      title = "Perbandingan IPM Aktual dan Prediksi SDM",
      subtitle = subtitle_sdm,
      x = "IPM Aktual",
      y = "IPM Prediksi",
      caption = "Garis putus-putus menunjukkan prediksi sempurna (Prediksi = Aktual)"
    ) +
    tema_grafik

  print(scatter_aktual_prediksi)
  simpan_plot(
    scatter_aktual_prediksi,
    "IPNYB_Scatter_Aktual_vs_Prediksi_SDM.png",
    lebar = 7,
    tinggi = 6
  )

  plot_titik_prediksi <- ggplot(data_actual_fitted, aes(x = observasi)) +
    geom_segment(
      aes(xend = observasi, y = aktual, yend = prediksi),
      color = "grey70",
      linewidth = 0.4
    ) +
    geom_point(aes(y = aktual, color = "Aktual"), size = 2) +
    geom_point(
      aes(y = prediksi, color = "Prediksi"),
      size = 2,
      shape = 4,
      stroke = 1
    ) +
    scale_color_manual(values = c("Aktual" = "#35608D", "Prediksi" = "#D1495B")) +
    labs(
      title = "Perbandingan IPM Aktual dan Prediksi SDM",
      subtitle = subtitle_sdm,
      x = "Observasi Kabupaten/Kota",
      y = "Indeks Pembangunan Manusia (IPM)",
      color = NULL
    ) +
    tema_grafik +
    theme(legend.position = "top")

  print(plot_titik_prediksi)

  plot_actual_fitted <- ggplot(data_actual_fitted, aes(x = observasi)) +
    geom_line(aes(y = prediksi, color = "Prediksi SDM"), linewidth = 0.8) +
    geom_point(aes(y = aktual, color = "Aktual"), size = 2, alpha = 0.8) +
    scale_color_manual(
      values = c("Aktual" = "#35608D", "Prediksi SDM" = "#D1495B")
    ) +
    labs(
      title = "Perbandingan IPM Aktual dan Prediksi SDM",
      subtitle = subtitle_sdm,
      x = "Observasi Kabupaten/Kota",
      y = "Indeks Pembangunan Manusia (IPM)",
      color = NULL,
      caption = paste0(
        "Titik menunjukkan nilai IPM aktual; ",
        "garis menunjukkan nilai prediksi Spatial Durbin Model."
      )
    ) +
    tema_grafik +
    theme(legend.position = "top")

  print(plot_actual_fitted)
  simpan_plot(
    plot_actual_fitted,
    "IPNYB_Actual_vs_Fitted_SDM.png",
    lebar = 11,
    tinggi = 6
  )
}
