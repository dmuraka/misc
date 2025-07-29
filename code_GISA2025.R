##############################################
##### sdmTMBパッケージによる時空間モデリング ######
##############################################

#install.packages(c("dplyr","ggplot2","sdmTMB","sf","sp"))
library(dplyr)
library(ggplot2)
library(sdmTMB)
library(sf)
library(sp)
library(visreg)

####### Pacific Cod（タラ）データの読み込み
pcod   <-data.frame(pcod[pcod$year>=2009, c("year","X","Y","depth","density","present")])
pcod[1:10,]

####### タラの在・不在のプロット
pcod_s <- st_as_sf(pcod,coords=c("X","Y"))
ggplot(pcod_s) + 
  geom_sf(aes(color = present)) +
  facet_wrap(~year) + 
  theme_minimal() +
  scale_color_gradient(low = "lightgrey", high = "red")

####### 空間メッシュの生成
n_site <- nrow(unique(pcod[,c("X","Y")]))
mesh   <- make_mesh(pcod, xy_cols = c("X", "Y"), n_knots =  min(n_site-1, 2000))
plot(mesh)

####### 時空間ロジスティック回帰（数分かかります）
mod    <- sdmTMB(present ~ s(depth) + as.factor( year ),
  data = pcod,mesh = mesh, family = binomial(link = "logit"),
  spatial = "on", time = "year", spatiotemporal = "ar1")

####### パラメータ推定
mod
tidy(mod, conf.int = TRUE) # スライド無し; 詳しくは後ほど
tidy(mod, effects = "ran_pars", conf.int = TRUE) # 同上
AIC(mod)

####### 水深の非線形効果
visreg(mod, xvar = "depth")

####### 予測
qcs_grid$year<-2017
p2   <- predict(mod, newdata = qcs_grid, type="response")
p2[1:3,]
p2_s <- st_as_sf(p2,coords=c("X","Y"))
ggplot(p2_s) + geom_sf(aes(color = est),pch=15,cex=0.5) +
  facet_wrap(~year) + theme_minimal() +
  scale_color_gradient(low = "lightgrey", high = "red")

####### 予測値の不確実性評価
p2_sim     <- predict(fit2, newdata = qcs_grid, type="response", nsim=500)
p2_s$est_sd<- apply(p2_sim, 1, sd)
ggplot(p2_s) + 
  geom_sf(aes(color = est_sd),pch=15,cex=0.5) +
  facet_wrap(~year) + 
  theme_minimal() +
  scale_color_gradient(low = "lightgrey", high = "red")

###############################################
##### 地理的制約を考慮する場合の時空間モデリング #####
###############################################

#install.packages(c("remotes","rnaturalearth"))
#remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE)
#remotes::install_github("ropensci/rnaturalearthhires")

library(rnaturalearthhires)
library(rnaturalearth)
library(sdmTMB)
library(sdmTMBextra)
library(sf)
library(ggplot2)

####### データの読み込みと地図化
data(pcod)                                     # Pacific codデータ
pcod      <-data.frame(pcod[pcod$year>2011, c("year","X","Y","depth","density","present","lat","lon")])
pcod$X    <- pcod$X*1000
pcod$Y    <- pcod$Y*1000
pcod_s    <- st_as_sf(pcod,coords=c("X","Y"), crs=3156)

####### 地形データ
canada_map<- ne_countries(scale = "large",country="canada")# カナダの陸地データ
pcod_map  <- st_crop(canada_map, c(xmin = -132, ymin = 50, xmax = -127, ymax = 53))#対象地域付近の陸地
pcod_map2 <- st_transform(pcod_map,crs=3156)

ggplot(pcod_map2) + geom_sf() + theme_minimal() + geom_sf(data = pcod_s, size = 1)

####### 地形を考慮した空間グラフの生成と地図化
n_site    <- nrow(unique(pcod[,c("X","Y")]))
mesh      <- make_mesh(pcod, xy_cols = c("X", "Y"), n_knots = min(n_site-1, 2000))
bmesh     <- add_barrier_mesh(spde_obj = mesh, barrier_sf = pcod_map2)
mesh_normal <- bmesh$mesh_sf[bmesh$normal_triangles , ]
mesh_barrier<- bmesh$mesh_sf[bmesh$barrier_triangles, ]
ggplot(pcod_map2) + geom_sf() +
  geom_sf(data = mesh_normal , size = 0.8, colour = "blue") +
  geom_sf(data = mesh_barrier, size = 2, colour = "red")

####### 空間ロジスティック回帰
fit       <- sdmTMB(present ~ s(depth),data = pcod, mesh = bmesh,
                    family = binomial(link = "logit"),spatial = "on")

####### 予測
data(qcs_grid)              # 2kmグリッドデータ
qcs_grid$X<- qcs_grid$X*1000# X座標をkmからm単位に変更
qcs_grid$Y<- qcs_grid$Y*1000# Y座標をkmからm単位に変更
p         <- predict(fit, newdata = qcs_grid, type="response")
p_s       <- st_as_sf(p,coords=c("X","Y"))
ggplot(p_s) + 
  geom_sf(aes(color = est),pch=15) +
  theme_minimal() +
  scale_color_gradient(low = "lightgrey", high = "red")

####### 予測値の不確実性評価
p2_sim     <- predict(fit, newdata = qcs_grid, type="response", nsim=500)
p2_s$est_sd<- apply(p2_sim, 1, sd)
ggplot(p2_s) + 
  geom_sf(aes(color = est_sd),pch=15,cex=0.5) +
  facet_wrap(~year) + 
  theme_minimal() +
  scale_color_gradient(low = "lightgrey", high = "red")
