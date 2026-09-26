library(sf)
library(spdep)
require(HMDHFDplus)
library(irlba)
library(tidyr)
library(ggplot2)
library(patchwork)
library(easyCODA)
library(readxl)
library(wwntests)

source('real_basic_function.R')




################
# Japanese data
################

state = c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", "Niigata",
          "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane",
          "Okayama", "Hiroshima", "Yamaguchi", "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")


population_vec <- rep(NA,length(state))

year_chose <- 2023

for(iwk in 1:length(state))
{
  index <- ifelse(iwk < 10, paste0("0", iwk), as.character(iwk))
  
  ## population
  population_table = readJMDweb(prefID = index, item = "Population", fixup = FALSE) 
  population_vec[iwk] <- round(sum(population_table[population_table$Year==year_chose,]$Total))
}


# Step 1: Download Japan prefecture boundaries (GADM level 1)
# This will download a GeoPackage (.gpkg) directly
##url <- "https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_JPN.gpkg"
##download.file(url, destfile = "japan.gpkg", mode = "wb")

# Step 2: Read in as sf object
japan_sf <- st_read("japan.gpkg", layer = "ADM_ADM_1")  # level 1 = prefectures

# Step 3: Inspect prefecture names
japan_sf$NAME_1

# Step 4: Build neighbors list
nb <- poly2nb(japan_sf, queen = TRUE)

# Step 5: Convert into a named list
pref_names <- japan_sf$NAME_1
neighbors_list <- setNames(lapply(nb, function(x) pref_names[x]), pref_names)


# initialize matrix
n <- length(state)
# initialize matrix with dimnames
W <- matrix(0, nrow = n, ncol = n, dimnames = list(state, state))

# fill in adjacency weights
for (i in state) {
  neigh <- neighbors_list[[i]]
  # only assign if neighbors exist
  if (length(neigh) > 0) {
    # make sure neighbors are valid city names
    neigh <- neigh[neigh %in% state]
    
    ## no population
    ##W[i, neigh] <- 1 / length(neigh)
    
    ## include population
    ind_j <- match(neigh,state)
    W[i, neigh] <- sqrt(population_vec[ind_j]) / sum(sqrt(population_vec[ind_j]))
  }
}

# check a row
##W["Aichi", ]

# Representative point of each prefecture
pref_centroid <- st_centroid(japan_sf)
dist_matrix <- st_distance(pref_centroid)


#### female ####
y <- matrix(NA,nrow = length(state),ncol = 111)
n <- length(state)


year_chose <- 2023

for(iwk in 1:length(state))
{
  index <- ifelse(iwk < 10, paste0("0", iwk), as.character(iwk))
  ## life table
  table = readJMDweb(prefID = index, item = "fltper_1x1", fixup = FALSE) ## female: fltper_1x1 male: mltper_1x1
  qx = table[table$Year==year_chose,]$qx
  
  # set radix
  start_pop = 10^5
  qx_accurate <- rep(0,length(qx)) 
  for(ik in 1:length(qx))
  {
    qx_accurate[ik] = qx[ik] * start_pop
    start_pop = start_pop - qx_accurate[ik]
  }
  dx <- 1 / length(qx)
  qx_normalize <- qx_accurate / (sum(qx_accurate) * dx)  ## normalize to density
  y[iwk,] <- qx_normalize/sqrt(sum(qx_normalize^2)) 
  ##y[iwk,] <- sqrt(qx_normalize) ## square root transformation 
}

##check,  square sum to 1
sum(y[1,]^2)

matplot(t(y^2) * 10^5, type = "l", lty = 1, col = 1:nrow(y),
        xlab = "Age", ylab = "Value",
        main = "47 Prefectures")

plot(0:110,(y[4,])^2 * 10^5, type = "l", lty = 1,
        xlab = "Age", ylab = "Value",
        main = "47 Prefectures")



#### estimation ####
fit_result <- PSSAR_fit_function_fast(y,W,dim_red='bootstrap')


rho_hat <- fit_result[1]
rho_hat ## male: 0.3161659, female: 0.3886082
## use population 0.3268027 for male
rej_result <- fit_result[2]
rej_result ## 1


## income
income_vec <- c(
  2742, 2507, 2841, 2945, 2697, 2897, 2943, 3327, 3479, 3283, 3047, 3116,
  5415, 3268, 2916, 3398, 3023, 3280, 3160, 3010, 2919, 3432, 3728, 3121,
  3318, 2983, 3190, 2968, 2632, 2913, 2515, 2667, 2769, 3109, 3199, 3092,
  3013, 2658, 2644, 2885, 2753, 2629, 2667, 2714, 2468, 2509, 2391
)

X <- matrix(NA,nrow = 47,ncol = 1)
X[,1] <- scale(income_vec)

fit_result <- SRMSAR_fit_function(y,W,X,dim_red='bootstrap')
rho_hat <- fit_result[1]
rho_hat ## 0.3939589
rej_result <- fit_result[2]
rej_result ## 1



## prediction
## consider use 1,...,n-1 to predict n
SSAR_predict_angle <- NULL
SSAREXR_predict_angle <- NULL
SSAR_predict_add <- NULL
SSAREXR_predict_add <- NULL


for (i in 1:n) {
  predict_combine <- predict_function(y,W,i,measure_set='JSD')
  SSAR_predict_angle <- c(SSAR_predict_angle,predict_combine[[1]])
  SSAR_predict_add <- c(SSAR_predict_add,predict_combine[[3]])
  
  
  ##SSAREXR_predict <- SSAREXR_predict_function(y,covariate,W,i,measure_set='MSE')
  ##SSAREXR_predict <- intr_SSAREXR_predict_function(y,covariate,W,i,measure_set='MSE')
  SSAREXR_predict <- new_SSAREXR_predict_function(y,X,W,i,measure_set='JSD')
  
  
  SSAREXR_predict_angle <- c(SSAREXR_predict_angle,SSAREXR_predict[[1]])
  SSAREXR_predict_add <- c(SSAREXR_predict_add,SSAREXR_predict[[2]])
  
  print(paste('index:',i,'SSAR:',round(predict_combine[[1]],4),'SSAREXR:',round(SSAREXR_predict[[1]],4)))
}

mean(SSAR_predict_angle) ##0.0340903
mean(SSAREXR_predict_angle) ## 0.03568196


# prediction error

# Combine into a data frame
df1 <- data.frame(
  index = rep(1:length(SSAR_predict_angle), 2),
  angle = c(SSAR_predict_angle, SSAREXR_predict_angle),
  group = rep(c("PSSAR", "SRMSAR"), each = length(SSAR_predict_angle))
)

df2 <- data.frame(
  index = rep(1:length(SSAR_predict_angle), 2),
  mse = c(SSAR_predict_add, SSAREXR_predict_add),
  group = rep(c("PSSAR", "SRMSAR"), each = length(SSAR_predict_angle))
)




p1 <- ggplot(df1, aes(x = index, y = angle, color = group, shape = group)) +
  geom_point(size = 2, fill = NA, stroke = 1) +
  scale_color_manual(values = c(
    "PSSAR" = "#D55E00",  # orange
    "SRMSAR" = "darkseagreen"
  )) +
  scale_shape_manual(values = c(
    "PSSAR" = 24,    # triangle
    "SRMSAR" = 23
  )) +
  labs(
    x = "Location Index",
    y = "Angle"
  ) +
  theme_classic()+
  theme(plot.title = element_text(face = "bold", hjust = 0.5,size=12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text=element_text(size=12),
        panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"))

p2 <- ggplot(df2, aes(x = index, y = mse, color = group, shape = group)) +
  geom_point(size = 2, fill = NA, stroke = 1) +
  scale_color_manual(values = c(
    "PSSAR" = "#D55E00",  # orange
    "SRMSAR" = "darkseagreen"
  )) +
  scale_shape_manual(values = c(
    "PSSAR" = 24,    # triangle
    "SRMSAR" = 23
  )) +
  labs(
    x = "Location Index",
    y = "MSE"
  ) +
  theme_classic()+
  theme(plot.title = element_text(face = "bold", hjust = 0.5,size=12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text=element_text(size=12),
        panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"))


p1 + p2 + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")


# prediction visualization
## i=1,2,202
ind <- 1
plot_list <- list()


for (i in c(1,16,17,22,41,47)) {
  predict_combine <- predict_function(y, W, i, measure_set = 'JSD')
  SSAREXR_predict <- new_SSAREXR_predict_function(y,X,W,i,measure_set='JSD')

  
  df <- data.frame(matrix(NA, nrow = 111, ncol = 4))
  colnames(df) <- c('Age',
                    'True Observation',
                    'PSSAR Prediction',
                    'SRMSAR Prediction')
  df$Age <- c(0:110)
  df[, 2] <- predict_combine[[5]] * 10^5
  df[, 3] <- predict_combine[[6]] * 10^5
  df[, 4] <- SSAREXR_predict[[4]] * 10^5

  
  df_long <- df %>%
    pivot_longer(
      cols = c(
        'True Observation',
        'PSSAR Prediction',
        'SRMSAR Prediction'
      ),
      names_to = "Model",
      values_to = "Value"
    )
  
  # Plot
  if (ind%in%c(1,4)){
    plot_list[[ind]] <- ggplot(df_long,
                               aes(
                                 x = Age,
                                 y = Value,
                                 color = Model,
                                 linetype = Model,
                                 group = Model
                               )) +
      geom_line(size = 1.2) +
      scale_linetype_manual(values = c(
        "True Observation" = "dotted",
        "PSSAR Prediction" = "dashed",
        "SRMSAR Prediction" = "dotdash"
      )) +
      scale_color_manual(values = c(
        "True Observation" = "pink",
        "PSSAR Prediction" = "#D55E00",  # orange
        "SRMSAR Prediction" = "darkseagreen" # blue
        ##"Extrinsic SAR Prediction 2" = "#009E73",  # green
        ##"Local Averaging" = "black" # blue
      )) +
      labs(
        title = state[i],
        x = "Age",
        y = "Life-Table Death Counts",
        color = "Model",
        shape = "Model"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 12
        ),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.key.width = unit(2, "cm")
      )
  } else{
    plot_list[[ind]] <- ggplot(df_long,
                               aes(
                                 x = Age,
                                 y = Value,
                                 color = Model,
                                 linetype = Model,
                                 group = Model
                               )) +
      geom_line(size = 1.2) +
      scale_linetype_manual(values = c(
        "True Observation" = "dotted",
        "PSSAR Prediction" = "dashed",
        "SRMSAR Prediction" = "dotdash"
      )) +
      scale_color_manual(values = c(
        "True Observation" = "pink",
        "PSSAR Prediction" = "#D55E00",  # orange
        "SRMSAR Prediction" = "darkseagreen" # blue
        ##"Extrinsic SAR Prediction 2" = "#009E73",  # green
        ##"Local Averaging" = "black" # blue
      )) +
      labs(
        title = state[i],
        x = "Age",
        y = "Life-Table Death Counts",
        color = "Model",
        shape = "Model"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 12
        ),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.key.width = unit(2, "cm")
      )
  }
  
  
  
  ind <- ind + 1
}


plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]] + plot_list[[6]] + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")


## prediction interval
PI_vec <- NULL
for (i in 1:n) {
  pred_int_results <- PI_function(y,W,i,dist_result,0.1)
  print(pred_int_results)
  PI_vec <- c(PI_vec,pred_int_results)
}

mean(PI_vec)



