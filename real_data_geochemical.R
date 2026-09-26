library(robCompositions)
library(sf)
library(ggplot2)
library(patchwork)
library(tidyr)
library(easyCODA)


source('real_basic_function.R')

###########################################################################
###########################################################################
###########################################################################
library(MASS)

# Please load MSAR code from https://github.com/XueningZhu/MSAR_code

###########################################################################
###########################################################################
###########################################################################


#### GEMAS DATA ####
data(gemas)
str(gemas)

SPA_data<-gemas[gemas$COUNTRY=='SPA',]

indx<-c(12,14,16,17,23)
SPA_data_now<-cbind(SPA_data[,2:3],SPA_data[,12:29])


## five major elements
for (i in 1:dim(SPA_data)[1]) {
  for (j in indx) {
    SPA_data_now[i,(j-9)]<-SPA_data[i,j]/sum(SPA_data[i,12:29])
  }
}

SPA_data_now<- subset(SPA_data_now, select = c(1,2,indx-9) )
SPA_data_now<-cbind(SPA_data_now,rep(0,dim(SPA_data)[1]))

for (i in 1:dim(SPA_data)[1]) {
  SPA_data_now[i,8]<-1-sum(SPA_data_now[i,3:7])
}

colnames(SPA_data_now)[8]<-'remaining'



SPA_data_now
SPA_data_now<-SPA_data_now[-12,]
dim(SPA_data_now)

location_sf <- st_as_sf(x = SPA_data_now[,1:2], 
                        coords = c("longitude", "latitude"),
                        ## crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
                        crs = 4326)


dist_result<-st_distance(location_sf)

first_location<-as.vector(as.matrix(dist_result)[1:dim(SPA_data_now)[1],1])
sort_value<-sort(first_location[first_location>0], decreasing=F)
first_neigh_ind<-which(first_location<=sort_value[8])


neigh_ind_matrix_save<-matrix(NA,202,9)

for (i in 1:202) {
  first_location<-as.vector(as.matrix(dist_result)[,i])
  first_neigh_ind<-which(first_location<=80000)
  ##sort_value<-sort(first_location[first_location>0], decreasing=F)
  ##first_neigh_ind<-which(first_location<=sort_value[3])
  
  num<-9
  if (length(first_neigh_ind)>9){
    num<-9
    print(i)
  }
  
  for (j in 1:min(num,length(first_neigh_ind))) {
    neigh_ind_matrix_save[i,j]<-first_neigh_ind[j]
  }
  
}


###########################################################################
###########################################################################

## 6 major
y <- as.matrix(sqrt(SPA_data_now[,3:8]))
## all
##y <- as.matrix(sqrt(SPA_data_now[,3:20]))

matplot(t(y), type = "l", lty = 1, col = 1:nrow(y),
        xlab = "Chemical", ylab = "Percentage",  xaxt = "n")

axis(
  side = 1,
  at = 1:6,
  labels = c("Al", "Ca", "Fe", "K", "Si",
             "Remaining")
)

matplot(t(y), type = "p", pch = 19, col = 1:nrow(y),
        xlab = "Composition", ylab = "Value",
        main = "Geochemical Data")

# initialize matrix
n <- nrow(y)

## first order neighbor
# initialize matrix with dimnames
W <- matrix(0, nrow = n, ncol = n)

# fill in adjacency weights
for (i in 1:n) {
  neigh <- setdiff(neigh_ind_matrix_save[i,],i)
  # only assign if neighbors exist
  if (length(neigh) > 0) {
    W[i, neigh] <- 1 / length(neigh)
  }
}

#############################################

#### estimation ####
system.time({ PSSAR_fit_function_fast(y,W,dim_red='No test') })
##user  system elapsed 
##0.138   0.002   0.140 

system.time({ clr_MSAR_fit_function(y,W) })


fit_result <- PSSAR_fit_function_fast(y,W,dim_red='bootstrap')

rho_hat <- fit_result[1]
rho_hat 
rej_result <- fit_result[2]
rej_result 




## prediction
## consider use 1,...,n-1 to predict n
SSAR_predict_angle <- NULL
SAR_predict_angle <- NULL
SSAR_predict_add <- NULL
SAR_predict_add <- NULL
CLR_predict_angle <- NULL
CLR_predict_add <- NULL


for (i in 1:n) {
  predict_combine <- predict_function(y,W,i,measure_set='MSE')
  SSAR_predict_angle <- c(SSAR_predict_angle,predict_combine[[1]])
  SAR_predict_angle <- c(SAR_predict_angle,predict_combine[[2]])
  SSAR_predict_add <- c(SSAR_predict_add,predict_combine[[3]])
  SAR_predict_add <- c(SAR_predict_add,predict_combine[[4]])
  
  ##clr_predict <- clr_predict_function(y,W,i,measure_set='MSE')
  clr_predict <- clr_MSAR_predict_function2(y,W,i,measure_set='MSE')
  CLR_predict_angle <- c(CLR_predict_angle,clr_predict[[1]])
  CLR_predict_add <- c(CLR_predict_add,clr_predict[[2]])

  print(paste('SSAR:',predict_combine[[1]],'MSAR:',clr_predict[[1]]))
}



df1 <- data.frame(
  index = rep(1:length(SSAR_predict_angle), 2),
  angle = c(SSAR_predict_angle, CLR_predict_angle),
  group = rep(c("SSAR", "MSAR"), each = length(SSAR_predict_angle))
)

df2 <- data.frame(
  index = rep(1:length(SSAR_predict_angle), 4),
  mse = c(SSAR_predict_add, CLR_predict_add),
  group = rep(c("SSAR", "MSAR"), each = length(SSAR_predict_angle))
)

result_list <- list(df1,df2)
##save(result_list,file = 'geochemical_prediction.RData')
result_list <- get(load('geochemical_prediction.RData'))
df1 <- result_list[[1]]
df1[df1$group=='SSAR',]$group <- 'PSSAR'
df2 <- result_list[[2]]
df2[df2$group=='SSAR',]$group <- 'PSSAR'

# prediction error

# Combine into a data frame



p1 <- ggplot(df1, aes(x = index, y = angle, color = group, shape = group)) +
  geom_point(size = 2, fill = NA, stroke = 1) +
  scale_color_manual(values = c(
    "PSSAR" = "#D55E00",  # orange
    "MSAR" = "#56B4E9", # blue
    "Extrinsic SAR 2" = "#009E73",  # green
    "Local Averaging" = "black"
  )) +
  scale_shape_manual(values = c(
    "PSSAR" = 24,    # triangle
    "MSAR" = 22, # square
    "Extrinsic SAR 2" = 18,  # diamond
    "Local Averaging" = 4
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
    "MSAR" = "#56B4E9", # blue
    "Extrinsic SAR 2" = "#009E73",  # green
    "Local Averaging" = "black"
  )) +
  scale_shape_manual(values = c(
    "PSSAR" = 24,    # triangle
    "MSAR" = 22, # square
    "Extrinsic SAR 2" = 18,  # diamond
    "Local Averaging" = 4
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
title_list <- c('Location 1',
                'Location 2',
                'Location 3',
                'Location 4',
                'Location 5',
                'Location 6')

for (i in c(4, 5, 7, 13, 27,46)) {
  predict_combine <- predict_function(y, W, i, measure_set = 'MSE')
  clr_predict <- clr_MSAR_predict_function2(y,W,i,measure_set='MSE')
  ##krg_predict <- krig_predict_function(y,W,i,measure_set='MSE')
  
  df <- data.frame(matrix(NA, nrow = 6, ncol = 4))
  colnames(df) <- c('Chemical',
                    'True Observation',
                    'PSSAR Prediction',
                    'MSAR Prediction')
                    ##'Extrinsic SAR Prediction 2',
                    ##'Local Averaging')
  

  df$Chemical <- factor(
    c('Al', 'Ca', 'Fe', 'K', 'Si', 'Remaining'),
    levels = c('Al', 'Ca', 'Fe', 'K', 'Si', 'Remaining')
  )
  df[, 2] <- predict_combine[[5]]
  df[, 3] <- predict_combine[[6]]
  ##df[, 4] <- predict_combine[[7]]
  df[, 4] <- clr_predict[[3]] 
  ##df[, 6] <- krg_predict[[3]] 
  
  df_long <- df %>%
    pivot_longer(
      cols = c(
        "True Observation",
        "PSSAR Prediction",
        "MSAR Prediction"
      ),
      names_to = "Model",
      values_to = "Value"
    )
  
  # Plot
  if (ind%in%c(1,4)){
    plot_list[[ind]] <- ggplot(df_long,
                               aes(
                                 x = Chemical,
                                 y = Value,
                                 color = Model,
                                 shape = Model,
                                 group = Model
                               )) +
      geom_point(size = 3, fill = NA, stroke = 1.4) +
      scale_color_manual(values = c(
        "True Observation" = "pink",
        "PSSAR Prediction" = "#D55E00",  # orange
        "MSAR Prediction" = "#56B4E9" # blue
        ##"Extrinsic SAR Prediction 2" = "#009E73",  # green
        ##"Local Averaging" = "black" # blue
      )) +
      scale_shape_manual(values = c(
        "True Observation" = 21,  # hollow circle
        "PSSAR Prediction" = 24,  # hollow triangle
        "MSAR Prediction" = 22    # hollow square
      )) +
      labs(
        title = title_list[[ind]],
        x = "Chemical",
        y = "Percentage",
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
        axis.title = element_text(face = "bold")
      )
  } else{
    plot_list[[ind]] <- ggplot(df_long,
                               aes(
                                 x = Chemical,
                                 y = Value,
                                 color = Model,
                                 shape = Model,
                                 group = Model
                               )) +
      geom_point(size = 3, fill = NA, stroke = 1.4) +
      scale_color_manual(values = c(
        "True Observation" = "pink",
        "PSSAR Prediction" = "#D55E00",  # orange
        "MSAR Prediction" = "#56B4E9" # blue
        ##"Extrinsic SAR Prediction 2" = "#009E73",  # green
        ##"Local Averaging" = "black" # blue
      )) +
      scale_shape_manual(values = c(
        "True Observation" = 21,  # hollow circle
        "PSSAR Prediction" = 24,  # hollow triangle
        "MSAR Prediction" = 22    # hollow square
      )) +
      labs(
        title = title_list[[ind]],
        x = "Chemical",
        y = "Percentage",
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
        axis.title = element_text(face = "bold")
      )
  }
  
  
  ind <- ind + 1
  print(ind)
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


#### SSAREXR ####

SPA_data<-SPA_data[-12,]


## covariate
covariate <- SPA_data[,c('soilclass','MeanTemp','AnnPrec')]
covariate[,1] <- as.numeric(as.factor(covariate[,1]))
X <- as.matrix(covariate)

dim(X)

fit_result <- SRMSAR_fit_function(y,W,X,dim_red='bootstrap')
rho_hat <- fit_result[1]
rho_hat 
rej_result <- fit_result[2]
rej_result 


## prediction
## consider use 1,...,n-1 to predict n
SSAR_predict_angle <- NULL
SSAREXR_predict_angle <- NULL
SSAR_predict_add <- NULL
SSAREXR_predict_add <- NULL


for (i in 1:n) {
  predict_combine <- predict_function(y,W,i,measure_set='MSE')
  SSAR_predict_angle <- c(SSAR_predict_angle,predict_combine[[1]])
  SSAR_predict_add <- c(SSAR_predict_add,predict_combine[[3]])
  
  
  SSAREXR_predict <- new_SSAREXR_predict_function(y,X,W,i,measure_set='MSE')
  
  
  SSAREXR_predict_angle <- c(SSAREXR_predict_angle,SSAREXR_predict[[1]])
  SSAREXR_predict_add <- c(SSAREXR_predict_add,SSAREXR_predict[[2]])
  
  print(paste('index:',i,'SSAR:',round(predict_combine[[1]],4),'SSAREXR:',round(SSAREXR_predict[[1]],4)))
}

mean(SSAR_predict_angle) ##0.2414158
mean(SSAREXR_predict_angle) 
## 0.2334643 using intrinsic mean, 0.2338368 using extrinsic mean
## 0.2352451 using new function



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



