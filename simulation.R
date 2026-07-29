library(irlba) 
library(movMF)
library(MASS)
library(Matrix)
library(easyCODA)
library(methods)
require(VGAM)
library(expm)
library(xtable)
library(ggplot2)
library(patchwork)

source('simulation_basic_function.R')

#### simulation run ####
#### PSSAR estimation and test ####
for (rep_ind_current in 1:200) {
  sample_set <- c(200,500,1000)
  rho_set <- c(-0.7,-0.3,0,0.1,0.4,0.9)
  neighbor_vec_set <- c(10,20,30)
  
  est_list <- list()
  test_list <- list()
  
  
  
  for (sample_current in sample_set) {
    list_ind <- which(sample_current==sample_set)
    
    est_results_matrix <- matrix(NA,nrow = 4,ncol = length(rho_set))
    test_results_matrix <- matrix(NA,nrow = 8,ncol = length(rho_set))
    
    est_results_r1 <- NULL
    est_results_r2 <- NULL
    est_results_r3 <- NULL
    est_results_r4 <- NULL
    test_results_r1 <- NULL
    test_results_r2 <- NULL
    test_results_r3 <- NULL
    test_results_r4 <- NULL
    test_results_r5 <- NULL
    test_results_r6 <- NULL
    test_results_r7 <- NULL
    test_results_r8 <- NULL
    
    
    for (rho in rho_set) {
      if (sample_current==200){
        ## row 1
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=6,neighbor_set=10, rho0=rho, dim_red='PCA')
        est_results_r1 <- c(est_results_r1,sim_results[1])
        test_results_r1 <- c(test_results_r1,sim_results[2])
        
        ## row 2
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=111,neighbor_set=10, rho0=rho, dim_red='PCA')
        est_results_r2 <- c(est_results_r2,sim_results[1])
        test_results_r2 <- c(test_results_r2,sim_results[2])
        
        ## row 5
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=6,neighbor_set=10, rho0=rho, dim_red='bootstrap')
        test_results_r5 <- c(test_results_r5,sim_results[2])
        
        ## row 6
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=111,neighbor_set=10, rho0=rho, dim_red='bootstrap')
        test_results_r6 <- c(test_results_r6,sim_results[2])
      } else{
        ## row 1
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=6,neighbor_set=10, rho0=rho, dim_red='PCA')
        est_results_r1 <- c(est_results_r1,sim_results[1])
        test_results_r1 <- c(test_results_r1,sim_results[2])
        
        ## row 2
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=111,neighbor_set=10, rho0=rho, dim_red='PCA')
        est_results_r2 <- c(est_results_r2,sim_results[1])
        test_results_r2 <- c(test_results_r2,sim_results[2])
        
        ## row 3
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=6,neighbor_set=neighbor_vec_set[list_ind], rho0=rho, dim_red='PCA')
        est_results_r3 <- c(est_results_r3,sim_results[1])
        test_results_r3 <- c(test_results_r3,sim_results[2])
        
        ## row 4
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=111,neighbor_set=neighbor_vec_set[list_ind], rho0=rho, dim_red='PCA')
        est_results_r4 <- c(est_results_r4,sim_results[1])
        test_results_r4 <- c(test_results_r4,sim_results[2])
        
        ## row 5
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=6,neighbor_set=10, rho0=rho, dim_red='bootstrap')
        test_results_r5 <- c(test_results_r5,sim_results[2])
        
        ## row 6
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=111,neighbor_set=10, rho0=rho, dim_red='bootstrap')
        test_results_r6 <- c(test_results_r6,sim_results[2])
        
        ## row 7
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=6,neighbor_set=neighbor_vec_set[list_ind], rho0=rho, dim_red='bootstrap')
        test_results_r7 <- c(test_results_r7,sim_results[2])
        
        ## row 8
        sim_results <- PSSAR_direct_sim_function(rep_ind_current,sample_size=sample_current,dim_set=111,neighbor_set=neighbor_vec_set[list_ind], rho0=rho, dim_red='bootstrap')
        test_results_r8 <- c(test_results_r8,sim_results[2])
      }
      
    }
    
    if (sample_current==200){
      est_results_matrix[1,]<- est_results_r1
      test_results_matrix[1,] <- test_results_r1
      est_results_matrix[2,]<- est_results_r2
      test_results_matrix[2,] <- test_results_r2
      test_results_matrix[5,] <- test_results_r5
      test_results_matrix[6,] <- test_results_r6
    } else{
      est_results_matrix[1,]<- est_results_r1
      test_results_matrix[1,] <- test_results_r1
      est_results_matrix[2,]<- est_results_r2
      test_results_matrix[2,] <- test_results_r2
      est_results_matrix[3,]<- est_results_r3
      test_results_matrix[3,] <- test_results_r3
      est_results_matrix[4,]<- est_results_r4
      test_results_matrix[4,] <- test_results_r4
      test_results_matrix[5,] <- test_results_r5
      test_results_matrix[6,] <- test_results_r6
      test_results_matrix[7,] <- test_results_r7
      test_results_matrix[8,] <- test_results_r8
    }
    
    
    
    est_list[[list_ind]] <- est_results_matrix
    test_list[[list_ind]] <- test_results_matrix
  }
  
  save(est_list,file=paste0('est_result_seed',rep_ind_current,'.RData'))
  save(test_list,paste0('test_result_seed',rep_ind_current,'.RData'))
}


#### estimation ####
files <- list.files(pattern = "^est_result.*\\.RData$", full.names = TRUE)


sample_set <- c(200,500,1000)
rho_set <- c(-0.7,-0.3,0,0.1,0.4,0.9)


est_bias_matrix_neighfix <- matrix(NA, nrow = 2*length(sample_set), ncol = length(rho_set))
est_mse_matrix_neighfix <- matrix(NA, nrow = 2*length(sample_set), ncol = length(rho_set))
est_bias_matrix_neighgrow <- matrix(NA, nrow = 2*length(sample_set), ncol = length(rho_set))
est_mse_matrix_neighgrow <- matrix(NA, nrow = 2*length(sample_set), ncol = length(rho_set))


for (sample_current in sample_set) {
  list_ind <- which(sample_current==sample_set)
  est_results_matrix_neighfix <- matrix(NA,nrow = 2*length(rho_set),ncol = length(files))
  est_results_matrix_neighgrow <- matrix(NA,nrow = 2*length(rho_set),ncol = length(files))
  for (rho in rho_set) {
    row_ind <- which(rho==rho_set)
    for (fileind in 1:length(files)) {
      
      result_current <- get(load(files[fileind]))
      
      fix_sample <- result_current[[list_ind]]
      
      est_results_matrix_neighfix[row_ind,fileind] <- fix_sample[1,row_ind]
      est_results_matrix_neighfix[(row_ind+length(rho_set)),fileind] <- fix_sample[2,row_ind]
      
      if (sample_current==200){
        est_results_matrix_neighgrow[row_ind,fileind] <- fix_sample[1,row_ind]
        est_results_matrix_neighgrow[(row_ind+length(rho_set)),fileind] <- fix_sample[2,row_ind]
      } else{
        est_results_matrix_neighgrow[row_ind,fileind] <- fix_sample[3,row_ind]
        est_results_matrix_neighgrow[(row_ind+length(rho_set)),fileind] <- fix_sample[4,row_ind]
      }
      
    }
  }
  
  
  for (rho in rho_set) {
    row_ind <- which(rho==rho_set)
    
    # --- neighfix ---
    
    # 1st parameter set Bias and RMSE
    est_bias_matrix_neighfix[list_ind,row_ind] <- mean(est_results_matrix_neighfix[row_ind,] - rho)
    est_mse_matrix_neighfix[list_ind,row_ind] <- sqrt(mean((est_results_matrix_neighfix[row_ind,] - rho)^2))
    
    # 2nd parameter set Bias and RMSE
    est_bias_matrix_neighfix[(list_ind+length(sample_set)),row_ind] <- mean(est_results_matrix_neighfix[(row_ind+length(rho_set)),] - rho)
    est_mse_matrix_neighfix[(list_ind+length(sample_set)),row_ind] <- sqrt(mean((est_results_matrix_neighfix[(row_ind+length(rho_set)),] - rho)^2))
    
    
    # --- neighgrow ---
    
    # 1st parameter set Bias and RMSE
    est_bias_matrix_neighgrow[list_ind,row_ind] <- mean(est_results_matrix_neighgrow[row_ind,] - rho)
    est_mse_matrix_neighgrow[list_ind,row_ind] <- sqrt(mean((est_results_matrix_neighgrow[row_ind,] - rho)^2))
    
    # 2nd parameter set Bias and RMSE
    est_bias_matrix_neighgrow[(list_ind+length(sample_set)),row_ind] <- mean(est_results_matrix_neighgrow[(row_ind+length(rho_set)),] - rho)
    est_mse_matrix_neighgrow[(list_ind+length(sample_set)),row_ind] <- sqrt(mean((est_results_matrix_neighgrow[(row_ind+length(rho_set)),] - rho)^2))
  }
  
}


## fix neighbor
table_df <- as.data.frame(matrix(NA,nrow = 2*length(rho_set), ncol = (2*length(sample_set)+1)))
table_df[,1] <- rep(rho_set,2)
table_df[1:6,c(2,4,6)] <- t(est_bias_matrix_neighfix[1:3,])
table_df[7:12,c(2,4,6)] <- t(est_bias_matrix_neighfix[4:6,])
table_df[1:6,c(3,5,7)] <- t(est_mse_matrix_neighfix[1:3,])
table_df[7:12,c(3,5,7)] <- t(est_mse_matrix_neighfix[4:6,])

# Print LaTeX table
formatted_df <- data.frame(lapply(table_df, function(x) sprintf("%.4f", x)))
print(xtable(formatted_df, align = c("l", rep("c", 7))), 
      include.rownames = FALSE, 
      sanitize.text.function = identity)

## growing neighbor
table_df <- as.data.frame(matrix(NA,nrow = 2*length(rho_set), ncol = (2*length(sample_set)+1)))
table_df[,1] <- rep(rho_set,2)
table_df[1:6,c(2,4,6)] <- t(est_bias_matrix_neighgrow[1:3,])
table_df[7:12,c(2,4,6)] <- t(est_bias_matrix_neighgrow[4:6,])
table_df[1:6,c(3,5,7)] <- t(est_mse_matrix_neighgrow[1:3,])
table_df[7:12,c(3,5,7)] <- t(est_mse_matrix_neighgrow[4:6,])

# Print LaTeX table
formatted_df <- data.frame(lapply(table_df, function(x) sprintf("%.4f", x)))
print(xtable(formatted_df, align = c("l", rep("c", 7))), 
      include.rownames = FALSE, 
      sanitize.text.function = identity)

#### test ####
files <- list.files(pattern = "^test_result.*\\.RData$", full.names = TRUE)


sample_set <- c(200,500,1000)
rho_set <- c(-0.7,-0.3,0,0.1,0.4,0.9)


test_matrix_neighfix <- matrix(NA, nrow = 4*length(sample_set), ncol = length(rho_set))
test_matrix_neighgrow <- matrix(NA, nrow = 4*length(sample_set), ncol = length(rho_set))

threshold_set <- qchisq(0.05, df=1, lower.tail=FALSE)

for (sample_current in sample_set) {
  list_ind <- which(sample_current==sample_set)
  test_results_matrix_neighfix <- matrix(NA,nrow = 4*length(rho_set),ncol = length(files))
  test_results_matrix_neighgrow <- matrix(NA,nrow = 4*length(rho_set),ncol = length(files))
  for (rho in rho_set) {
    row_ind <- which(rho==rho_set)
    for (fileind in 1:length(files)) {
      
      result_current <- get(load(files[fileind]))
      
      fix_sample <- result_current[[list_ind]]
      
      test_results_matrix_neighfix[row_ind,fileind] <- fix_sample[1,row_ind]
      test_results_matrix_neighfix[(row_ind+length(rho_set)),fileind] <- fix_sample[5,row_ind]
      test_results_matrix_neighfix[(row_ind+2*length(rho_set)),fileind] <- fix_sample[2,row_ind]
      test_results_matrix_neighfix[(row_ind+3*length(rho_set)),fileind] <- fix_sample[6,row_ind]
      
      if (sample_current==200){
        test_results_matrix_neighgrow[row_ind,fileind] <- fix_sample[1,row_ind]
        test_results_matrix_neighgrow[(row_ind+length(rho_set)),fileind] <- fix_sample[5,row_ind]
        test_results_matrix_neighgrow[(row_ind+2*length(rho_set)),fileind] <- fix_sample[2,row_ind]
        test_results_matrix_neighgrow[(row_ind+3*length(rho_set)),fileind] <- fix_sample[6,row_ind]
      } else{
        test_results_matrix_neighgrow[row_ind,fileind] <- fix_sample[3,row_ind]
        test_results_matrix_neighgrow[(row_ind+length(rho_set)),fileind] <- fix_sample[7,row_ind]
        test_results_matrix_neighgrow[(row_ind+2*length(rho_set)),fileind] <- fix_sample[4,row_ind]
        test_results_matrix_neighgrow[(row_ind+3*length(rho_set)),fileind] <- fix_sample[8,row_ind]
      }
      
    }
  }
  for (rho in rho_set) {
    row_ind <- which(rho==rho_set)
    test_matrix_neighfix[list_ind,row_ind] <- mean(test_results_matrix_neighfix[row_ind,])
    test_matrix_neighfix[(list_ind+length(sample_set)),row_ind] <- mean(test_results_matrix_neighfix[(row_ind+length(rho_set)),])
    test_matrix_neighfix[(list_ind+2*length(sample_set)),row_ind] <- mean(test_results_matrix_neighfix[(row_ind+2*length(rho_set)),])
    test_matrix_neighfix[(list_ind+3*length(sample_set)),row_ind] <- mean(test_results_matrix_neighfix[(row_ind+3*length(rho_set)),])
    test_matrix_neighgrow[list_ind,row_ind] <- mean(test_results_matrix_neighgrow[row_ind,])
    test_matrix_neighgrow[(list_ind+length(sample_set)),row_ind] <- mean(test_results_matrix_neighgrow[(row_ind+length(rho_set)),])
    test_matrix_neighgrow[(list_ind+2*length(sample_set)),row_ind] <- mean(test_results_matrix_neighgrow[(row_ind+2*length(rho_set)),])
    test_matrix_neighgrow[(list_ind+3*length(sample_set)),row_ind] <- mean(test_results_matrix_neighgrow[(row_ind+3*length(rho_set)),])
  }
}

test_matrix_neighfix
test_matrix_neighgrow


## dim 6 fix neigh
title_list <- list('(a) n = 200','(b) n = 500', '(c) n = 1000')
plot_list <- list()

for (sample_current in sample_set) {
  list_ind <- which(sample_current==sample_set)
  
  df_plot <- data.frame(matrix(0,nrow = 2*length(rho_set),ncol = 3))
  df_plot[,1] <- rep(c(0,1,2,3,4,5),2)
  df_plot[,2] <- c(rep('PCA',length(rho_set)),rep('bootstrap',length(rho_set)))
  df_plot[1:length(rho_set),3] <- test_matrix_neighfix[list_ind,c(3,4,2,5,1,6)]
  df_plot[(length(rho_set)+1):(2*length(rho_set)),3] <- test_matrix_neighfix[list_ind+length(sample_set),c(3,4,2,5,1,6)]
  
  colnames(df_plot) <- c('Signal','Method','RejectResults')
  df_plot$RejectResults <- 100 * df_plot$RejectResults
  if (list_ind==1){ 
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  } else{
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  }
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]]  + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

## dim 111 fix neigh
title_list <- list('(a) n = 200','(b) n = 500', '(c) n = 1000')
plot_list <- list()

for (sample_current in sample_set) {
  list_ind <- which(sample_current==sample_set)
  
  df_plot <- data.frame(matrix(0,nrow = 2*length(rho_set),ncol = 3))
  df_plot[,1] <- rep(c(0,1,2,3,4,5),2)
  df_plot[,2] <- c(rep('PCA',length(rho_set)),rep('bootstrap',length(rho_set)))
  df_plot[1:length(rho_set),3] <- test_matrix_neighfix[list_ind+2*length(sample_set),c(3,4,2,5,1,6)]
  df_plot[(length(rho_set)+1):(2*length(rho_set)),3] <- test_matrix_neighfix[list_ind+3*length(sample_set),c(3,4,2,5,1,6)]
  
  colnames(df_plot) <- c('Signal','Method','RejectResults')
  df_plot$RejectResults <- 100 * df_plot$RejectResults
  if (list_ind==1){ 
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  } else{
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  }
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]]  + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")


## dim 6 grow neigh
title_list <- list('(a) n = 200','(b) n = 500', '(c) n = 1000')
plot_list <- list()

for (sample_current in sample_set) {
  list_ind <- which(sample_current==sample_set)
  
  df_plot <- data.frame(matrix(0,nrow = 2*length(rho_set),ncol = 3))
  df_plot[,1] <- rep(c(0,1,2,3,4,5),2)
  df_plot[,2] <- c(rep('PCA',length(rho_set)),rep('bootstrap',length(rho_set)))
  df_plot[1:length(rho_set),3] <- test_matrix_neighgrow[list_ind,c(3,4,2,5,1,6)]
  df_plot[(length(rho_set)+1):(2*length(rho_set)),3] <- test_matrix_neighgrow[list_ind+length(sample_set),c(3,4,2,5,1,6)]
  
  colnames(df_plot) <- c('Signal','Method','RejectResults')
  df_plot$RejectResults <- 100 * df_plot$RejectResults
  if (list_ind==1){ 
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  } else{
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  }
  
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]]  + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

## dim 111 grow neigh
title_list <- list('(a) n = 200','(b) n = 500', '(c) n = 1000')
plot_list <- list()

for (sample_current in sample_set) {
  list_ind <- which(sample_current==sample_set)
  
  df_plot <- data.frame(matrix(0,nrow = 2*length(rho_set),ncol = 3))
  df_plot[,1] <- rep(c(0,1,2,3,4,5),2)
  df_plot[,2] <- c(rep('PCA',length(rho_set)),rep('bootstrap',length(rho_set)))
  df_plot[1:length(rho_set),3] <- test_matrix_neighgrow[list_ind+2*length(sample_set),c(3,4,2,5,1,6)]
  df_plot[(length(rho_set)+1):(2*length(rho_set)),3] <- test_matrix_neighgrow[list_ind+3*length(sample_set),c(3,4,2,5,1,6)]
  
  colnames(df_plot) <- c('Signal','Method','RejectResults')
  df_plot$RejectResults <- 100 * df_plot$RejectResults
  if (list_ind==1){ 
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  } else{
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  }
  
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]]  + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")


#### PSSAR prediction ####

for (rep_ind_current in 1:200) {
  sample_set <- c(200,400,600,800,1000)
  neigh_vec <- c(10,20,30,40,50)
  
  ## sim_results_PSSAR_matrix
  ## row 1 --- d=6, neigh=10, ##% (n=200,400,600,800,1000)
  ## row 2 --- d=6, neigh=grow, ##% (n=200,400,600,800,1000)
  ## row 3 --- d=111, neigh=10, ##% (n=200,400,600,800,1000)
  ## row 4 --- d=111, neigh=grow, ##% (n=200,400,600,800,1000)
  
  
  ## coverage_matrix  & PI_width_matrix
  ## row 1 --- d=6, neigh=10, 95% (n=200,400,600,800,1000)
  ## row 2 --- d=6, neigh=10, 90% (n=200,400,600,800,1000)
  ## row 3 --- d=6, neigh=10, 80% (n=200,400,600,800,1000)
  ## row 4 --- d=111, neigh=10, 95% (n=200,400,600,800,1000)
  ## row 5 --- d=111, neigh=10, 90% (n=200,400,600,800,1000)
  ## row 6 --- d=111, neigh=10, 80% (n=200,400,600,800,1000)
  ## row 7 --- d=6, neigh=grow, 95% (n=200,400,600,800,1000)
  ## row 8 --- d=6, neigh=grow, 90% (n=200,400,600,800,1000)
  ## row 9 --- d=6, neigh=grow, 80% (n=200,400,600,800,1000)
  ## row 10 --- d=111, neigh=grow, 95% (n=200,400,600,800,1000)
  ## row 11 --- d=111, neigh=grow, 90% (n=200,400,600,800,1000)
  ## row 12 --- d=111, neigh=grow, 80% (n=200,400,600,800,1000)
  
  
  sim_results_PSSAR_matrix <- matrix(NA,nrow = 4,ncol = length(sample_set))
  coverage_matrix <- matrix(NA,nrow = 12,ncol = length(sample_set))
  PI_width_matrix <- matrix(NA,nrow = 12,ncol = length(sample_set))
  
  
  for (sample_current in sample_set) {
    col_ind <- which(sample_current==sample_set)
    
    pred_error <- PSSAR_sim_prediction_int_function(rep_ind_current,sample_size=sample_current,d=6,neigh_set = 10, alpha = c(0.05,0.1,0.15))
    sim_results_PSSAR_matrix[1,col_ind] <- pred_error[1]
    coverage_matrix[c(1,2,3),col_ind] <- pred_error[c(2,4,6)]
    PI_width_matrix[c(1,2,3),col_ind] <- pred_error[c(3,5,7)]
    
    if (sample_current==200){
      sim_results_PSSAR_matrix[2,col_ind] <- pred_error[1]
      coverage_matrix[c(7,8,9),col_ind] <- pred_error[c(2,4,6)]
      PI_width_matrix[c(7,8,9),col_ind] <- pred_error[c(3,5,7)]
    } else{
      pred_error <- PSSAR_sim_prediction_int_function(rep_ind_current,sample_size=sample_current,d=6,neigh_set = neigh_vec[col_ind], alpha = c(0.05,0.1,0.15))
      sim_results_PSSAR_matrix[2,col_ind] <- pred_error[1]
      coverage_matrix[c(7,8,9),col_ind] <- pred_error[c(2,4,6)]
      PI_width_matrix[c(7,8,9),col_ind] <- pred_error[c(3,5,7)]
    }
    
    
    pred_error <- PSSAR_sim_prediction_int_function(rep_ind_current,sample_size=sample_current,d=111,neigh_set = 10, alpha = c(0.05,0.1,0.15))
    sim_results_PSSAR_matrix[3,col_ind] <- pred_error[1]
    coverage_matrix[c(4,5,6),col_ind] <- pred_error[c(2,4,6)]
    PI_width_matrix[c(4,5,6),col_ind] <- pred_error[c(3,5,7)]
    
    if (sample_current==200){
      sim_results_PSSAR_matrix[4,col_ind] <- pred_error[1]
      coverage_matrix[c(10,11,12),col_ind] <- pred_error[c(2,4,6)]
      PI_width_matrix[c(10,11,12),col_ind] <- pred_error[c(3,5,7)]
    } else{
      pred_error <- PSSAR_sim_prediction_int_function(rep_ind_current,sample_size=sample_current,d=111,neigh_set = neigh_vec[col_ind], alpha = c(0.05,0.1,0.15))
      sim_results_PSSAR_matrix[4,col_ind] <- pred_error[1]
      coverage_matrix[c(10,11,12),col_ind] <- pred_error[c(2,4,6)]
      PI_width_matrix[c(10,11,12),col_ind] <- pred_error[c(3,5,7)]
    }
    
  }
  
  result_list <- list(sim_results_PSSAR_matrix,coverage_matrix,PI_width_matrix)
  
  save(result_list,paste0('pred_result_seed',rep_ind_current,'.RData'))
}


files <- list.files(pattern = "^pred_result.*\\.RData$", full.names = TRUE)


sample_set <-  c(200,400,600,800,1000)
result_matrix_neighfix <- matrix(NA,nrow = 14,ncol = 5)
result_matrix_neighgrow <- matrix(NA,nrow = 14,ncol = 5)


for (sample_current in sample_set) {
  col_ind <- which(sample_current==sample_set)
  
  pred_results_d6_neighfix <- rep(NA,length(files))
  coverage_results_d6_neighfix <- matrix(NA,nrow = 3,ncol = length(files))
  width_results_d6_neighfix <- matrix(NA,nrow = 3,ncol = length(files))
  
  pred_results_d111_neighfix <- rep(NA,length(files))
  coverage_results_d111_neighfix <- matrix(NA,nrow = 3,ncol = length(files))
  width_results_d111_neighfix <- matrix(NA,nrow = 3,ncol = length(files))
  
  pred_results_d6_neighgrow <- rep(NA,length(files))
  coverage_results_d6_neighgrow <- matrix(NA,nrow = 3,ncol = length(files))
  width_results_d6_neighgrow <- matrix(NA,nrow = 3,ncol = length(files))
  
  pred_results_d111_neighgrow <- rep(NA,length(files))
  coverage_results_d111_neighgrow <- matrix(NA,nrow = 3,ncol = length(files))
  width_results_d111_neighgrow <- matrix(NA,nrow = 3,ncol = length(files))
  
  for (fileind in 1:length(files)) {
    
    result_current <- get(load(files[fileind]))
    
    SSAR_sample <- result_current[[1]]
    coverage_sample <- result_current[[2]]
    width_sample <- result_current[[3]]
    
    pred_results_d6_neighfix[fileind] <- SSAR_sample[1,col_ind]
    coverage_results_d6_neighfix[,fileind] <- coverage_sample[1:3,col_ind]
    width_results_d6_neighfix[,fileind] <- width_sample[1:3,col_ind]
    
    pred_results_d111_neighfix[fileind] <- SSAR_sample[3,col_ind]
    coverage_results_d111_neighfix[,fileind] <- coverage_sample[4:6,col_ind]
    width_results_d111_neighfix[,fileind] <- width_sample[4:6,col_ind]
    
    pred_results_d6_neighgrow[fileind] <- SSAR_sample[2,col_ind]
    coverage_results_d6_neighgrow[,fileind] <- coverage_sample[7:9,col_ind]
    width_results_d6_neighgrow[,fileind] <- width_sample[7:9,col_ind]
    
    pred_results_d111_neighgrow[fileind] <- SSAR_sample[4,col_ind]
    coverage_results_d111_neighgrow[,fileind] <- coverage_sample[10:12,col_ind]
    width_results_d111_neighgrow[,fileind] <- width_sample[10:12,col_ind]
    
  }
  coverage_mean_d6_neighfix <- rowMeans(coverage_results_d6_neighfix)
  coverage_mean_d111_neighfix <- rowMeans(coverage_results_d111_neighfix)
  width_mean_d6_neighfix <- rowMeans(width_results_d6_neighfix)
  width_mean_d111_neighfix <- rowMeans(width_results_d111_neighfix)
  
  coverage_mean_d6_neighgrow <- rowMeans(coverage_results_d6_neighgrow)
  coverage_mean_d111_neighgrow <- rowMeans(coverage_results_d111_neighgrow)
  width_mean_d6_neighgrow <- rowMeans(width_results_d6_neighgrow)
  width_mean_d111_neighgrow <- rowMeans(width_results_d111_neighgrow)
  
  
  result_matrix_neighfix[,col_ind] <- c(
    mean(pred_results_d6_neighfix),coverage_mean_d6_neighfix[1],width_mean_d6_neighfix[1],
    coverage_mean_d6_neighfix[2],width_mean_d6_neighfix[2],
    coverage_mean_d6_neighfix[3],width_mean_d6_neighfix[3],
    mean(pred_results_d111_neighfix),coverage_mean_d111_neighfix[1],width_mean_d111_neighfix[1],
    coverage_mean_d111_neighfix[2],width_mean_d111_neighfix[2],
    coverage_mean_d111_neighfix[3],width_mean_d111_neighfix[3]
  )
  
  result_matrix_neighgrow[,col_ind] <- c(
    mean(pred_results_d6_neighgrow),coverage_mean_d6_neighgrow[1],width_mean_d6_neighgrow[1],
    coverage_mean_d6_neighgrow[2],width_mean_d6_neighgrow[2],
    coverage_mean_d6_neighgrow[3],width_mean_d6_neighgrow[3],
    mean(pred_results_d111_neighgrow),coverage_mean_d111_neighgrow[1],width_mean_d111_neighgrow[1],
    coverage_mean_d111_neighgrow[2],width_mean_d111_neighgrow[2],
    coverage_mean_d111_neighgrow[3],width_mean_d111_neighgrow[3]
  )
  
}


# Print LaTeX table
formatted_df <- round(result_matrix_neighfix, 4)
print(
  xtable(formatted_df,
         align = c("l", rep("c", 5)),
         digits = c(0, rep(4, 5))  # 4 digits for each column
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity
)


formatted_df <- round(result_matrix_neighgrow, 4)
print(
  xtable(formatted_df,
         align = c("l", rep("c", 5)),
         digits = c(0, rep(4, 5))  # 4 digits for each column
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity
)


## compare with MSAR

for (rep_ind_current in 1:200) {
  sample_set <- c(200,400,600,800,1000)

  sim_results_PSSAR_matrix <- rep(0,length(sample_set))
  sim_results_MSAR_matrix <- rep(0,length(sample_set))
  
  
  for (sample_current in sample_set) {
    col_ind <- which(sample_current==sample_set)
    
    pred_error <- PSSAR_trans_sim_prediction_function2(rep_ind_current,sample_size=sample_current)
    sim_results_PSSAR_matrix[col_ind] <- pred_error[1]
    sim_results_MSAR_matrix[col_ind] <- pred_error[2]
    
  }
  
  result_list <- list(sim_results_PSSAR_matrix,sim_results_MSAR_matrix)
  
  save(result_list,paste0('pred_result_seed',rep_ind_current,'.RData'))
}


files <- list.files(pattern = "^pred_result.*\\.RData$", full.names = TRUE)



sample_set <-  c(200,400,600,800,1000)
title_list <- list('(a) n=200','(b) n=400','(c) n=600','(d) n=800','(e) n=1000')
plot_list <- list()

for (sample_current in sample_set) {
  col_ind <- which(sample_current==sample_set)
  pred_results_matrix <- matrix(NA,nrow = 2,ncol = length(files))
  for (fileind in 1:length(files)) {
    
    result_current <- get(load(files[fileind]))
    
    
    SSAR_sample <- result_current[[1]]
    MSAR_sample <- result_current[[2]]
    
    pred_results_matrix[1,fileind] <- SSAR_sample[col_ind]
    pred_results_matrix[2,fileind] <- MSAR_sample[col_ind]
    
  }
  # Row 1 = SSAR errors, Row 2 = MSAR errors
  n <- ncol(pred_results_matrix)
  
  df <- data.frame(
    Model = rep(c("PSSAR", "MSAR"), each = n),
    Error = c(pred_results_matrix[1, ], pred_results_matrix[2, ])
  )
  
  print(round(rowMeans(pred_results_matrix),4))
  
  # Draw violin plot
  
  if(col_ind==1){
    plot_list[[col_ind]] <- ggplot(df, aes(x = Model, y = Error, fill = Model)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4) +
      scale_fill_manual(values = c(
        "PSSAR" = "#D55E00",  # orange
        "MSAR" = "#56B4E9" # blue
      )) +
      scale_y_continuous(limits = c(0, 2.5)) +   # force y-axis to start at 0
      labs(title = title_list[[col_ind]],
           y = "Prediction Error") +
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
    plot_list[[col_ind]] <- ggplot(df, aes(x = Model, y = Error, fill = Model)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4) +
      scale_fill_manual(values = c(
        "PSSAR" = "#D55E00",  # orange
        "MSAR" = "#56B4E9" # blue
      )) +
      scale_y_continuous(limits = c(0, 2.5)) +   # force y-axis to start at 0
      labs(title = title_list[[col_ind]],
           y = "Prediction Error") +
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
  
  
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]]  + plot_layout(ncol = 5, guides = "collect") &
  theme(legend.position = "bottom")

#### SRMSAR test ####
for (rep_ind_current in 1:200) {
  sample_set_default <- c(200,500,1000)
  sample_set <- c(200,500,1000)
  rho_set <- c(-0.7,-0.3,0,0.1,0.4,0.9)
  neighbor_vec_set <- c(10,20,30)
  
  test_list <- list()

  for (sample_current in sample_set) {
    
    log_msg(paste("Starting sample_current =", sample_current))
    
    list_ind <- which(sample_current == sample_set_default)
    
    est_results_matrix <- matrix(NA, nrow = 4, ncol = length(rho_set))
    test_results_matrix <- matrix(NA, nrow = 8, ncol = length(rho_set))
    
    est_results_r1 <- NULL; test_results_r1 <- NULL
    est_results_r2 <- NULL; test_results_r2 <- NULL
    est_results_r3 <- NULL; test_results_r3 <- NULL
    est_results_r4 <- NULL; test_results_r4 <- NULL
    test_results_r5 <- NULL; test_results_r6 <- NULL
    test_results_r7 <- NULL; test_results_r8 <- NULL
    
    for (rho in rho_set) {
      
      if (sample_current == 200) {
        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 6,
          neighbor_set = 10, rho0 = rho, dim_red = 'PCA'
        )
        est_results_r1 <- c(est_results_r1, sim_results[1])
        test_results_r1 <- c(test_results_r1, sim_results[2])

        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 111,
          neighbor_set = 10, rho0 = rho, dim_red = 'PCA'
        )
        est_results_r2 <- c(est_results_r2, sim_results[1])
        test_results_r2 <- c(test_results_r2, sim_results[2])

        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 6,
          neighbor_set = 10, rho0 = rho, dim_red = 'bootstrap'
        )
        test_results_r5 <- c(test_results_r5, sim_results[2])
        log_msg("Finished row 5")
        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 111,
          neighbor_set = 10, rho0 = rho, dim_red = 'bootstrap'
        )
        test_results_r6 <- c(test_results_r6, sim_results[2])

        
      } else {
        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 6,
          neighbor_set = 10, rho0 = rho, dim_red = 'PCA'
        )
        est_results_r1 <- c(est_results_r1, sim_results[1])
        test_results_r1 <- c(test_results_r1, sim_results[2])

        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 111,
          neighbor_set = 10, rho0 = rho, dim_red = 'PCA'
        )
        est_results_r2 <- c(est_results_r2, sim_results[1])
        test_results_r2 <- c(test_results_r2, sim_results[2])

        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 6,
          neighbor_set = neighbor_vec_set[list_ind], rho0 = rho, dim_red = 'PCA'
        )
        est_results_r3 <- c(est_results_r3, sim_results[1])
        test_results_r3 <- c(test_results_r3, sim_results[2])
        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 111,
          neighbor_set = neighbor_vec_set[list_ind], rho0 = rho, dim_red = 'PCA'
        )
        est_results_r4 <- c(est_results_r4, sim_results[1])
        test_results_r4 <- c(test_results_r4, sim_results[2])
        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 6,
          neighbor_set = 10, rho0 = rho, dim_red = 'bootstrap'
        )
        test_results_r5 <- c(test_results_r5, sim_results[2])
        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 111,
          neighbor_set = 10, rho0 = rho, dim_red = 'bootstrap'
        )
        test_results_r6 <- c(test_results_r6, sim_results[2])
        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 6,
          neighbor_set = neighbor_vec_set[list_ind], rho0 = rho, dim_red = 'bootstrap'
        )
        test_results_r7 <- c(test_results_r7, sim_results[2])
        
        sim_results <- SSARR_tangent_space_sim_function(
          rep_ind_current, sample_size = sample_current, dim_set = 111,
          neighbor_set = neighbor_vec_set[list_ind], rho0 = rho, dim_red = 'bootstrap'
        )
        test_results_r8 <- c(test_results_r8, sim_results[2])
      }
    }
    
    
    # Store results
    if (sample_current == 200) {
      est_results_matrix[1,] <- est_results_r1
      test_results_matrix[1,] <- test_results_r1
      est_results_matrix[2,] <- est_results_r2
      test_results_matrix[2,] <- test_results_r2
      test_results_matrix[5,] <- test_results_r5
      test_results_matrix[6,] <- test_results_r6
    } else {
      est_results_matrix[1,] <- est_results_r1
      test_results_matrix[1,] <- test_results_r1
      est_results_matrix[2,] <- est_results_r2
      test_results_matrix[2,] <- test_results_r2
      est_results_matrix[3,] <- est_results_r3
      test_results_matrix[3,] <- test_results_r3
      est_results_matrix[4,] <- est_results_r4
      test_results_matrix[4,] <- test_results_r4
      test_results_matrix[5,] <- test_results_r5
      test_results_matrix[6,] <- test_results_r6
      test_results_matrix[7,] <- test_results_r7
      test_results_matrix[8,] <- test_results_r8
    }
    
    
    est_list[[list_ind]] <- est_results_matrix
    test_list[[list_ind]] <- test_results_matrix
    
  }
  save(test_list,paste0('test_result_seed',rep_ind_current,'.RData'))
}



files <- list.files(pattern = "^test_result.*\\.RData$", full.names = TRUE)


sample_set_default <- c(200,500,1000)
sample_set <- c(200,500,1000)
rho_set <- c(-0.7,-0.3,0,0.1,0.4,0.9)


test_matrix_neighfix <- matrix(NA, nrow = 4*length(sample_set_default), ncol = length(rho_set))
test_matrix_neighgrow <- matrix(NA, nrow = 4*length(sample_set_default), ncol = length(rho_set))

threshold_set <- qchisq(0.05, df=1, lower.tail=FALSE)

for (sample_current in sample_set) {
  list_ind <- which(sample_current==sample_set_default)
  test_results_matrix_neighfix <- matrix(NA,nrow = 4*length(rho_set),ncol = length(files))
  test_results_matrix_neighgrow <- matrix(NA,nrow = 4*length(rho_set),ncol = length(files))
  for (rho in rho_set) {
    row_ind <- which(rho==rho_set)
    for (fileind in 1:length(files)) {
      
      result_current <- get(load(files[fileind]))
      
      fix_sample <- result_current[[list_ind]]
      
      test_results_matrix_neighfix[row_ind,fileind] <- fix_sample[1,row_ind]
      test_results_matrix_neighfix[(row_ind+length(rho_set)),fileind] <- fix_sample[5,row_ind]
      test_results_matrix_neighfix[(row_ind+2*length(rho_set)),fileind] <- fix_sample[2,row_ind]
      test_results_matrix_neighfix[(row_ind+3*length(rho_set)),fileind] <- fix_sample[6,row_ind]
      
      if (sample_current==200){
        test_results_matrix_neighgrow[row_ind,fileind] <- fix_sample[1,row_ind]
        test_results_matrix_neighgrow[(row_ind+length(rho_set)),fileind] <- fix_sample[5,row_ind]
        test_results_matrix_neighgrow[(row_ind+2*length(rho_set)),fileind] <- fix_sample[2,row_ind]
        test_results_matrix_neighgrow[(row_ind+3*length(rho_set)),fileind] <- fix_sample[6,row_ind]
      } else{
        test_results_matrix_neighgrow[row_ind,fileind] <- fix_sample[3,row_ind]
        test_results_matrix_neighgrow[(row_ind+length(rho_set)),fileind] <- fix_sample[7,row_ind]
        test_results_matrix_neighgrow[(row_ind+2*length(rho_set)),fileind] <- fix_sample[4,row_ind]
        test_results_matrix_neighgrow[(row_ind+3*length(rho_set)),fileind] <- fix_sample[8,row_ind]
      }
      
    }
  }
  for (rho in rho_set) {
    row_ind <- which(rho==rho_set)
    test_matrix_neighfix[list_ind,row_ind] <- mean(test_results_matrix_neighfix[row_ind,])
    test_matrix_neighfix[(list_ind+length(sample_set_default)),row_ind] <- mean(test_results_matrix_neighfix[(row_ind+length(rho_set)),])
    test_matrix_neighfix[(list_ind+2*length(sample_set_default)),row_ind] <- mean(test_results_matrix_neighfix[(row_ind+2*length(rho_set)),])
    test_matrix_neighfix[(list_ind+3*length(sample_set_default)),row_ind] <- mean(test_results_matrix_neighfix[(row_ind+3*length(rho_set)),])
    test_matrix_neighgrow[list_ind,row_ind] <- mean(test_results_matrix_neighgrow[row_ind,])
    test_matrix_neighgrow[(list_ind+length(sample_set_default)),row_ind] <- mean(test_results_matrix_neighgrow[(row_ind+length(rho_set)),])
    test_matrix_neighgrow[(list_ind+2*length(sample_set_default)),row_ind] <- mean(test_results_matrix_neighgrow[(row_ind+2*length(rho_set)),])
    test_matrix_neighgrow[(list_ind+3*length(sample_set_default)),row_ind] <- mean(test_results_matrix_neighgrow[(row_ind+3*length(rho_set)),])
  }
}



test_matrix_neighfix
test_matrix_neighgrow


## dim 6 fix neigh
title_list <- list('(a) n = 200','(b) n = 500', '(c) n = 1000')
plot_list <- list()

for (sample_current in sample_set_default) {
  list_ind <- which(sample_current==sample_set_default)
  
  df_plot <- data.frame(matrix(0,nrow = 2*length(rho_set),ncol = 3))
  df_plot[,1] <- rep(c(0,1,2,3,4,5),2)
  df_plot[,2] <- c(rep('PCA',length(rho_set)),rep('bootstrap',length(rho_set)))
  df_plot[1:length(rho_set),3] <- test_matrix_neighfix[list_ind,c(3,4,2,5,1,6)]
  df_plot[(length(rho_set)+1):(2*length(rho_set)),3] <- test_matrix_neighfix[list_ind+length(sample_set_default),c(3,4,2,5,1,6)]
  
  colnames(df_plot) <- c('Signal','Method','RejectResults')
  df_plot$RejectResults <- 100 * df_plot$RejectResults
  
  if (list_ind==1){ 
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  } else{
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  }
  
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]]  + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

## dim 111 fix neigh
title_list <- list('(a) n = 200','(b) n = 500', '(c) n = 1000')
plot_list <- list()

for (sample_current in sample_set_default) {
  list_ind <- which(sample_current==sample_set_default)
  
  df_plot <- data.frame(matrix(0,nrow = 2*length(rho_set),ncol = 3))
  df_plot[,1] <- rep(c(0,1,2,3,4,5),2)
  df_plot[,2] <- c(rep('PCA',length(rho_set)),rep('bootstrap',length(rho_set)))
  df_plot[1:length(rho_set),3] <- test_matrix_neighfix[list_ind+2*length(sample_set_default),c(3,4,2,5,1,6)]
  df_plot[(length(rho_set)+1):(2*length(rho_set)),3] <- test_matrix_neighfix[list_ind+3*length(sample_set_default),c(3,4,2,5,1,6)]
  
  colnames(df_plot) <- c('Signal','Method','RejectResults')
  df_plot$RejectResults <- 100 * df_plot$RejectResults
  
  if (list_ind==1){
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  } else{
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  }
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]]  + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")


## dim 6 grow neigh
title_list <- list('(a) n = 200','(b) n = 500', '(c) n = 1000')
plot_list <- list()

for (sample_current in sample_set_default) {
  list_ind <- which(sample_current==sample_set_default)
  
  df_plot <- data.frame(matrix(0,nrow = 2*length(rho_set),ncol = 3))
  df_plot[,1] <- rep(c(0,1,2,3,4,5),2)
  df_plot[,2] <- c(rep('PCA',length(rho_set)),rep('bootstrap',length(rho_set)))
  df_plot[1:length(rho_set),3] <- test_matrix_neighgrow[list_ind,c(3,4,2,5,1,6)]
  df_plot[(length(rho_set)+1):(2*length(rho_set)),3] <- test_matrix_neighgrow[list_ind+length(sample_set_default),c(3,4,2,5,1,6)]
  
  colnames(df_plot) <- c('Signal','Method','RejectResults')
  df_plot$RejectResults <- 100 * df_plot$RejectResults
  
  if (list_ind==1){ 
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
  } else{
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (/%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  }
  
  
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]]  + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

## dim 111 grow neigh
title_list <- list('(a) n = 200','(b) n = 500', '(c) n = 1000')
plot_list <- list()

for (sample_current in sample_set_default) {
  list_ind <- which(sample_current==sample_set_default)
  
  df_plot <- data.frame(matrix(0,nrow = 2*length(rho_set),ncol = 3))
  df_plot[,1] <- rep(c(0,1,2,3,4,5),2)
  df_plot[,2] <- c(rep('PCA',length(rho_set)),rep('bootstrap',length(rho_set)))
  df_plot[1:length(rho_set),3] <- test_matrix_neighgrow[list_ind+2*length(sample_set_default),c(3,4,2,5,1,6)]
  df_plot[(length(rho_set)+1):(2*length(rho_set)),3] <- test_matrix_neighgrow[list_ind+3*length(sample_set_default),c(3,4,2,5,1,6)]
  
  colnames(df_plot) <- c('Signal','Method','RejectResults')
  df_plot$RejectResults <- 100 * df_plot$RejectResults
  
  if (list_ind==1){
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  } else{
    plot_list[[list_ind]] <- ggplot(df_plot,aes(x = Signal, y = RejectResults, color = Method, shape = Method, linetype = Method)) + 
      geom_point() +
      geom_line() +
      geom_hline(yintercept = 5, linetype="dotted") + 
      scale_y_continuous(limits = c(0, 100)) +
      labs(x='Signal strength',y='Percentage of Rejections (%)',title = title_list[[list_ind]])+
      theme_classic()+
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            legend.text=element_text(size=14),
            plot.margin=margin(5,15,5,5),
            panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=10),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))
    
  }
  
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]]  + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")


#### SRMSAR prediction ####


for (rep_ind_current in 1:200) {
  sample_set <- c(200,400,600,800,1000)
  neigh_vec <- c(10,20,30,40,50)
  
  sample_set <- c(200,400,600,800,1000)
  neigh_vec <- c(10,20,30,40,50)

  
  sim_results_PSSAR_matrix <- matrix(NA,nrow = 4,ncol = length(sample_set))
  coverage_matrix <- matrix(NA,nrow = 12,ncol = length(sample_set))
  PI_width_matrix <- matrix(NA,nrow = 12,ncol = length(sample_set))
  
  
  for (sample_current in sample_set) {
    col_ind <- which(sample_current==sample_set)
    
    pred_error <- SSARR_sim_prediction_int_function(rep_ind_current,sample_size=sample_current,d=6,neigh_set = 10, alpha = c(0.05,0.1,0.2))
    sim_results_PSSAR_matrix[1,col_ind] <- pred_error[1]
    coverage_matrix[c(1,2,3),col_ind] <- pred_error[c(2,4,6)]
    PI_width_matrix[c(1,2,3),col_ind] <- pred_error[c(3,5,7)]
    
    if (sample_current==200){
      sim_results_PSSAR_matrix[2,col_ind] <- pred_error[1]
      coverage_matrix[c(7,8,9),col_ind] <- pred_error[c(2,4,6)]
      PI_width_matrix[c(7,8,9),col_ind] <- pred_error[c(3,5,7)]
    } else{
      pred_error <- SSARR_sim_prediction_int_function(rep_ind_current,sample_size=sample_current,d=6,neigh_set = neigh_vec[col_ind], alpha = c(0.05,0.1,0.15))
      sim_results_PSSAR_matrix[2,col_ind] <- pred_error[1]
      coverage_matrix[c(7,8,9),col_ind] <- pred_error[c(2,4,6)]
      PI_width_matrix[c(7,8,9),col_ind] <- pred_error[c(3,5,7)]
    }
    
    
    pred_error <- SSARR_sim_prediction_int_function(rep_ind_current,sample_size=sample_current,d=111,neigh_set = 10, alpha = c(0.05,0.1,0.15))
    sim_results_PSSAR_matrix[3,col_ind] <- pred_error[1]
    coverage_matrix[c(4,5,6),col_ind] <- pred_error[c(2,4,6)]
    PI_width_matrix[c(4,5,6),col_ind] <- pred_error[c(3,5,7)]
    
    if (sample_current==200){
      sim_results_PSSAR_matrix[4,col_ind] <- pred_error[1]
      coverage_matrix[c(10,11,12),col_ind] <- pred_error[c(2,4,6)]
      PI_width_matrix[c(10,11,12),col_ind] <- pred_error[c(3,5,7)]
    } else{
      pred_error <- SSARR_sim_prediction_int_function(rep_ind_current,sample_size=sample_current,d=111,neigh_set = neigh_vec[col_ind], alpha = c(0.05,0.1,0.15))
      sim_results_PSSAR_matrix[4,col_ind] <- pred_error[1]
      coverage_matrix[c(10,11,12),col_ind] <- pred_error[c(2,4,6)]
      PI_width_matrix[c(10,11,12),col_ind] <- pred_error[c(3,5,7)]
    }
    
  }
  
  result_list <- list(sim_results_PSSAR_matrix,coverage_matrix,PI_width_matrix)
  
  
  save(result_list,paste0('pred_result_seed',rep_ind_current,'.RData'))
}

files <- list.files(pattern = "^pred_result.*\\.RData$", full.names = TRUE)


sample_set <-  c(200,400,600,800,1000)
result_matrix_neighfix <- matrix(NA,nrow = 14,ncol = 5)
result_matrix_neighgrow <- matrix(NA,nrow = 14,ncol = 5)





for (sample_current in sample_set) {
  col_ind <- which(sample_current==sample_set)
  
  pred_results_d6_neighfix <- rep(NA,length(files))
  coverage_results_d6_neighfix <- matrix(NA,nrow = 3,ncol = length(files))
  width_results_d6_neighfix <- matrix(NA,nrow = 3,ncol = length(files))
  
  pred_results_d111_neighfix <- rep(NA,length(files))
  coverage_results_d111_neighfix <- matrix(NA,nrow = 3,ncol = length(files))
  width_results_d111_neighfix <- matrix(NA,nrow = 3,ncol = length(files))
  
  pred_results_d6_neighgrow <- rep(NA,length(files))
  coverage_results_d6_neighgrow <- matrix(NA,nrow = 3,ncol = length(files))
  width_results_d6_neighgrow <- matrix(NA,nrow = 3,ncol = length(files))
  
  pred_results_d111_neighgrow <- rep(NA,length(files))
  coverage_results_d111_neighgrow <- matrix(NA,nrow = 3,ncol = length(files))
  width_results_d111_neighgrow <- matrix(NA,nrow = 3,ncol = length(files))
  
  for (fileind in 1:length(files)) {
    
    result_current <- get(load(files[fileind]))
    
    SSAR_sample <- result_current[[1]]
    coverage_sample <- result_current[[2]]
    width_sample <- result_current[[3]]
    
    pred_results_d6_neighfix[fileind] <- SSAR_sample[1,col_ind]
    coverage_results_d6_neighfix[,fileind] <- coverage_sample[1:3,col_ind]
    width_results_d6_neighfix[,fileind] <- width_sample[1:3,col_ind]
    
    pred_results_d111_neighfix[fileind] <- SSAR_sample[3,col_ind]
    coverage_results_d111_neighfix[,fileind] <- coverage_sample[4:6,col_ind]
    width_results_d111_neighfix[,fileind] <- width_sample[4:6,col_ind]
    
    pred_results_d6_neighgrow[fileind] <- SSAR_sample[2,col_ind]
    coverage_results_d6_neighgrow[,fileind] <- coverage_sample[7:9,col_ind]
    width_results_d6_neighgrow[,fileind] <- width_sample[7:9,col_ind]
    
    pred_results_d111_neighgrow[fileind] <- SSAR_sample[4,col_ind]
    coverage_results_d111_neighgrow[,fileind] <- coverage_sample[10:12,col_ind]
    width_results_d111_neighgrow[,fileind] <- width_sample[10:12,col_ind]
    
  }
  coverage_mean_d6_neighfix <- rowMeans(coverage_results_d6_neighfix)
  coverage_mean_d111_neighfix <- rowMeans(coverage_results_d111_neighfix)
  width_mean_d6_neighfix <- rowMeans(width_results_d6_neighfix)
  width_mean_d111_neighfix <- rowMeans(width_results_d111_neighfix)
  
  coverage_mean_d6_neighgrow <- rowMeans(coverage_results_d6_neighgrow)
  coverage_mean_d111_neighgrow <- rowMeans(coverage_results_d111_neighgrow)
  width_mean_d6_neighgrow <- rowMeans(width_results_d6_neighgrow)
  width_mean_d111_neighgrow <- rowMeans(width_results_d111_neighgrow)
  
  
  result_matrix_neighfix[,col_ind] <- c(
    mean(pred_results_d6_neighfix),coverage_mean_d6_neighfix[1],width_mean_d6_neighfix[1],
    coverage_mean_d6_neighfix[2],width_mean_d6_neighfix[2],
    coverage_mean_d6_neighfix[3],width_mean_d6_neighfix[3],
    mean(pred_results_d111_neighfix),coverage_mean_d111_neighfix[1],width_mean_d111_neighfix[1],
    coverage_mean_d111_neighfix[2],width_mean_d111_neighfix[2],
    coverage_mean_d111_neighfix[3],width_mean_d111_neighfix[3]
  )
  
  result_matrix_neighgrow[,col_ind] <- c(
    mean(pred_results_d6_neighgrow),coverage_mean_d6_neighgrow[1],width_mean_d6_neighgrow[1],
    coverage_mean_d6_neighgrow[2],width_mean_d6_neighgrow[2],
    coverage_mean_d6_neighgrow[3],width_mean_d6_neighgrow[3],
    mean(pred_results_d111_neighgrow),coverage_mean_d111_neighgrow[1],width_mean_d111_neighgrow[1],
    coverage_mean_d111_neighgrow[2],width_mean_d111_neighgrow[2],
    coverage_mean_d111_neighgrow[3],width_mean_d111_neighgrow[3]
  )
  
}


# Print LaTeX table
formatted_df <- round(result_matrix_neighfix, 4)
print(
  xtable(formatted_df,
         align = c("l", rep("c", 5)),
         digits = c(0, rep(4, 5))  # 4 digits for each column
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity
)


formatted_df <- round(result_matrix_neighgrow, 4)
print(
  xtable(formatted_df,
         align = c("l", rep("c", 5)),
         digits = c(0, rep(4, 5))  # 4 digits for each column
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity
)

## compare with MSAR

for (rep_ind_current in 1:200) {
  sample_set <- c(200,400,600,800,1000)
  
  sim_results_PSSAR_matrix <- rep(0,length(sample_set))
  sim_results_MSAR_matrix <- rep(0,length(sample_set))
  
  
  for (sample_current in sample_set) {
    col_ind <- which(sample_current==sample_set)
    
    pred_error <- SSARR_trans_sim_prediction_function2(rep_ind_current,sample_size=sample_current)
    sim_results_PSSAR_matrix[col_ind] <- pred_error[1]
    sim_results_MSAR_matrix[col_ind] <- pred_error[2]
    
  }
  
  result_list <- list(sim_results_PSSAR_matrix,sim_results_MSAR_matrix)
  
  save(result_list,paste0('pred_result_seed',rep_ind_current,'.RData'))
}


files <- list.files(pattern = "^pred_result.*\\.RData$", full.names = TRUE)


sample_set <-  c(200,400,600,800,1000)
title_list <- list('(a) n=200','(b) n=400','(c) n=600','(d) n=800','(e) n=1000')
plot_list <- list()

for (sample_current in sample_set) {
  col_ind <- which(sample_current==sample_set)
  pred_results_matrix <- matrix(NA,nrow = 2,ncol = length(files))
  for (fileind in 1:length(files)) {
    
    result_current <- get(load(files[fileind]))
    
    
    SSAR_sample <- result_current[[1]]
    MSAR_sample <- result_current[[2]]
    
    pred_results_matrix[1,fileind] <- SSAR_sample[col_ind]
    pred_results_matrix[2,fileind] <- MSAR_sample[col_ind]
    
  }
  # Row 1 = SSAR errors, Row 2 = MSAR errors
  n <- ncol(pred_results_matrix)
  
  df <- data.frame(
    Model = rep(c("SRMSAR", "MSAR"), each = n),
    Error = c(pred_results_matrix[1, ], pred_results_matrix[2, ])
  )
  
  print(round(rowMeans(pred_results_matrix),4))
  
  # Draw violin plot
  
  if(col_ind==1){
    plot_list[[col_ind]] <- ggplot(df, aes(x = Model, y = Error, fill = Model)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4) +
      scale_fill_manual(values = c(
        "SRMSAR" = "darkseagreen",  # orange
        "MSAR" = "#56B4E9" # blue
      )) +
      scale_y_continuous(limits = c(0, 3.5)) +   # force y-axis to start at 0
      labs(title = title_list[[col_ind]],
           y = "Prediction Error") +
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
    plot_list[[col_ind]] <- ggplot(df, aes(x = Model, y = Error, fill = Model)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4) +
      scale_fill_manual(values = c(
        "SRMSAR" = "darkseagreen",  # orange
        "MSAR" = "#56B4E9" # blue
      )) +
      scale_y_continuous(limits = c(0, 3.5)) +   # force y-axis to start at 0
      labs(title = title_list[[col_ind]],
           y = "Prediction Error") +
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
  
  
  
}

plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]]  + plot_layout(ncol = 5, guides = "collect") &
  theme(legend.position = "bottom")
