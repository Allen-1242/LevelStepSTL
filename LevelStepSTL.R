#LevelStepSTL

#Novel Method to decompose a levelShifted time series 

#library import
library(strucchange)
library(ggplot2)
library(tseries)
library(data.table)
library(digest)
library(Kendall)

#Testing with a toy dataset 
data1 <- read.csv("C:\\Users\\sunny\\Downloads\\Nile_DataSet.csv", header=TRUE, stringsAsFactors=FALSE)
Data <- data1$x
#Version checking

ds_series <- function(Data, frequency = 52, break_level = 0.05, mean_check = 1.5, median_level = 1.5, level_length = 12, conf_level = 0.1, window_len = 12, mk_test = FALSE, Alpha = 0.01, tau = 0.8, plot = FALSE)
{
  #Writing data to another value 
  y <- Data
  
  #Definition of external vectors 
  anomalies <- vector()
  
  #Inital Sanity checks
  if(!is.numeric(frequency)  || !is.numeric(break_level) || !is.numeric(mean_check) || !is.numeric(median_level) || !is.numeric(level_length) || !is.numeric(conf_level) || !is.numeric(window_len))
  {
    stop(print('Value needs to be numeric'))
  }
  
  if(!is.logical(plot))  
  {
    stop(print('Value needs to be boolean'))
  }  
  
  #If the time series contains NA
  if(any(is.na(y)))
  {
    stop(print('Interpolation of time series needed , recommend the zoo package'))
  }
  
  
  #Removing leading and lagging zeroes 
  non_zero_length <- y[min(which(y != 0)) : length(y)]
  
  #Updating the frame
  z <- (length(y) - length(non_zero_length)) + 1 
  y <- y[c(z : length(y))]
  
  #Strucchange values
  mvalue <- NA
  bp <- NA 
  
  tryCatch(
    {
      mvalue = breakpoints( y ~ 1, h = break_level)
    }, error = function(e){bp <<- NA}
  )
  
  if(!(is.na(mvalue$breakpoints)))
  {
    bp = mvalue$breakpoints
  }
  
  if(is.na(bp)){print('Change break value , min segment size must be larger than the number of regressors')}
  
  #Writing the breakpoints 
  t <- vector()
  
  if(is.na(bp))
  {
    t <- c(0, length(y))
  }else
  {
    t <- c(0, bp, length(y))
  }
  
  #Seasonal check for level changes 
  difference_mat <- outer(t,t, "-")
  difference_table <- which(difference_mat == frequency , arr.ind = TRUE)
  difference_table <- data.frame(difference_table)
  
  p1 <- vector()
  p2 <- vector()
  
  if(dim(difference_table)[1] != 0)
  {
    for(l in seq (from = 1 , to = c(dim(difference_table)[1]), by = 1))
    {
      p1 <- c(t[difference_table['row'][l,]], t[difference_table['col'][l,]])
      p2 <- c(p2, p1)
    }
  }
  
  
  t <- t[!(t %in% p2)]
  
  #Median cleaning of breakpoints
  med_flag <- FALSE
  n <- length(y)
  k <- vector()
  i <- 1
  j <- 3
  
  while(i < n)
  {
    while(j < n+1)
    {
      L <- t[i]
      R <- t[j]
      Lr <- t[i+1]
      
      if((is.na(R)))
      {
        med_flag = TRUE
        break
      }
      
      mid1 <- median(y[c(L+1) : c(Lr)])
      mid2 <- median(y[c(Lr + 1)] : (R))
      
      if(mid2 == 0)
      {
        mid2 <- 0.0000001
      }
      
      if((abs(mid1/mid2) > median_level) || abs(mid2/mid1 > median_level))
      {
        k <- c(k , Lr)
        break
      }else
      {
        t <- t[t!=Lr]
      }
    }
    
    i <- i + 1
    j <- j + 1
    L <- t[i]
    
    if(med_flag == TRUE)
    {
      break
    }
      
  }

  i <- vector()
  j <- vector()
  
  #Mean point check
  q <- vector()
  for(x in seq(0, c(length(k)), 1))
  {
    if(x == length(k))
    {
      break
    }else
    {
      before = y[c(k[[x+1]] - 4):c(k[[x+1]]) - 1]
      rownames(before) <- NULL
      
      after = y[c(k[[x+1]] + 1):c(k[[x+1]]) + 4]
      rownames(after) <- NULL
      
      #Getting the mean values
      M1 = mean(before)
      M2 = mean(after)
      
      #Checking the mean ratio
      if((abs(M1/M2)) > mean_check || ((abs(M2/M1)) > mean_check))
      {
        q <- c(q, k[x+1])
      }
      
      
      
    }
  }
  
  #Min level check
  q <- c(0, q, length(y))
  if(length(q != 2))
  {
    breaks <- vector()
    for(i in seq(from = 1, to = c(length(q) - 1), by = 1))
    {
      if(((q[i+1] - q[i] >= level_length)))
      {
        breaks <- c(q[i+1], breaks)
      }
    }
    
    breaks <- unique(sort(breaks))
    breaks <- c(0, breaks, length(y))
    
    if(length(breaks) == 2)
    {
      break
    }else if(c(length(y) - tail(breaks, 2)[1] <= level_length))
    {
      breaks <- breaks[-match(tail(breaks, 2)[1], breaks)]
    }
    
    q <- breaks 
    
  }else
  {
    q <- c(0, q, length(q))
    q <- unique(sort(q))
  }
  
  #Anomaly check
  #Overall logic is to remove the peridic median value from the series 
  #Identify anomalies using the MAD cutoff
  
  #Initalizations
  y_med <- y
  temp1 <- y
  s <- q
  
  window_medians <- vector()
  outliers <- vector()
  outliers_new <- vector()
  anomalies <- vector()
  
  for(i in seq(from = 1, to = c(length(s) - 1), by = 1))
  {
    window_len <- 14
    win_len <- ceiling(length(y[c(s[i]+1) : c(s[i+1])])/window_len)
    
    for(w in seq(1, (win_len)))
    {
      if(s[i] + window_len >= s[i+1])
      {
        window_len <- s[i+1] - s[i]
      }
      
      med_1 <- median(y_med[c(s[i] + 1) : c(s[i] + window_len)])
      y_med[c(s[i] + 1) : c(s[i] + window_len)] <- med_1
      
      #Storing the medians
      window_medians <- c(window_medians, rep(med_1, window_len))
      
      #Updating length 
      s[i] <- s[i] + window_len
      
    }
  }
  
  #Adjusting window median length
  win_len_new <- (length(window_medians) - length(y_med)) + 1
  window_medians <- window_medians[c(win_len_new : length(window_medians))]
  
  #Subtracting the median values 
  temp1 <- y - y_med
  
  #First pass identifying outliers that are not seasonal in nature 
  outlier_deason <- vector()
  for(i in seq(from = 1, to = c(length(q) - 1), by = 1))
  {
    med_1 <- mad(temp1[c(q[i] + 1) : c(q[i+1])], center = median(c(temp1[c(q[i] + 1) : c(q[i+1])])), constant = 1.4)
    median_val <- median(c(temp1[c(q[i] + 1) : c(q[i+1])]))
    
    for(p in seq(1, (length(temp1[c(q[i] + 1) : c(q[i+1])]))))
    {
      if(temp1[c(q[i] + 1) : c(q[i+1])][p] > c(median_val + (med_1 * conf_level)) || temp1[c(q[i] + 1) : c(q[i+1])][p] < c(median_val - (med_1 * conf_level)))
      {
        #Writing the anomalies out
        outlier_deason <- c(outlier_deason, c(p+q[i]))
      }
    }
  }
  
  #52 Frequency check
  dif_mat_new <- outer(outlier_deason, outlier_deason, '-')
  graph_new <- which(dif_mat_new == frequency, arr.ind = TRUE)
  graph_new <- data.frame(graph_new)
  
  pairs_1_new <- vector()
  pairs_new <- vector()
  
  if(dim(graph_new)[1] != 0)
  {
    for(t in seq(from = 1 , to = c(dim(graph_new)[1]), by = 1))
    {
      pairs_new <- c(outlier_deason[graph_new['row'][t,]], outlier_deason[graph_new['col'][t,]])
      pairs_1_new <- c(pairs_1_new, pairs_new)
    }
  }
  
  #No outliers detected
  outlier_deason <- outlier_deason[(outlier_deason %in% pairs_1_new)]
  
  
  #################
  anomalies <- sort(c(outlier_deason))
  
  for(i in seq(from = 1, to = c(length(q) - 1), by = 1))
  {
    med_1 <- mad(temp1[c(q[i]+1) : c(q[i+1])], center = median(c(temp1[c(q[i]) : c(q[i+1])])), constant = 1.4)
    median_val <- median(c(temp1[c(q[i]+1) : c(q[i+1])]))
    
    #Extracting the anomalies 
    anom_new <- vector()
    anom_new <- anomalies[between(anomalies, q[i], q[i+1])]
    
    if(i == 1)
    {
      anom_new <- anom_new - c(q[i])
    }else
    {
      anom_new <- anom_new - c(q[i] - 1)
    }
    
    
    for(p in seq(1, (length(anom_new))))
    {
      if(i == 1)
      {
        if(is.na(temp1[c(q[i] + 1) : c(q[i+1])][anom_new[p]]))
        {
          break
        }
        
        if(temp1[c(q[i] + 1) : c(q[i+1])][anom_new[p]] > 0)
        {
          temp1[c(q[i] + 1) : c(q[i+1])][anom_new[p]] <- median_val + (med_1 * conf_level)
        }else
        {
          temp1[c(q[i] + 1) : c(q[i+1])][anom_new[p]] <- median_val - (med_1 * conf_level)
        }
      }else
      {
        if(is.na(temp1[c(q[i] + 1) : c(q[i+1])][anom_new[p]]))
        {
           break
        }
        
        if(is.na(temp1[c(q[i] + 1) : c(q[i+1])][anom_new[p+1]]))
        {
          break
        }
        
        if(temp1[c(q[i] + 1) : c(q[i+1])][anom_new[p+1]] > 0)
        {
          temp1[c(q[i] + 1) : c(q[i+1])][anom_new[p+1]] <- median_val + (med_1 * conf_level)
        }else
        {
          temp1[c(q[i] + 1) : c(q[i+1])][anom_new[p+1]] <- median_val - (med_1 * conf_level)
        }

      }
    }
    
  }  
    
    
  final_x <- (temp1 + window_medians)
  y <- final_x
  
  anomalies <- sort(unique(anomalies))
  
  #Smoothening the series
  k <- vector()
  k <- q
  k <- sort(k)
  
  trend_line <- vector()
  if(length(k) == 0)
  {
    q = lowess(Data)
    trend_line = c(trend_line, q$y) 
  }else
  {
    trend_line <- vector()
    pointer <- vector()
    
    for(i in seq(from = 0, to = c(length(k) - 2), by = 1))
    {
      if(i == 0)
      {
       pointer <- 0
      }else
      {
        pointer <- 1
      }
      
      v <- y[(c(k[i+1]+pointer) : k[i+2])]
      q = lowess(v)
      trend_line = c(trend_line, q$y)
    }
  }
  
  #Detrending the series 
  de_trend <- c(Data - trend_line)
  
  #Trasnforming into a series
  ts.QTY1 = ts(data = as.vector(t(de_trend)), frequency = frequency)
  
  decomposed <- NA
  tryCatch(
    {
      decomposed <- stl(ts.QTY1, s.window = 'periodic')
    }, error = function(e){ts_QTY1 <<- NULL}
  )
  
  if(length(decomposed) != 1)
  {
    seasonal <- decomposed$time.series[,1]
    trend <- decomposed$time.series[,2]
    remainder <- decomposed$time.series[,3]
    
    #Removing seasonality 
    ts_QTY1 <- ts.QTY1 - seasonal
  }else
  {
    seasonal <- NULL
    trend <- NULL
    remainder <- NULL
  }
  
  #Plotting function for anomalies
  if(plot)
  {
    #Borrows from twitters anomaly visualization
    
    #Anomalies plot 
    b <- Data
    
    #Cleaning the data
    temp <- data.frame(b)
    
    colnames(temp) <- NULL
    rownames(temp) <- NULL
    
    temp$name <- rownames(temp)
    
    temp <- temp[,c(2,1)]
    
    #Getting time stamps as single columns 
    #Setting row and column names
    
    if(any((names(temp) == c("timestamp", "count")) == FALSE))
    {
      colnames(temp) <- c("timestamp", "count")
    }
    
    #Changing for a numeric type 
    temp$timestamp <- as.numeric(temp$timestamp)
    
    #Plotting based on anomalies
    anomalies_data <- temp[temp$timestamp %in% anomalies, ]
    
    #Final plot function
    visual <- ggplot(temp, aes(timestamp, count, group = 1), title = 'Anomalies') + geom_line(size = 0.5) +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            plot.background = element_rect(fill = 'white'),
            plot.background = element_rect(fill = 'white'),
            panel.background = element_rect(fill = ' white'), 
            axis.line = element_line(colour = 'grey'))
    plot(visual)
    
    return(plot = visual)
    
  }
  
  newList <- list('anomalies' = anomalies, 'adjusted_series' = y, 'trend_line' = trend_line, 'ds_series' = ts_QTY1, 
                  'Trend' = trend, 'Seasonality' = seasonal,'residual' = remainder, 'breakpoints' = k)
  
  return(newList)
  
  
}

