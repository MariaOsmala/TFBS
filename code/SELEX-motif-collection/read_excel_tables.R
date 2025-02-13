# utility function to get rle as a named vector
vec_rle <- function(v){
  temp <- rle(v)
  out <- temp$values
  names(out) <- temp$lengths
  return(out)
}

# utility function to map table with their columns/rows in a bigger table
make_df_index <- function(v){
  table_rle <- vec_rle(v)
  divide_points <- c(0,cumsum(names(table_rle)))
  table_index <- map2((divide_points + 1)[1:length(divide_points)-1],
                      divide_points[2:length(divide_points)],
                      ~.x:.y)
  return(table_index[table_rle])
}

# split a large table in one direction if there are blank columns or rows
split_direction <- function(df,direction = "col"){
  if(direction == "col"){
    col_has_data <- unname(map_lgl(df,~!all(is.na(.x))))
    df_mapping <- make_df_index(col_has_data)
    out <- map(df_mapping,~df[,.x])
  } else if(direction == "row"){
    row_has_data <- df %>% 
      mutate_all(~!is.na(.x)) %>%
      as.matrix() %>% 
      apply(1,any)
    df_mapping <- make_df_index(row_has_data)
    out <- map(df_mapping,~df[.x,])
  }
  return(out)
}

# sometimes different tables are close to each other without empty space in between
# in thoses cases we need to manually insert rows and columns in the dataframe in
# order to make split_df() function to work
insert_row <- function(df,at,on = "below"){
  df <- as_tibble(df)
  if(!is.integer(at)) stop("at should be a integer")
  if(on == "below"){
    
  } else if(on == "above"){
    at <- at - 1
  } else {
    stop("on should be either 'below' or above")
  }
  
  if(at > 0){
    df1 <- df[1:at,]
    df2 <- df[at+1:nrow(df)-at,]
    out <- df1 %>%
      add_row() %>%
      bind_rows(df2)  
  } else if(at == 0) {
    df1 <- df[0,]
    df2 <- df[at+1:nrow(df)-at,]
    out <- df1 %>%
      add_row() %>%
      bind_rows(df2)
  } else {
    stop("Positon 'at' should be at least 1")
  }
  return(out)
}

insert_column <- function(df,at,on = "right"){
  df <- as_tibble(df)
  if(!is.integer(at)) stop("at should be a integer")
  if(on == "right"){
    
  } else if(on == "left"){
    at <- at - 1
  } else {
    stop("on should be either 'right' or left")
  }
  
  if(at > 0){
    df1 <- df[,1:at]
    df1[,1+ncol(df1)] <- NA
    df2 <- df[,at+1:ncol(df)-at]
    out <- df1 %>%
      bind_cols(df2)  
  } else if(at == 0) {
    df1 <- df[,0]
    df1[,1] <- NA
    df2 <- df
    out <- df1 %>%
      bind_cols(df2)
  } else {
    stop("Positon 'at' should be at least 1")
  }
  return(out)
}


# split a large table into smaller tables if there are blank columns or rows
# if you still see entire rows or columns missing. Please increase complexity
split_df <- function(df,showWarnig = TRUE,complexity = 1){
  if(showWarnig){
    warning("Please don't use first row as column names.")
  }
  
  out <- split_direction(df,"col")
  
  for(i in 1 :complexity){
    out <- out %>%
      map(~split_direction(.x,"row")) %>%
      flatten() %>%
      map(~split_direction(.x,"col")) %>%
      flatten()
  }
  return(out)
  
}

# set first row as column names for a data frame and remove the original first row
set_1row_colname <- function(df){
  colnames(df) <- as.character(df[1,])
  out <- df[-1,]
  return(out)
}

#display the rough shape of table in a sheet with multiple tables
display_table_shape <- function(df){
  colnames(df) <- 1:ncol(df)
  
  out <- df %>%
    map_df(~as.numeric(!is.na(.x))) %>%
    gather(key = "x",value = "value") %>%
    mutate(x = as.numeric(x)) %>%
    group_by(x) %>%
    mutate(y = -row_number()) %>%
    ungroup() %>%
    dplyr::filter(value == 1) %>%
    ggplot(aes(x = x, y = y,fill = value)) +
    geom_tile(fill = "skyblue3") +
    scale_x_continuous(position = "top") +
    theme_void() +
    theme(legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))
  return(out)
}