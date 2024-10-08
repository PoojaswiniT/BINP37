
#####

split_ID = function(transcript_matrix){ #Splitting the ID
  # Split based on dot
  id_split = strsplit(as.character(transcript_matrix$transcript_id), "\\.")

  # Extract values before and after the dot
  transcript_matrix$transcript_id = sapply(id_split, `[`, 1)
  transcript_matrix$version = sapply(id_split, `[`, 2)

  # Re-arrange the order
  transcript_matrix = transcript_matrix %>%
    dplyr::select(transcript_id, version, everything())
  return(transcript_matrix)
}

#####