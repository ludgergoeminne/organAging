preprocess_SVG_file <- function(file.in, suffix = "", replace.Greek.letters = TRUE){
  # file.in <- "//lispnas1/common/Users/goeminne/AmyloAge30/results/correlations_aggregates_alpha_crystallins.svg"
  file.out <- gsub("(.+)\\.svg", paste0("\\1", suffix, ".svg"), file.in) # suffix <- "_new"
  
  xml.obj <- readChar(file.in, file.info(file.in)$size)
  
  # After it has been opened in Inkscape:
  # xml.obj <- gsub('\n.*?textLength=\".*?\"', '', xml.obj, perl = TRUE)
  # xml.obj <- gsub('\n.*?lengthAdjust=\".*?\"', '', xml.obj, perl = TRUE)
  
  # Before it has been opened in Inkscape:
  xml.obj <- gsub(" textLength='.*?'", "", xml.obj, perl = TRUE)
  xml.obj <- gsub(" lengthAdjust='.*?'", "", xml.obj, perl = TRUE)
  xml.obj <- gsub(" clip-path='url(.*?)'", "", xml.obj, perl = TRUE)
  
  if(replace.Greek.letters){
    xml.obj <- gsub("[A|a]lpha", "α", xml.obj, perl = TRUE)
    xml.obj <- gsub("[B|b]eta", "β", xml.obj, perl = TRUE)
    xml.obj <- gsub("[D|d]elta", "δ", xml.obj, perl = TRUE)
    xml.obj <- gsub("DELTA", "Δ", xml.obj, perl = TRUE)
    xml.obj <- gsub("rho", "ρ", xml.obj, perl = TRUE)
    xml.obj <- gsub("cirρsis", "cirrhosis", xml.obj, perl = TRUE) # Don't add rho for cirrhosis!
    xml.obj <- gsub("RHO", "Ρ", xml.obj, perl = TRUE)
  }
  
  writeChar(xml.obj, file.out)
}