fml_maker <- function(class, C, db){
  case_when(
    class == "CAR" ~ paste0("C", C + 7, "H", C*2 - (2*db) + 13, "NO4"),
    
    class == "LPC" | class == "LPC-Na" ~ paste0("C", C + 8, "H", C*2 - (2*db) + 18, "NO7P"),
    class == "LPI" ~ paste0("C", C + 9, "H", C*2 - (2*db) + 17, "O12P"),
    
    class == "PA" ~ paste0("C", C + 3, "H", C*2 - (2*db) + 5, "O8P"),
    class == "mPA"~ paste0("C", C + 4, "H", C*2 - (2*db) + 7, "O8P"),
    class == "PC" ~ paste0("C", C + 8, "H", C*2 - (2*db) + 16, "NO8P"),
    class == "PE" ~ paste0("C", C + 5, "H", C*2 - (2*db) + 10, "NO8P"),
    class == "PG" ~ paste0("C", C + 6, "H", C*2 - (2*db) + 11, "O10P"),
    class == "PI" ~ paste0("C", C + 9, "H", C*2 - (2*db) + 15, "O13P"),
    class == "PS" ~ paste0("C", C + 6, "H", C*2 - (2*db) + 10, "NO10P"),
    
    class == "CER" | class == "CER-FA" ~ paste0("C", C, "H", C*2 - (2*db) + 1, "NO3"),
    class == "CERhex" | class == "CERhex-FA" ~ paste0("C", C + 6, "H", C*2 - (2*db) + 11, "NO8"),
    class == "CERlact" ~ paste0("C", C + 12, "H", C*2 - (2*db) + 21, "NO13"),

    class == "MGDG" | class == "MGDG-Na" ~ paste0("C", C + 9, "H", C*2 - (2*db) + 14, "O10"),
    class == "DGDG" | class == "DGDG-Na" ~ paste0("C", C + 15, "H", C*2 - (2*db) + 24, "O15"),
    
    class == "DAG" ~ paste0("C", C + 3, "H", C*2 - (2*db) + 4, "O5"),
    class == "TAG" ~ paste0("C", C + 3, "H", C*2 - (2*db) + 2, "O6")
  )
}

.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}

extract_best <- function(x){
  x <- x[x==max(x)]
  out <- paste(names(x), collapse = "-")
  out
}
