library(diagram)
openplotmat()
elpos <- coordinates(c(2, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 2, 1))
fromto <- matrix(ncol = 2, byrow = TRUE,
                 data = c(1, 3, 1, 4, 2, 5, 2, 6, 
                          3, 8, 4, 9, 5, 10, 6, 11,
                          8, 13, 9, 14, 10, 15, 11, 16,
                          13, 18, 14, 19, 15, 20, 16, 21,
                          18, 23, 19, 24, 20, 25, 21, 26,
                          23, 28, 24, 29, 25, 30, 26, 31,
                          28, 32, 29, 32, 30, 33, 31, 33,
                          32, 34, 33, 34))
nr <- nrow(fromto)
arrpos <- matrix(ncol = 2, nrow = nr)
for (i in 1:nr)
  arrpos[i, ] <- straightarrow (to = elpos[fromto[i, 2], ],
                                from = elpos[fromto[i, 1], ],
                                lwd = 2, arr.pos = 0.65, arr.length = 0.3)

textellipse(elpos[1,], 0.025, cex = 0.8, lab = "tissues")
textellipse(elpos[2,], 0.025, cex = 0.8, lab = "ripening")
textellipse(elpos[3,], 0.025, cex = 0.8, lab = "POS")
textellipse(elpos[4,], 0.025, cex = 0.8, lab = "NEG")
textellipse(elpos[5,], 0.025, cex = 0.8, lab = "POS")
textellipse(elpos[6,], 0.025, cex = 0.8, lab = "NEG")
textrect(elpos[7,], 0.45, 0.015, cex = 0.8, lab = "Detected features in MS1 data")
textellipse(elpos[8,], 0.025, cex = 0.8, lab = "3558")
textellipse(elpos[9,], 0.025, cex = 0.8, lab = "2158")
textellipse(elpos[10,], 0.025, cex = 0.8, lab = "2395")
textellipse(elpos[11,], 0.025, cex = 0.8, lab = "1858")
textrect(elpos[12,], 0.45, 0.015, cex = 0.8, lab = "Non-isotopic MS1 features")
textellipse(elpos[13,], 0.025, cex = 0.8, lab = "2790")
textellipse(elpos[14,], 0.025, cex = 0.8, lab = "1713")
textellipse(elpos[15,], 0.025, cex = 0.8, lab = "1728")
textellipse(elpos[16,], 0.025, cex = 0.8, lab = "1458")
textrect(elpos[17,], 0.45, 0.015, cex = 0.8, lab = "Non-isotopic MS1 features with at least one MS2 spectra")
textellipse(elpos[18,], 0.025, cex = 0.8, lab = "716")
textellipse(elpos[19,], 0.025, cex = 0.8, lab = "525")
textellipse(elpos[20,], 0.025, cex = 0.8, lab = "515")
textellipse(elpos[21,], 0.025, cex = 0.8, lab = "458")
textrect(elpos[22,], 0.45, 0.015, cex = 0.8, lab = "Annotated features")
textellipse(elpos[23,], 0.025, cex = 0.8, lab = "231")
textellipse(elpos[24,], 0.025, cex = 0.8, lab = "165")
textellipse(elpos[25,], 0.025, cex = 0.8, lab = "178")
textellipse(elpos[26,], 0.025, cex = 0.8, lab = "121")
textrect(elpos[27,], 0.45, 0.015, cex = 0.8, lab = "Unique compounds")
textellipse(elpos[28,], 0.025, cex = 0.8, lab = "201")
textellipse(elpos[29,], 0.025, cex = 0.8, lab = "139")
textellipse(elpos[30,], 0.025, cex = 0.8, lab = "171")
textellipse(elpos[31,], 0.025, cex = 0.8, lab = "107")
textellipse(elpos[32,], 0.025, cex = 0.8, lab = "188")
textellipse(elpos[33,], 0.025, cex = 0.8, lab = "168")
textellipse(elpos[34,], 0.025, cex = 0.8, lab = "216")

