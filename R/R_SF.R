

#' RSF_convert Function
#'
#' Do basic transformation on *.rsf file (rsf = wind ReSource File): Crop, convert in 12 sector, create new height levels, compute and export shear table.
#' @param rsf_file1: Input *.rsf file name (name should be surounded by "").
#' @param output_name: Output *.rsf file name " (name should be surounded by "").
#' @param crop: Should the rsf be croped ?  TRUE/FALSE. Default value = FALSE.
#' @param layout_csv: Input *.csv file name (name should be surounded by "")  containing the WTG positions (X,Y).Default value = NULL. This input is used to crop and benchmark the production (12 Sectors versus 36 sectors).
#' @param crop_buffer: Integer value: Distance in meters around the WTG postions that should be crop (default value = 700m)
#' @param twelve_S: Should the rsf be converted into 12 sectors ?  TRUE/FALSE. Default value = FALSE.
#' @param rsf_file2: Input *.rsf file name (name should be surounded by ""). Second rsf file used for shear and 3 dimensionals rsf creation. Default value = NULL.
#' @param shear_out: Should the shear should be extracted and exported (csv file) ?  TRUE/FALSE. Default value = FALSE.
#' @param three_D: Should a three dimensional rsf created  ? TRUE/FALSE. Default value = FALSE.
#' @param layer_H: vector containing  the level of the 3 dimensional rsf (used if three_D = TRUE). if not filled and three_D = TRUE, the level multiple of 10 would be created. Default value = NULL.
#' @keywords rsf, read, convert
#' @export
#' @examples
#' crop the "Aldermyrberget 166m.rsf" 700m (default value) around the layout "layout V14A.csv":
#' RSF_convert(rsf_file1 = "Aldermyrberget 166m.rsf",crop = TRUE,layout_csv = "layout V14A.csv" )
#'
#' Crop 500m and convert into 12 sector "Aldermyrberget 166m.rsf"
#' RSF_convert(rsf_file1 = "Aldermyrberget 166m.rsf",crop = TRUE,layout_csv = "layout V14A.csv", crop_buffer = 500, twelve_S = TRUE )
#'
#' create a 3 dimensional rsf on the levels (112,137,142).
#' RSF_convert(rsf_file1 = "Aldermyrberget 166m.rsf",three_D = TRUE, rsf_file2 = "Aldermyrberget 103m.rsf", layer_H = c(112,137,142))
#'
#' create a 3 dimensional rsf (default level) and export the shear table.
#' RSF_convert(rsf_file1 = "Aldermyrberget 166m.rsf",three_D = TRUE, rsf_file2 = "Aldermyrberget 103m.rsf",shear_out = TRUE)
#'
#' crop, convert to 12 sector, create a 3 dimensional rsf (containing the 112,137 ad 142m levels) and also export the shear table.
#' RSF_convert(rsf_file1 = "Aldermyrberget 166m.rsf",crop = TRUE,layout_csv = "layout V14A.csv",twelve_S = TRUE,layer_H = c(112,137,142), three_D = TRUE, rsf_file2 = "Aldermyrberget 103m.rsf",shear_out = TRUE)
#'
RSF_convert <- function(rsf_file1,
                         output_name = rsf_file1,
                         crop = FALSE,
                         layout_csv = NULL,
                         crop_buffer = 700,
                         twelve_S = FALSE,
                         rsf_file2 = NULL,
                         shear_out = FALSE,
                         three_D = FALSE,
                         layer_H = NULL) {

        library("data.table")

        # read and parse data -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        rsf_H1 <- Read_RSF(rsf_file1)
        H1 <- unique(rsf_H1[, Height])

        if (!is.null(rsf_file2)) {
            rsf2 <- TRUE
            rsf_H2 <- Read_RSF(rsf_file2)
            H2 <- unique(rsf_H2[, Height])

            if (length(rsf_H1) != length(rsf_H2)) {
                stop("The inputs rsf file dimension do not match, Check the coherence of the two *.rsf files")
            } else if (abs(H1 - H2) < 20) {
                warning(paste("the heigth difference is below 20m", H1, H2))
            }
        } else {
            rsf2 <- FALSE
        }



        if (!is.null(layout_csv)) {
            layout <- read.csv(layout_csv)
            names(layout) <- c("X", "Y")
            layout <- layout[complete.cases(layout), ]

            if (length(layout) != 2) {
                stop(
                    "Check that the layout file is a *.csv file with two columns using ",
                    " as delimiter"
                )
            } else if (min(layout[, 1]) < rsf_H1[, min(X)] |
                       max(layout[, 1]) > rsf_H1[, max(X)] |
                       min(layout[, 2]) < rsf_H1[, min(Y)] |
                       max(layout[, 2]) > rsf_H1[, max(Y)]) {
                stop("At least one position of the layout seems outside the rsf area: Check your positions and the coherence of the coordinate systems")
            }
        }


        # Crop files -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*


        if (crop == TRUE) {
            Xmax <- max(layout[, 1] + crop_buffer)
            Xmin <- min(layout[, 1] - crop_buffer)
            Ymax <- max(layout[, 2] + crop_buffer)
            Ymin <- min(layout[, 2] - crop_buffer)

            rsf_H1 <- Crop_RSF(rsf_H1, Xmax, Xmin, Ymax, Ymin)

            if (rsf2 == TRUE) {
                rsf_H2 <- Crop_RSF(rsf_H2, Xmax, Xmin, Ymax, Ymin)
            }

            output_name <- paste("[crop_", crop_buffer, "]", output_name, sep = "")
        }

        # 12 sectors conversion if required -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        if (twelve_S == TRUE) {
            rsf_H1_12 <- S36_to_S12_RSF(rsf_H1)

            bench_table <- Bensh_RSF(rsf_H1, rsf_H1_12, layout, Gen_power_curve())

            Ratio <- bench_table[, mean(ratio12_36)]

            write.csv(bench_table,
                      paste("benchmark_table_", output_name, ".csv", sep = ""))

            pdf(file = paste("benchmark_graph_", output_name, ".pdf", sep = ""))
            hist(
                bench_table$ratio12_36,
                10,
                col = "green",
                main = "AEP ratio: 12 sectors versus 36 sectors",
                xlab = "Ratio",
                breaks = seq(0.9, 1.1, 0.005)
            )
            abline(v = mean(bench_table$ratio12_36), lwd = 3)
            rug(bench_table$ratio12_36)
            dev.off()

            rsf_H1 <- rsf_H1_12
            rm(rsf_H1_12)

            if (rsf2 == TRUE) {
                rsf_H2_12 <- S36_to_S12_RSF(rsf_H2)
                rsf_H2 <- rsf_H2_12
                rm(rsf_H2_12)
            }

            output_name <-
                paste("[12S_",
                      round(Ratio * 100 - 100, 1),
                      "]",
                      output_name,
                      sep = "")

        }


        # Shear table -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        if (shear_out == TRUE) {
            shear_table_file_name <- paste("[shear_table]", output_name, ".csv", sep = "")
            shear_table <- Shear_RSF(rsf_H2, rsf_H1)
            ShearTable(shear_table, shear_table_file_name)  #create shear file
        }


        # 3D rsf -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        if (three_D == TRUE) {
            # creation of the vector containing the different calculation heights
            H1 <- unique(rsf_H1[, Height])
            H2 <- unique(rsf_H2[, Height])

            if (is.null(layer_H)) {     #if layer_H is not define, create a layer every 10 meters
                layer_H <- seq(0, 200, 10)
                layer_H <- layer_H[layer_H > min(H1, H2) & layer_H < max(H1, H2)]
            }

            # creation of the 3 dimensionals rsf
            rsf_H1 <- Interpol_RSF(rsf_H1, rsf_H2, layer_H)
            output_name <- paste("[3D]", output_name, sep = "")
        }


        # Write to file -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

        Write_RSF(rsf_H1 , output_name)

    }


#' Read_RSF Function
#'
#' This function allows you to read an *.rsf file. it reconises if it is 36 or 12 sectors, load it accordingly, with columns names.
#' rsf = wind ReSource File. it is a wind mapping output from wind flow model. the new object will be a data.table
#' @param Input_file: "*.rsf" file ()
#' @keywords rsf, read
#' @export
#' @examples
#' load into memory an rsf file into the object rsf_166:
#' rsf_166 <- Read_RSF("RSF-windresource-CFD_Aldermyrberget 166m.rsf")
Read_RSF <- function(Input_file) {
    ## function read an *.rsf file. it reconises if it is 36 or 12 sectors, load it accordingly, with columns names

    library("readr")              ## use for the read_fwf function much more faster than read.fwf
    library("data.table")         ## use forthe output format

    con <- file(Input_file, "r")
    line <-
        readLines(con, 1)      ## sample line to check if it is a 36 sectors or 12 sectors file
    close(con)

    if (nchar(line) > 400) {
        ## it means it should be a 36 sectors (actually 36S should be exactly 540)
        RSF_Col_Format <- c(10, 10, 10, 8, 5, 6, 5, 15, 3, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4,
                            4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4,
                            4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4,
                            4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4,
                            4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4,
                            4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5)


        ## creation of colomn name vector
        RSF_Col_Name <- c("Label",
                          "X",
                          "Y",
                          "Z",
                          "Height",
                          "Aave",
                          "k",
                          "Blank",
                          "Sector")
        RSF_Col_Name <- c(RSF_Col_Name, paste(c("F36_", "A36_", "k36_"), rep(seq(0, 350, 10), each = 3), sep = ""))



    } else{
        ##if less than 400, then it means it should be a 12 sectors (actually 12S should be exactly 228)
        RSF_Col_Format <- c(10, 10, 10, 8, 5, 6, 5, 15, 3, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4,
                            4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4,
                            4, 5, 4, 4, 5)

        RSF_Col_Name <-
            c("Label",
              "X",
              "Y",
              "Z",
              "Height",
              "Aave",
              "k",
              "Blank",
              "Sector")
        RSF_Col_Name <-
            c(RSF_Col_Name, paste(c("F12_", "A12_", "k12_"), rep(seq(0, 330, 30), each = 3), sep = ""))
    }

    RSF_table <-
        data.table(read_fwf(
            Input_file,
            fwf_widths(RSF_Col_Format, col_names = RSF_Col_Name)
        ))

    if (is.na(RSF_table$Label[1])) {RSF_table$Label <- "RRR"} # to avoid issue in other functions

    RSF_table <- RSF_table[order(RSF_table[, Y], RSF_table[, X])]   #Order the file: Y and X. Windpro is not always coherent when creating an rsf file from WasP.

    RSF_table

}

#' Crop_RSF Function
#'
#' Crop an *.rsf file around a rectangle define by Xmax, Xmin, Ymax, Ymin.
#' @param rsf: "*.rsf" file ()
#' @param Xmax: X maximum limit
#' @param Xmax: X minimum limit
#' @param Ymax: Y maximum limit
#' @param Ymin: Y minimum limit
#' @keywords rsf, crop
#' @export
#' @examples
#' load into memory an rsf file into the object rsf_166:
#' rsf_166 <- Read_RSF("RSF-windresource-CFD_Aldermyrberget 166m.rsf")
Crop_RSF <- function(rsf, Xmax, Xmin, Ymax, Ymin) {
        # (A,B,C,D, could be read as "Xmax, Xmin, Ymax, Ymin", or just "Layout_Table, Buffer in [m]")

        library("data.table")

        if (!is(rsf, "data.table")) {
            stop("the Object pass in Crop_RSF is not a data.table")
        } else{
            Sub_rsf <- rsf[X <= Xmax & X >= Xmin & Y <= Ymax & Y >= Ymin]
        }
    }

#' S36_to_S12_RSF Function
#'
#' Tranform rsf from 36 sectors to 12. it combines the neighboring sectors 3 by 3 by doing a ponderate average.
#' @param RSF36: 36 sectors rsf (already loaded in memory by the use fo Read_RSF() ).
#' @keywords rsf, 36, convertion
#' @export
#' @examples
#' convert rsf87_36S from 36 to 12 sector and store the result in rsf87_12S:
#' rsf87_12S <- S36_to_S12_RSF(rsf87_36S)
S36_to_S12_RSF <- function(RSF36) {
    # Combine the neighbouring sectors 3 by 3 to end-up with a 12 sectors RSF.
    # For the moment, it uses columns names as reference: beter for understanding, but this mean that the input RSF should have the corrects columns names
    # Could maybe be change by using columns numbers & and using a loop on the 12 sectors (it would be less easy to understand)

    if (length(RSF36) != 117) {
        stop("The object passed in S36_to_S12_RSF function is not a 36 Sector RSF (lenght of RSF36 should be 117)")
    } else if (!is(RSF36, "data.table")) {
        stop("the Object pass in S36_to_S12_RSF function is not a data.table")
    }


    ## creating the RSF 12 sector
    RSF12 <- data.table(RSF36)


    ## Create new columns with 12 sectors values (Frequency, Average Wind Speed, k)
    RSF12[, F12_0 := F36_350 + F36_0 + F36_10]
    RSF12[, A12_0 := round((A36_350 * F36_350 + A36_0 * F36_0 + A36_10 * F36_10) /
                               F12_0, 0)]
    RSF12[, k12_0 := round((k36_350 * F36_350 + k36_0 * F36_0 + k36_10 * F36_10) /
                               F12_0, 0)]

    RSF12[, F12_30 := F36_20 + F36_30 + F36_40]
    RSF12[, A12_30 := round((A36_20 * F36_20 + A36_30 * F36_30 + A36_40 *
                                 F36_40) / F12_30, 0)]
    RSF12[, k12_30 := round((k36_20 * F36_20 + k36_30 * F36_30 + k36_40 *
                                 F36_40) / F12_30, 0)]

    RSF12[, F12_60 := F36_50 + F36_60 + F36_70]
    RSF12[, A12_60 := round((A36_50 * F36_50 + A36_60 * F36_60 + A36_70 *
                                 F36_70) / F12_60, 0)]
    RSF12[, k12_60 := round((k36_50 * F36_50 + k36_60 * F36_60 + k36_70 *
                                 F36_70) / F12_60, 0)]

    RSF12[, F12_90 := F36_80 + F36_90 + F36_100]
    RSF12[, A12_90 := round((A36_80 * F36_80 + A36_90 * F36_90 + A36_100 *
                                 F36_100) / F12_90, 0)]
    RSF12[, k12_90 := round((k36_80 * F36_80 + k36_90 * F36_90 + k36_100 *
                                 F36_100) / F12_90, 0)]

    RSF12[, F12_120 := F36_110 + F36_120 + F36_130]
    RSF12[, A12_120 := round((A36_110 * F36_110 + A36_120 * F36_120 + A36_130 *
                                  F36_130) / F12_120,0)]
    RSF12[, k12_120 := round((k36_110 * F36_110 + k36_120 * F36_120 + k36_130 *
                                  F36_130) / F12_120,0)]

    RSF12[, F12_150 := F36_140 + F36_150 + F36_160]
    RSF12[, A12_150 := round((A36_140 * F36_140 + A36_150 * F36_150 + A36_160 *
                                  F36_160) / F12_150,0)]
    RSF12[, k12_150 := round((k36_140 * F36_140 + k36_150 * F36_150 + k36_160 *
                                  F36_160) / F12_150,0)]

    RSF12[, F12_180 := F36_170 + F36_180 + F36_190]
    RSF12[, A12_180 := round((A36_170 * F36_170 + A36_180 * F36_180 + A36_190 *
                                  F36_190) / F12_180,0)]
    RSF12[, k12_180 := round((k36_170 * F36_170 + k36_180 * F36_180 + k36_190 *
                                  F36_190) / F12_180,0)]

    RSF12[, F12_210 := F36_200 + F36_210 + F36_220]
    RSF12[, A12_210 := round((A36_200 * F36_200 + A36_210 * F36_210 + A36_220 *
                                  F36_220) / F12_210,0)]
    RSF12[, k12_210 := round((k36_200 * F36_200 + k36_210 * F36_210 + k36_220 *
                                  F36_220) / F12_210,0)]

    RSF12[, F12_240 := F36_230 + F36_240 + F36_250]
    RSF12[, A12_240 := round((A36_230 * F36_230 + A36_240 * F36_240 + A36_250 *
                                  F36_250) / F12_240,0)]
    RSF12[, k12_240 := round((k36_230 * F36_230 + k36_240 * F36_240 + k36_250 *
                                  F36_250) / F12_240,0)]

    RSF12[, F12_270 := F36_260 + F36_270 + F36_280]
    RSF12[, A12_270 := round((A36_260 * F36_260 + A36_270 * F36_270 + A36_280 *
                                  F36_280) / F12_270,0)]
    RSF12[, k12_270 := round((k36_260 * F36_260 + k36_270 * F36_270 + k36_280 *
                                  F36_280) / F12_270,0)]

    RSF12[, F12_300 := F36_290 + F36_300 + F36_310]
    RSF12[, A12_300 := round((A36_290 * F36_290 + A36_300 * F36_300 + A36_310 *
                                  F36_310) / F12_300,0)]
    RSF12[, k12_300 := round((k36_290 * F36_290 + k36_300 * F36_300 + k36_310 *
                                  F36_310) / F12_300,0)]

    RSF12[, F12_330 := F36_320 + F36_330 + F36_340]
    RSF12[, A12_330 := round((A36_320 * F36_320 + A36_330 * F36_330 + A36_340 *
                                  F36_340) / F12_330,0)]
    RSF12[, k12_330 := round((k36_320 * F36_320 + k36_330 * F36_330 + k36_340 *
                                  F36_340) / F12_330,0)]

    RSF12[, Sector := 12]


    #remove the 36 sectors columns
    RSF12[, paste(c("F36_", "A36_", "k36_"), rep(seq(0, 350, 10), each = 3), sep = "") := NULL]

}

#' Write_RSF Function
#'
#' Write an rsf object from memory to a file (*.rsf). This file can be used with Windpro.
#' @param RSF: rsf object already loaded into memory
#' @param output_file_name: name for the output file. it has to contain the extension .rsf.
#' @keywords rsf, write
#' @export
#' @examples
#' write from memory a rsf (data table) into the file "rsf_80m_Aldermyberget.rsf":
#' Write_RSF(rsf80, "rsf_80m_Aldermyberget.rsf")
Write_RSF <- function(RSF, output_file_name) {
    # write an rsf table ("RSF") in memory to a *.rsf file ("output_file_name")

    library("gdata") ## use for the write.fwf function (write.fwf writes object in *f*ixed *w*idth *f*ormat )

    Num_Sectors <- RSF[, unique(Sector)]
    RSF_Col_Format <- c(c(10, 10, 10, 8, 5, 5, 6, 15, 3), rep(c(4, 4, 5), Num_Sectors))

    # write the results of the 12 sectors RSF in a file (Fixe Wild File)
    write.fwf(
        RSF,
        output_file_name,
        colnames = FALSE,
        width = RSF_Col_Format,
        sep = "",
        na = "40"
    )    #na need a value (other wise windpro will trow an error)

}

#' Shear Function
#'
#' Create a Shear map (rsf format) based on the two rsf files (rsf_H1,rsf_H2). It is recommended to have abs(H1-H2)>20m. the two rsf shoud be coherent (domain, and number of setors). The output is also in a rsf format.
#' @param rsf_H1: rsf object already loaded into memory
#' @param rsf_H2: rsf object already loaded into memory
#' @keywords rsf, shear
#' @export
#' @examples
#' creat a shear map from the two rsf rsf87 and rsf137:
#' rsf_shear <- Shear_RSF(rsf87,rsf137)
Shear_RSF <- function(rsf_H1, rsf_H2) {
    #creation of a "Shear" RSF, that can be visualised with windPro (need column manipulation and correcttion factor to revert back the A --> Vave correction)

    H1 <- unique(rsf_H1[, Height])
    H2 <- unique(rsf_H2[, Height])
    N_Sector <- unique(rsf_H2[, Sector])

    if (N_Sector == 12) {
        sector_size = 30
    } else if (N_Sector == 36) {
        sector_size = 10
    }

    Shear <-
        rsf_H1[, .(Label, X, Y, Z, Height, Aave, k, Blank, Sector)]  #  initialization shear <= first columns of the rsf

    Shear[, ':=' (
        Height = 999,
        Aave = log(rsf_H1[, Aave] / rsf_H2[, Aave]) / log(H1 / H2),
        Sector = 0 )] #change of the All sector (999) part
    Buffer <- Shear[, .(Label, X, Y, Z, Height, Aave, k, Blank, Sector)]

    for (i in 1:N_Sector) {
        # fil sectorwise wind shear

        j <-
            i * 3 + 7 # column number in the rsf data.tables for the concerned sector

        Buffer[, "Aave" := (log(rsf_H1[, j + 1, with = FALSE] / rsf_H2[, j + 1, with = FALSE]) / log(H1 / H2)), with = FALSE]  # note that here because of the With = FALSE columns name has to be passed as a Charactere string : with " "
        Buffer[, Height := (i - 1) * sector_size]
        Shear <- rbind(Shear, Buffer)

    }

    #Shear[,k:=2]                                   #  set k to constant = 2 (rsf "Hack")
    #correction_Factor = 1/gamma(1+1/2)
    #Shear[,Aave:=round(Aave*correction_Factor,2)]  #  reverse function from Aave to Vave ...but with "shear values" as Vave (rsf "Hack")

    names(Shear)[names(Shear) == "Aave"] <- "ShearAve"
    names(Shear)[names(Shear) == "Height"] <- "Direction"
    Shear <- Shear[complete.cases(Shear)]
}

#' Prod_RSF Function
#'
#' Compute the production [MWh] for one turbine type (power_curve) at the unique position (Xp, Yp) based on the rsf (rsf).
#' @param rsf: rsf object already loaded into memory
#' @param Xp: turbine position (x)
#' @param Yp: turbine position (y)
#' @param power_curve: power curve used to compute production. the format of power_curve is a data frame with two columns (Wind speed [m/s], power [kW]). can be generate by Gen_power_curve() function.
#' @keywords rsf, production
#' @export
#' @examples
#' Compute AEP for V117 at the position (3756122,1625482) and store it into AEP_T1:
#' pc <- Gen_power_curve("V117-3.45")
#' AEP_T1 <- Prod_RSF(rsf87,3756122,1625482,pc)
Prod_RSF <- function(rsf, Xp, Yp, power_curve) {
    ## compute the production [MWh] for one turbine type(power curve) on one position.

    library(data.table)

    if (!is(rsf, "data.table")) {
        stop("the rsf Object pass in Prod_RSF function is not a data.table")
    } else if (Xp < rsf[, min(X)] |
               Xp > rsf[, max(X)] |
               Yp < rsf[, min(Y)] |
               Yp > rsf[, max(Y)]) {
        stop("The position Xp,Yp pass in Prod_RSF function, is outside the rsf Area")
    } else if (sum(dim(power_curve) == c(81, 2)) != 2) {
        stop("the power_curve Object pass in Prod_RSF function is not in the correct dimension (81x2)")
    }


    #node <- rsf[abs(X-Xp) < 12.5 & abs(Y-Yp) < 12.5, ]   #closest node to the turbine position (reason of 12.5 ==> rsf resolution = 25m)
    node <- rsf[which.min(abs(X - Xp) + abs(Y - Yp))]     # have to be change the rsf grid is not always constant (25m)

    ws_ms <- seq(0, 40, 0.5)

    DT <- data.table(matrix(0, 81, (length(rsf) - 9) / 3 + 1))   #intialisation de la data table

    DT[, 1 := ws_ms, with = FALSE]

    selection <- grep("F", names(node), value = TRUE)
    tot_freq <- sum(node[, selection, with = FALSE])                  # correct the cfd rsf rouding in the frequency sum ( should be equal to 1000....i dont know if it has an impact in Windpro Calc)

    #print(tot_freq)

    for (i in 1:((length(node) - 9) / 3)) {   #migh be optimized with mapply and grep into coloumn names...but it is a little table
        freq <- as.numeric(node[, 7 + i * 3, with = FALSE]) * 1000 / tot_freq
        aave <- as.numeric(node[, 8 + i * 3, with = FALSE] / 10)
        kave <- as.numeric(node[, 9 + i * 3, with = FALSE] / 100)
        DT[, i + 1 := freq * dweibull(ws_ms, kave, aave) / 2, with = FALSE]  # selection de la colonnes par son numero, on divise par deux car on est en step de 0.5 m/s
    }

    DT[, sum := rowSums(.SD), by = V1]
    sum(power_curve[, 2] * DT[, sum]) * 8760 / (1000 * 1000)     # 8760 for hours in a year, /1000 kW => MW,  /1000 freq table is in per 1000.
}

#' Bensh_RSF Function
#'
#' Benshark wind turbine production of two different rsf (wind mapping). It is mainly use to benchmark 36 versus 12 sectors rsf.
#' @param rsf36: first rsf object already loaded into memory
#' @param rsf12: second rsf object already loaded into memory
#' @param layout: data frame of two columns containing the turbine positions (X,Y).
#' @param power_curve: power curve used to compute production. the format of power_curve is a data frame with two columns (Wind speed [m/s], power [kW]). can be generate by Gen_power_curve function.
#' @keywords rsf, production, benshmark
#' @export
#' @examples
#' Compute AEP for V117 on the position (3756122,1625482):
#' layout <- read.csv("layout_file.csv")
#' pc <- Gen_power_curve("V117-3.45")
#' bench_table <- Bensh_RSF(rsf_H1_36,rsf_H1_12, layout, pc)
Bensh_RSF <- function(rsf36, rsf12, layout, power_curve) {
    if (max(rsf12[, Aave] != rsf36[, Aave])) {
        print("Warning: the input RSF files pass into Bensh_RSF does not seems coherent - Aave columns is not the same")
    } else if (min(layout[, 1]) < rsf36[, min(X)] |
               max(layout[, 1]) > rsf36[, max(X)] |
               min(layout[, 2]) < rsf36[, min(Y)] |
               max(layout[, 2]) > rsf36[, max(Y)]) {
        stop("at least one position of The layout pass in Bensh_RSF function, is outside the rsf Area")
    } else if (sum(dim(power_curve) == c(81, 2)) != 2) {
        stop("the power_curve Object pass in Bensh_RSF function is not in the correct dimension (81x2)")
    }

    names(layout) <- c("X", "Y")
    results <- data.table(layout)

    for (i in 1:nrow(layout)) {
        results[i, ':=' (
            AEP_36 = Prod_RSF(rsf36, X, Y, power_curve),
            AEP_12 = Prod_RSF(rsf12, X, Y, power_curve)
        )]
    }

    results[, ratio12_36 := AEP_12 / AEP_36]

}

#' Gen_power_curve Function
#'
#' Generate a power curve. the format of the generated power_curve is a data frame with two columns (Wind speed [m/s], power [kW])
#' @param turbine: choose the turbine type. (for the moment only "V117-3.45" is availbale)
#' @keywords rsf, production, power curve
#' @export
#' @examples
#' generate V117 power curve and store it in "pc" object:
#' pc <- Gen_power_curve("V117-3.45")
Gen_power_curve <- function(turbine = "V117-3.45") {
    # genrate a power curve data frame, that can be used in other function. to
    # match the dim of other function, the power curve need to go from 0ms to
    # 40mS by 0.5ms. (to dev: add other turbine types)
    PC <- data.frame(
        "wind_speed[m/s]" = seq(0, 40, 0.5),
        "Power[kW]" = c(
            0,
            0,
            0,
            0,
            0,
            0,
            22,
            78,
            150,
            237,
            340,
            466,
            617,
            796,
            1006,
            1247,
            1522,
            1831,
            2178,
            2544,
            2905,
            3201,
            3374,
            3435,
            3448,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            3450,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0
        )
    )
}

#' Interpol_RSF Function
#'
#' Based on the interpolation of rsf at H1 and H2, this function creates a combined rsf file for H!, H2 and the height(s) provided in the vector "layer_H". it assumes a shear profile of the wind speed, and linear interpolation for the other parameters (Freq, and k). The function can actually extrapolate if needed. It is recommended that the difference H1-H2 is at least 20m.
#' @param rsf_H1: first rsf object already loaded into memory (wind mapping at height H1)
#' @param rsf_H2: second rsf object already loaded into memory (wind mapping at height H2)
#' @param layer_H: vector containing target Heights [m] to add (height above the ground)
#' @keywords rsf, interpolation
#' @export
#' @examples
#' Generate rsf at 87m, 105m, 137m based on the rsf_100 and rsf150 inputs (respectively at 100m and 150m).
#' rsf_3D <- Interpol_RSF(rsf_100, rsf_150,c(87,105,137))
Interpol_RSF <- function(rsf_H1, rsf_H2, layer_H) {
    # create a rsf file at the heigth(s) provided in the vector "layer_H", based on the interpolation of rsf at H1 and H2. it assumes a shear profile of the wind speed.

    # check the coherence of the grids (X,Y) of the two RSF. In order to work with column operations, grids should be exactly similar.
    if (min(rsf_H1[, X] == rsf_H2[, X]) == 0 |
        min(rsf_H1[, Y] == rsf_H2[, Y]) == 0) {
        stop("the two rsf Objects pass in Interpol_RSF function do not have the same meshing/grid")
    }

    H1 <- unique(rsf_H1[, Height])
    H2 <- unique(rsf_H2[, Height])
    N_Sector <- unique(rsf_H2[, Sector])


    #interpol_table is a table with the vertial "gradient of each parameters: Freq (linear), Avae (power law), k (linear)
    interpol_table <- rsf_H1[, .(Label, X, Y, Z, Height, Aave, k, Blank, Sector)]  #initialisation of the fist columns


    for (i in 1:N_Sector) {  # fill sectorwise gradients

        j <- i * 3 + 7 # column number in the rsf data.tables for the concerned sector

        # vertical "gradient" for Freq : assumed linear
        interpol_table[, paste("Slope_F", i - 1, sep = "") := (rsf_H1[, j, with = FALSE] - rsf_H2[, j, with = FALSE]) / (H1 - H2), with = FALSE]

        # vertical "gradient" for Aave: assumed power law
        interpol_table[, paste("Shear_", i - 1, sep = "") := (log(rsf_H1[, j + 1, with = FALSE] / rsf_H2[, j + 1, with = FALSE]) / log(H1 / H2)), with = FALSE]

        # vertical "gradient" for k : assumed linear
        interpol_table[, paste("Slope_k", i - 1, sep = "") := (rsf_H1[, j + 2, with = FALSE] - rsf_H2[, j + 2, with = FALSE]) / (H1 - H2), with = FALSE]
    }


    interpol_table[, ':=' (
        Aave = log(rsf_H1[, Aave] / rsf_H2[, Aave]) / log(H1 / H2),
        k = (rsf_H1[, k] - rsf_H2[, k]) / (H1 - H2))]   #the overall gradient columns (not sectorwise)

    setnames(interpol_table, "Aave", "Shear")
    setnames(interpol_table, "k", "Slope_k")

    rsd_3D <- data.table(rsf_H1)

    for (H3 in layer_H) {
        # to merge with the code above where layer_H = one height only)

        rsf_H3 <- data.table(rsf_H1)
        rsf_H3[, ':='(
            Aave = round(rsf_H1[, Aave] * (H3 / H1) ^ interpol_table[, Shear], 2),
            k = round(rsf_H1[, k] + interpol_table[, Slope_k] * (H3 - H1), 3),
            Height = H3)]

        for (i in 1:N_Sector) {
            j <- i * 3 + 7 # column number in the rsf data.tables for the concerned sector

            rsf_H3[, j := round(rsf_H1[, j, with = FALSE] + interpol_table[, j, with = FALSE] * (H3 - H1), 0), with = FALSE]          #Freq
            rsf_H3[, j + 1 := round(rsf_H1[, j + 1, with = FALSE] * (H3 / H1) ^ interpol_table[, j + 1, with = FALSE], 0), with = FALSE]     #Aave
            rsf_H3[, j + 2 := round(rsf_H1[, j + 2, with = FALSE] + interpol_table[, j + 2, with = FALSE] * (H3 - H1), 0), with = FALSE]    #k

        }

        rsf_H3[1, 1] <-"sheared"                   #tag ad to the first line of the RSF to know that it has been interpolate by this function
        rsf_H3 <- rsf_H3[complete.cases(rsf_H3)]   #remove some artefact line (because of neg value or 0 in the input rsf)
        rsd_3D <- rbind(rsd_3D, rsf_H3)

    }

    rsd_3D <- rbind(rsd_3D, rsf_H2)

}

#' Shear2png Function
#'
#' Export the "RSF_shear" rsf to a serie of georefenced images (*.png and *.pgw). For more information about the *.pgw (word file): https://en.wikipedia.org/wiki/World_file). those images can be directly used in windpro as background maps.
#' @param RSF_shear: "shear rsf" already loaded into memory (could be generated by the function Shear_RSF())
#' @keywords rsf, shear, png
#' @export
#' @examples
#' generate a shear rsf from rsf_87 and rsf137 and export it into png pictures.
#' shear_ <- Shear_RSF(rsf87,rsf137)
#' Shear2png(shear_)
Shear2png <- function(RSF_shear) {
    #Save the shear RSF to a serie (average and for all sectors) of georefenced images (*.png and *.pgw).
    #for more information about the *.pgw (word file): https://en.wikipedia.org/wiki/World_file)


    library(ggplot2)
    library(cowplot)   # just for the use a function plot_grid() to remove the white border of the plot
    library(png)       # to insert the legend on the graph (pictures)

    mypng <- readPNG('map_legend_shear.png')

    #Parameters initialisation for creating the picture ----------------------------------------------

    # Parameters for the png
    pxl_size <- abs(unique(RSF_shear$X)[1] - unique(RSF_shear$X)[2])    #pixel size [m]
    img_size_x <- length(unique(RSF_shear$X))   #image size X [pxl]
    img_size_y <- length(unique(RSF_shear$Y))   #image size Y [pxl]

    # Parameters for the pgw file *-*-*-*-*-*-*-*-*-
    first_pxl_x <- min(RSF_shear$X)
    first_pxl_y <- max(RSF_shear$Y)


    # Creation of pictures ----------------------------------------------------------------------------

    N_Sector <- as.integer(unique(RSF_shear$Direction))

    for (s in N_Sector) {
        png_name <- paste("shear_", s, ".png", sep = "")
        pwg_name <-
            paste("shear_", s, ".pgw", sep = "")          #same name than the PNG (but with pwg extension)
        cat(
            pxl_size,
            "\n0\n0\n-",
            pxl_size,
            "\n",
            first_pxl_x,
            "\n",
            first_pxl_y,
            file = pwg_name,
            sep = ""
        ) #content in pwg is always the same whatever the sector

        shear <- subset(RSF_shear, Direction == as.integer(s))

        plot <- ggplot(shear, aes(X, Y, z = ShearAve)) +
            geom_raster(aes(fill = ShearAve)) +
            scale_fill_distiller(palette = "Spectral", limits = c(-0.1, 0.8)) +
            labs(x = NULL, y = NULL) +
            theme_nothing() +
            stat_contour(aes(colour = ..level..),
                         breaks = c(0, 0.6),
                         show.legend = NA) + guides(colour = FALSE) +
            annotation_raster(
                mypng,
                ymin = first_pxl_y - pxl_size * 99 ,  # insert legend
                ymax = first_pxl_y,
                xmin = first_pxl_x,
                xmax = first_pxl_x + pxl_size * 46
            )

        print(plot_grid(plot, scale = 1.1))   # remove those bloody borders

        dev.copy(png,
                 file = png_name,
                 width = img_size_x,
                 height = img_size_y)
        dev.off()
    }
}

#' Shear2pngV2 Function
#'
#' Save (export) the "RSF_shear" rsf to a serie (average and for all sectors) of georefenced images (*.png and *.pgw). for more information about the *.pgw (word file): https://en.wikipedia.org/wiki/World_file). those images can be directly used in windpro as background maps.
#' @param RSF_shear: "Shear rsf" already loaded into memory
#' @keywords rsf, shear, png
#' @export
#' @examples
#' generate a shear rsf from rsf_87 and rsf137 and export it into png pictures.
#' shear_ <- Shear_RSF(rsf87,rsf137)
#' Shear2pngV2(shear_)
Shear2pngV2 <- function(RSF_shear) {
    #Save the shear RSF to a serie (average and for all sectors) of georefenced images (*.png and *.pgw).
    #for more information about the *.pgw (word file): https://en.wikipedia.org/wiki/World_file)


    library(ggplot2)
    library(cowplot)   # just for the use a function plot_grid() to remove the white border of the plot
    library(png)       # to insert the legend on the graph (pictures)

    mypng <- readPNG('map_legend_shear.png')

    #Parameters initialisation for creating the picture ----------------------------------------------

    # Parameters for the png
    pxl_size <- abs(unique(RSF_shear$X)[1] - unique(RSF_shear$X)[2])    #pixel size [m]
    img_size_x <- length(unique(RSF_shear$X))   #image size X [pxl]
    img_size_y <- length(unique(RSF_shear$Y))   #image size Y [pxl]

    # Parameters for the pgw file *-*-*-*-*-*-*-*-*-
    first_pxl_x <- min(RSF_shear$X)
    first_pxl_y <- max(RSF_shear$Y)


    # Creation of pictures ----------------------------------------------------------------------------

    N_Sector <- as.integer(unique(RSF_shear$Direction))

    for (s in N_Sector) {
        png_name <- paste("shear_", s, ".png", sep = "")
        pwg_name <-
            paste("shear_", s, ".pgw", sep = "")          #same name than the PNG (but with pwg extension)
        cat(
            pxl_size,
            "\n0\n0\n-",
            pxl_size,
            "\n",
            first_pxl_x,
            "\n",
            first_pxl_y,
            file = pwg_name,
            sep = ""
        ) #content in pwg is always the same whatever the sector

        shear <- subset(RSF_shear, Direction == as.integer(s))

        plot <- ggplot(shear, aes(X, Y, z = ShearAve)) +
            geom_raster(aes(fill = ShearAve)) +
            scale_fill_distiller(palette = "Spectral", limits = c(-0.1, 0.8)) +
            labs(x = NULL, y = NULL) +
            theme_nothing() +
            stat_contour(aes(colour = ..level..),
                         breaks = c(0, 0.6),
                         show.legend = NA) + guides(colour = FALSE) +
            annotation_raster(
                mypng,
                ymin = first_pxl_y - pxl_size * 99 , # insert legend (existing picture load in memory)
                ymax = first_pxl_y,
                xmin = first_pxl_x,
                xmax = first_pxl_x + pxl_size * 46
            )

        print(plot_grid(plot, scale = 1.1))   # remove those bloody borders

        dev.copy(png,
                 file = png_name,
                 width = img_size_x,
                 height = img_size_y)
        dev.off()
    }
}

#' ShearTable Function
#'
#' Re-format the shear table (result from Shear_RSF) from (x,y,sector) table (3 x n) --> ((x,y) ~ sector) table (2+n_sector x n).it then writes the result in a csv file (can be used afterward to extrapolate to new HH or extract the shear value on some positions)
#' @param RSF_shear: "Shear rsf" already loaded into memory
#' @param output_name: name for the output file. it has to contain the extension "*.csv"
#' @keywords rsf, shear, csv
#' @export
#' @examples
#' generate a shear rsf from rsf_87_12S and rsf137_12S and export it into png pictures.
#' shear_ <- Shear_RSF(rsf87_12S,rsf137_12S)
#' ShearTable(shear_,"Sheartable.csv")
ShearTable <- function(RSF_shear, output_name = "sheartable.csv") {
    # Re-format the shear table (result from Shear_RSF) from (x,y,sector) table (3 x n) --> ((x,y) ~ sector) table (2+n_sector x n)
    # write the result in a csv file (can be used afterward to extrapolate to new HH or extract the shear value on some positions)
    # Note for later: should be directly implemented into Shear_RSF (need to change the shear2png function)

    library(tidyr)    # use for the spread function

    mat <- spread(data = RSF_shear,
               key = Direction,
               value = ShearAve)
    mat <- subset(mat, select = -c(Label, Z, k, Blank, Sector))            #remove unnecessary columns
    mat <- round(x = mat, digits = 2)
    write.csv2(x = mat,
               file = output_name,
               row.names = FALSE)
    mat

}

#' ShearExtract Function
#'
#'Extract shear at turbine locations, and create a *.shr file for VSC (use as CFD input)
#' @param shear_csv: shear file in CSV format (output of ShearTable() function)
#' @param layout_csv: layout file in csv format (X,Y)
#' @param output_name: name of the output file,it has to contain the extension "*.shr"
#' @keywords rsf, shear, extraction
#' @export
#' @examples
#' Exract directional shear for each turbine location present in the file "layout_file.csv".
#' shear_ <- Shear_RSF(rsf87,rsf137)
#' ShearExtract("Sheartable.csv","layout_file.csv","shear_VSC_input.shr")
ShearExtract <- function(shear_csv, layout_csv, output_name = "shear_VSC_input.shr") {
        # Extract shear at turbine locations, and creat a *.shr file for VSC (use as CFD input)
        # shear_csv: shear file in CSV format (output of "ShearTable" function)
        # layout_csv: layout file in csv format (X,Y)
        # output_name: name of the output file (*.shr)

        library(data.table)

        shear <- read.csv2(shear_csv)
        shear <- data.table(shear)
        shear <- subset(shear, select = -c(X999))

        layout <- read.csv(layout_csv)
        layout <- layout[complete.cases(layout), ]

        con <- file(output_name, "w")
        writeLines("Wind Shear deducted by CFD-RSF difference", con = con)
        close(con)

        layout_shear <- data.table(NULL)

        for (i in 1:nrow(layout)) {
            Xp <- layout[i, 1]
            Yp <- layout[i, 2]
            buffer <-
                shear[which.min(abs(X - Xp) + abs(Y - Yp))]
            layout_shear <- rbind(layout_shear, buffer)
        }

        layout_shear[, 1] <- layout[, 1]  # change the X to the exact X of the turbine
        layout_shear[, 2] <- layout[, 2]  # change the Y to the exact Y of the turbine

        write.table(
            layout_shear,
            file = output_name,
            append = TRUE,
            row.names = FALSE,
            sep = ";"
        )

    }

#' Vizu_RSF Function
#'
#' Interactive plot of a rsf object (direction, Height, parameters...)
#' @param rsf: first rsf object already loaded into memory (wind mapping)
#' @keywords rsf, plot, vizualisation
#' @export
#' @examples
#' Vizu_RSF(rsf_87)
Vizu_RSF <- function(rsf) {
    # dynamic plot(raster map) of rsf variables (A, k, freq, shear).

    library(ggplot2)
    library(manipulate)    # for dynamic plot (you choose the variable you want to plot)

    names_col <- sort(names(rsf))    # creation of the list of variable availabe in the plot
    names_col <- as.list(names_col)
    picker_list <- names_col[names_col != "Label"]  # remove usless variables for the plot
    picker_list <- picker_list[picker_list != "X"]
    picker_list <- picker_list[picker_list != "Y"]
    picker_list <- picker_list[picker_list != "Z"]
    picker_list <- picker_list[picker_list != "Height"]
    picker_list <- picker_list[picker_list != "Blank"]
    picker_list <- picker_list[picker_list != "Sector"]

    p <- ggplot(rsf, aes(X, Y))

    manipulate(
        p + geom_raster(aes(fill = subset(rsf, select = variable)[[1]])) +  #[[1]] --> needed because of input format condition of geom_raster
            scale_fill_distiller(palette = "Spectral") +
            labs(x = NULL, y = NULL),
        variable = picker(picker_list)
    )
}
