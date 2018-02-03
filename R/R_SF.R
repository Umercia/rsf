

#' RSF_convert Function
#'
#' Do basic transformation on *.rsf file (rsf = wind ReSource File): Crop, convert in 12 sector, create new hight levels, compute and export shear table.
#' @param rsf_file1: Input *.rsf file name (name should be surrounded by "").
#' @param output_name: Output *.rsf file name " (name should be surrounded by "").
#' @param crop: Should the rsf be cropped ?  TRUE/FALSE. Default value = FALSE.
#' @param layout_csv: Input *.csv file name (name should be surrounded by "")  containing the WTG positions (X,Y).Default value = NULL. This input is used to crop and benchmark the production (12 Sectors versus 36 sectors).
#' @param crop_buffer: Integer value: Distance in meters around the WTG positions that should be cropped (default value = 700m)
#' @param twelve_S: Should the rsf be converted into 12 sectors ?  TRUE/FALSE. Default value = FALSE.
#' @param rsf_file2: Input *.rsf file name (name should be surrounded by ""). Second rsf file used for shear and 3 dimensional rsf creation. Default value = NULL.
#' @param shear_out: Should the shear should be extracted and exported (csv file) ?  TRUE/FALSE. Default value = FALSE.
#' @param three_D: Should a 3 dimensional rsf created  ? TRUE/FALSE. Default value = FALSE.
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

        if(crop == TRUE & is.null(layout_csv)) {
            stop("to crop the rsf file you need to provide a layout file")
        } else if(twelve_S == TRUE & is.null(layout_csv)) {
            warning("you did not provide any layout file, the benchmark will be bypassed")
        } else if( (shear_out == TRUE | three_D == TRUE) & is.null(rsf_file2)) {
            stop("you need to provide a second rsf file in order to create a 3 Dimensional rsf or a shear matrix")
        }else if( !is.null(layer_H) & three_D == FALSE) {
            stop("To extrapolate a signal to a new specific height you need to pass three_D parameter as TRUE")
        }


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

            if (is.null(layout_csv)){
                stop("Layout (*.csv) is needed to crop the map around it.")
                }

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

            if (!is.null(layout_csv)) {

                bench_table <- Bensh_RSF(rsf_H1, rsf_H1_12, layout, Gen_power_curve())

                Ratio <- bench_table[, mean(ratio12_36)]
                Ratio <- round(Ratio * 100 - 100, 1) # extract the last unit and digit.

                write.csv(bench_table,
                          paste("benchmark_table_", output_name, ".csv", sep = ""))
                message(paste("benchmark_table_", output_name, ".csv created in your working directory: ", getwd(), sep = ""))

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
                message(paste("benchmark_graph_", output_name, ".pdf created in your working directory: ", getwd(), sep = ""))

            } else {

                warning("!!! No 36S/12S benchmark performed because no layout (*.csv) provided.")
                Ratio <- "XX"

            }

            rsf_H1 <- rsf_H1_12
            rm(rsf_H1_12)

            if (rsf2 == TRUE) {
                rsf_H2_12 <- S36_to_S12_RSF(rsf_H2)
                rsf_H2 <- rsf_H2_12
                rm(rsf_H2_12)
            }

            output_name <-
                paste("[12S_",
                      Ratio,
                      "]",
                      output_name,
                      sep = "")

        }


        # Shear table -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        if (shear_out == TRUE) {
            shear_table_file_name <- paste("[shear_table]", output_name, ".csv", sep = "")
            shear_table <- Shear_RSF(rsf_H2, rsf_H1)
            ShearTable(shear_table, shear_table_file_name)  #create shear file
            message(paste(shear_table_file_name, " created in your working directory: ", getwd(), sep = ""))
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

            # creation of the 3 dimensional rsf
            rsf_H1 <- Interpol_RSF(rsf_H1, rsf_H2, layer_H)
            output_name <- paste("[3D]", output_name, sep = "")
        }



        # Name of the output file: tag the heights of the rsf3D -------------------------------

        heights_n <- paste(unique(rsf_H1$Height)[order(unique(rsf_H1$Height))])
        heights_n <- paste("[" , heights_n, "]",sep = "")

        tag <- as.character()

        for (i in 1:length(heights_n)){

            tag <- paste(tag,heights_n[i],sep = "")

        }

        output_name <- gsub(pattern = ".rsf|.RSF",replacement = "",x = output_name)  # remove the .rsf extension

        # Write to file -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

        Write_RSF(rsf_H1 , output_name)
        message(paste(output_name, " created in your working directory: ", getwd(), sep = ""))

    }


#' RSF_add_layer
#'
#' This function allows to add some layers (height) to a 3 dimensions rsf (meaning a rsf file with already two heights).
#' The new rsf file will be created in the same folder, with the extra height(s) tagged in the file name.
#' Note that the maximum and minimum height levels are used to interpolate or extrapolate. In case of extrapolation, at least 20m difference is recommended in-between min and max level height
#' @param rsf3D_file: 3 dimensional"*.rsf" file name (file name should be surrounded by "")
#' @param layer_H: vector containing  the level height to add.
#' @param output_name: Output *.rsf file name " (name should be surrounded by "").
#' @keywords rsf, layer, 3D
#' @export
#' @examples
#' Add height 30 and 217 to the rsf file "[3D][12S_-0.1]RSF-CFD_Douglas_West.rsf" :
#' RSF_add_layer(rsf3D_file = "[3D][12S_-0.1]RSF-CFD_Douglas_West.rsf",layer_H = c(30,217))
RSF_add_layer <- function(rsf3D_file, layer_H, output_name = rsf3D_file){

    # read inputs -------------------------------------------------------------------------
    rsf3D <- Read_RSF(rsf3D_file)
    maxHeight <- max(unique(rsf3D$Height))
    minHeight <- min(unique(rsf3D$Height))

    if (abs(maxHeight - minHeight) < 20) {
        warning(paste("the heigth difference is below 20m", maxHeight, minHeight))
    }

    # creation of new levels --------------------------------------------------------------
    rsf_new_height <- Interpol_RSF(rsf_H1 = rsf3D[Height == maxHeight,],
                                   rsf_H2 = rsf3D[Height == minHeight,],
                                   layer_H = layer_H)

    rsf_new_height <- rsf_new_height[Height %in% layer_H,]

    rsf3D <- rbind(rsf3D,rsf_new_height)

    rm(rsf_new_height)

    # Name of the output file: tag the heights of the rsf3D -------------------------------

    heights_n <- paste(unique(rsf3D$Height)[order(unique(rsf3D$Height))])

    heights_n <- paste("[" , heights_n, "]",sep = "")

    tag <- as.character()

    for (i in 1:length(heights_n)){

        tag <- paste(tag,heights_n[i],sep = "")

    }
    output_name <- gsub(pattern = ".rsf|.RSF",replacement = "",x = output_name) # remove the .rsf extension

    output_name <- paste(output_name,tag,sep = "")


    # Write the results into a new rsf file ------------------------------------------------

    Write_RSF(rsf3D, output_file_name = output_name)
}

#' Read_RSF Function
#'
#' This function allows you to read a *.rsf file. it recognises if it is 36 or 12 sectors, load it accordingly, with columns names.
#' rsf = wind ReSource File. it is a wind mapping output from wind flow model. the new object will be a data.table
#' @param Input_file: "*.rsf" file ()
#' @keywords rsf, read
#' @export
#' @examples
#' load into memory an rsf file into the object rsf_166:
#' rsf_166 <- Read_RSF("RSF-windresource-CFD_Aldermyrberget 166m.rsf")
Read_RSF <- function(Input_file) {
        ## function read an *.rsf file. it recognises if it is 36 or 12 sectors, load it accordingly, with columns names
        library("readr")              ## use for the read_fwf function much more faster than read.fwf
        library("data.table")         ## use forthe output format

        file_type = substr(x = Input_file, start = nchar(Input_file)-2,stop = nchar(Input_file))

        if (file_type == "rsf" | file_type == "wrg") {
                print(paste("Read_RSF, detected input file format:", file_type))
        } else {
                stop(paste("Error Read_RSF - Unknow file format:", file_type))
        }

        con <- file(Input_file, "r")
        line <- readLines(con, 1)      ## sample line to check if it is a 36 sectors or 12 sectors file
        line <- readLines(con, 1)      ## read the second line: skip the first line in case of a *.wrg input
        close(con)

        if (nchar(line) > 400) {
                ## it means it should be a 36 sectors (actually 36S should be exactly 540)
                RSF_Col_Format <- c(10, 10, 10, 8, 5, 5, 6, 15, 3, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4,
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
                RSF_Col_Name <- c(RSF_Col_Name, paste(c("F36_", "A36_", "k36_"), rep(c("000","010","020","030","040","050","060","070","080","090","100","110","120","130","140","150","160","170","180","190","200","210","220","230","240","250","260","270","280","290","300","310","320","330","340","350"), each = 3), sep = ""))

        } else{
                ##if less than 400, then it means it should be a 12 sectors (actually 12S should be exactly 228)
                RSF_Col_Format <- c(10, 10, 10, 8, 5, 5, 6, 15, 3, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4,
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
                RSF_Col_Name <- c(RSF_Col_Name, paste(c("F12_", "A12_", "k12_"), rep(c("000","030","060","090","120","150","180","210","240","270","300","330"), each = 3), sep = ""))
        }



        if (file_type == "rsf") {
                RSF_table <- data.table(read_fwf(
                        Input_file,
                        fwf_widths(RSF_Col_Format, col_names = RSF_Col_Name)
                ))
        } else if (file_type == "wrg") {
                RSF_table <- data.table(read.table(file = Input_file, skip = 1, col.names = RSF_Col_Name))

        }


        if (is.na(RSF_table$Label[1])) {RSF_table$Label <- "RRR"} # to avoid issue in other functions

        RSF_table <- RSF_table[order(RSF_table[, Height],RSF_table[, Y], RSF_table[, X])]   #Order the file: Height, Y and X. Windpro is not always coherent when creating an rsf file from WasP.

        for(col in grep(pattern = "A[0-9]",x = names(RSF_table), value = TRUE)){    # set 2 m/s as the minimum A ("0m/s" creates issue in the extrapolation)
                set(RSF_table, i=which(RSF_table[[col]] < 20), j=col, value= 20)

        }

        RSF_table

}

#' Crop_RSF Function
#'
#' Crop an *.rsf file around a rectangle define by Xmax, Xmin, Ymax, Ymin.
#' @param rsf: rsf object already loaded into memory
#' @param Xmax: X maximum limit
#' @param Xmin: X minimum limit
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
    RSF12[, F12_000 := F36_350 + F36_000 + F36_010]
    RSF12[, A12_000 := round((A36_350 * F36_350 + A36_000 * F36_000 + A36_010 * F36_010) /
                               F12_000, 0)]
    RSF12[, k12_000 := round((k36_350 * F36_350 + k36_000 * F36_000 + k36_010 * F36_010) /
                               F12_000, 0)]

    RSF12[, F12_030 := F36_020 + F36_030 + F36_040]
    RSF12[, A12_030 := round((A36_020 * F36_020 + A36_030 * F36_030 + A36_040 *
                                 F36_040) / F12_030, 0)]
    RSF12[, k12_030 := round((k36_020 * F36_020 + k36_030 * F36_030 + k36_040 *
                                 F36_040) / F12_030, 0)]

    RSF12[, F12_060 := F36_050 + F36_060 + F36_070]
    RSF12[, A12_060 := round((A36_050 * F36_050 + A36_060 * F36_060 + A36_070 *
                                 F36_070) / F12_060, 0)]
    RSF12[, k12_060 := round((k36_050 * F36_050 + k36_060 * F36_060 + k36_070 *
                                 F36_070) / F12_060, 0)]

    RSF12[, F12_090 := F36_080 + F36_090 + F36_100]
    RSF12[, A12_090 := round((A36_080 * F36_080 + A36_090 * F36_090 + A36_100 *
                                 F36_100) / F12_090, 0)]
    RSF12[, k12_090 := round((k36_080 * F36_080 + k36_090 * F36_090 + k36_100 *
                                 F36_100) / F12_090, 0)]

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
    RSF12[, paste(c("F36_", "A36_", "k36_"), rep(c("000","010","020","030","040","050","060","070","080","090","100","110","120","130","140","150","160","170","180","190","200","210","220","230","240","250","260","270","280","290","300","310","320","330","340","350"), each = 3), sep = "") := NULL]

}

#' Write_RSF Function
#'
#' Write a rsf object from memory to a file (*.rsf). This file can be used with Windpro.
#' @param RSF: rsf object already loaded into memory
#' @param output_file_name: name for the output file. it has to contain the extension .rsf.
#' @keywords rsf, write
#' @export
#' @examples
#' write from memory a rsf (data table) into the file "rsf_80m_Aldermyberget.rsf":
#' Write_RSF(rsf80, "rsf_80m_Aldermyberget.rsf")
Write_RSF <- function(RSF, output_file_name) {
    # write an rsf table ("RSF") in memory to a *.rsf file ("output_file_name").

    library("gdata") ## use for the write.fwf function (write.fwf writes object in *f*ixed *w*idth *f*ormat )


    output_file_name <- paste(output_file_name, ".rsf", sep= "")  #add rsf extension to the name

    Num_Sectors <- RSF[, unique(Sector)]
    RSF_Col_Format <- c(c(10, 10, 10, 8, 5, 5, 6, 15, 3), rep(c(4, 4, 5), Num_Sectors))


    for(col in grep(pattern = "A[0-9]",x = names(RSF), value = TRUE)){    # set 2 m/s as the minimum A (remove extrapolation artefacts)
            set(RSF, i=which(RSF[[col]] < 20), j=col, value= 20)

    }

    for(col in grep(pattern = "A[0-9]",x = names(RSF), value = TRUE)){    # set 20 m/s as the maximum A (remove extrapolation artefacts)
            set(RSF, i=which(RSF[[col]] > 200), j=col, value= 200)

    }



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
#' Benchmark wind turbine production of two different rsf (wind mapping). It is mainly used to benchmark 36 versus 12 sectors rsf.
#' @param rsf36: first rsf object already loaded into memory
#' @param rsf12: second rsf object already loaded into memory
#' @param layout: data frame of two columns containing the turbine positions (X,Y).
#' @param power_curve: power curve used to compute production. the format of power_curve is a data frame with two columns (Wind speed [m/s], power [kW]). can be generate by Gen_power_curve function.
#' @keywords rsf, production, benchmark
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
#' @param turbine: choose the turbine type. (for the moment only "V117-3.45" is available)
#' @keywords rsf, production, power curve
#' @export
#' @examples
#' generate V117 power curve and store it in "pc" object:
#' pc <- Gen_power_curve("V117-3.45")
Gen_power_curve <- function(turbine = "V117-3.45") {
    # Genrate a power curve data frame, that can be used in other functions. to
    # match the dim of other functions, the power curve need to go from 0 m/s to
    # 40 m/s by 0.5m/s. (to dev: add other turbine types)
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
#' @param layer_H: vector containing target heights [m] to add (height above the ground)
#' @keywords rsf, interpolation
#' @export
#' @examples
#' Generate rsf at 87m, 105m, 137m based on the rsf_100 and rsf150 inputs (respectively at 100m and 150m).
#' rsf_3D <- Interpol_RSF(rsf_100, rsf_150,c(87,105,137))
Interpol_RSF <- function(rsf_H1, rsf_H2, layer_H) {
    # Create a rsf file at the height(s) provided in the vector "layer_H", based on the interpolation of rsf at H1 and H2. it assumes a shear profile of the wind speed.

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
#' Export the "RSF_shear" rsf to a series of georefenced images (*.png and *.pgw). For more information about the *.pgw (word file): https://en.wikipedia.org/wiki/World_file). Those images can be directly used in windpro as background maps.
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
#' Save (export) the "RSF_shear" rsf to a series (average and for all sectors) of georefenced images (*.png and *.pgw). for more information about the *.pgw (word file): https://en.wikipedia.org/wiki/World_file). Those images can be directly used in windpro as background maps.
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
#' Re-format the shear table (result from Shear_RSF) from (x,y,sector) table (3 x n) --> ((x,y) ~ sector) table (2+n_sector x n).it then writes the result in a csv file (can be used afterwards to extrapolate to new HH or extract the shear value on some positions)
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
    # write the result in a csv file (can be used afterwards to extrapolate to new HH or extract the shear value on some positions)
    # Notes for later: should be directly implemented into Shear_RSF (need to change the shear2png function)

    library(tidyr)    # use for the spread function

    mat <- spread(data = RSF_shear,
               key = Direction,
               value = ShearAve)
    mat <- subset(mat, select = -c(Label, k, Blank, Sector))            #remove unnecessary columns
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
        # Extract shear at turbine locations, and create a *.shr file for VSC (use as CFD input)
        # shear_csv: shear file in CSV format (output of "ShearTable" function)
        # layout_csv: layout file in csv format (X,Y)
        # output_name: name of the output file (*.shr)

        library(data.table)

        shear <- read.csv2(shear_csv)
        shear <- data.table(shear)
        shear <- subset(shear, select = -c(X999, Z))

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

#' RSF_plot Function
#'
#' Interactive plot (raster/map) of an rsf object ( direction versus (Height, A, k, Wind speed) )
#' @param rsf_file: Input *.rsf file name (name should be surrounded by "").
#' @param layout_csv: Input *.csv file name (name should be surrounded by "") containing the WTG positions (X,Y).Default value = NULL. if provided, turbine positions would be plotted on the map.
#' @keywords rsf, plot, vizualisation
#' @export
#' @examples
#' RSF_plot(rsf_file = "Aldermyrberget 166m.rsf", layout_csv = "layout V14A.csv")
RSF_plot <- function(rsf_file,layout_csv = NULL){

    library(ggplot2)
    library(manipulate)


# Reading Inputs *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

    rsf <- Read_RSF(rsf_file)

    Num_Sectors <- rsf[, unique(Sector)]

    if (!is.null(layout_csv)) {
        layout <- read.csv(layout_csv)
        layout <- layout[complete.cases(layout), ]
        names(layout) <- c("X", "Y")
    }


# unit conversion (k*100, A*10, average wind speed compute) *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    if(Num_Sectors == 12) {

        rsf[, ':=' (k12_000 = k12_000/100,
                    k12_030 = k12_030/100,
                    k12_060 = k12_060/100,
                    k12_090 = k12_090/100,
                    k12_120 = k12_120/100,
                    k12_150 = k12_150/100,
                    k12_180 = k12_180/100,
                    k12_210 = k12_210/100,
                    k12_240 = k12_240/100,
                    k12_270 = k12_270/100,
                    k12_300 = k12_300/100,
                    k12_330 = k12_330/100)]

        rsf[, ':='(A12_000 = A12_000/10,
                   A12_030 = A12_030/10,
                   A12_060 = A12_060/10,
                   A12_090 = A12_090/10,
                   A12_120 = A12_120/10,
                   A12_150 = A12_150/10,
                   A12_180 = A12_180/10,
                   A12_210 = A12_210/10,
                   A12_240 = A12_240/10,
                   A12_270 = A12_270/10,
                   A12_300 = A12_300/10,
                   A12_330 = A12_330/10)]

        rsf[, ':='(Vave = Aave*gamma(1+1/k),
                   V12_000 = A12_000*gamma(1+1/k12_000),
                   V12_030 = A12_030*gamma(1+1/k12_030),
                   V12_060 = A12_060*gamma(1+1/k12_060),
                   V12_090 = A12_090*gamma(1+1/k12_090),
                   V12_120 = A12_120*gamma(1+1/k12_120),
                   V12_150 = A12_150*gamma(1+1/k12_150),
                   V12_180 = A12_180*gamma(1+1/k12_180),
                   V12_210 = A12_210*gamma(1+1/k12_210),
                   V12_240 = A12_240*gamma(1+1/k12_240),
                   V12_270 = A12_270*gamma(1+1/k12_270),
                   V12_300 = A12_300*gamma(1+1/k12_300),
                   V12_330 = A12_330*gamma(1+1/k12_330))]
    }

    if(Num_Sectors == 36) {

        rsf[, ':=' (k36_000 = k36_000/100,
                    k36_010 = k36_010/100,
                    k36_020 = k36_020/100,
                    k36_030 = k36_030/100,
                    k36_040 = k36_040/100,
                    k36_050 = k36_050/100,
                    k36_060 = k36_060/100,
                    k36_070 = k36_070/100,
                    k36_080 = k36_080/100,
                    k36_090 = k36_090/100,
                    k36_100 = k36_100/100,
                    k36_110 = k36_110/100,
                    k36_120 = k36_120/100,
                    k36_130 = k36_130/100,
                    k36_140 = k36_140/100,
                    k36_150 = k36_150/100,
                    k36_160 = k36_160/100,
                    k36_170 = k36_170/100,
                    k36_180 = k36_180/100,
                    k36_190 = k36_190/100,
                    k36_200 = k36_200/100,
                    k36_210 = k36_210/100,
                    k36_220 = k36_220/100,
                    k36_230 = k36_230/100,
                    k36_240 = k36_240/100,
                    k36_250 = k36_250/100,
                    k36_260 = k36_260/100,
                    k36_270 = k36_270/100,
                    k36_280 = k36_280/100,
                    k36_290 = k36_290/100,
                    k36_300 = k36_300/100,
                    k36_310 = k36_310/100,
                    k36_320 = k36_320/100,
                    k36_330 = k36_330/100,
                    k36_340 = k36_340/100,
                    k36_350 = k36_350/100)]

        rsf[, ':=' (A36_000 = A36_000/10,
                    A36_010 = A36_010/10,
                    A36_020 = A36_020/10,
                    A36_030 = A36_030/10,
                    A36_040 = A36_040/10,
                    A36_050 = A36_050/10,
                    A36_060 = A36_060/10,
                    A36_070 = A36_070/10,
                    A36_080 = A36_080/10,
                    A36_090 = A36_090/10,
                    A36_100 = A36_100/10,
                    A36_110 = A36_110/10,
                    A36_120 = A36_120/10,
                    A36_130 = A36_130/10,
                    A36_140 = A36_140/10,
                    A36_150 = A36_150/10,
                    A36_160 = A36_160/10,
                    A36_170 = A36_170/10,
                    A36_180 = A36_180/10,
                    A36_190 = A36_190/10,
                    A36_200 = A36_200/10,
                    A36_210 = A36_210/10,
                    A36_220 = A36_220/10,
                    A36_230 = A36_230/10,
                    A36_240 = A36_240/10,
                    A36_250 = A36_250/10,
                    A36_260 = A36_260/10,
                    A36_270 = A36_270/10,
                    A36_280 = A36_280/10,
                    A36_290 = A36_290/10,
                    A36_300 = A36_300/10,
                    A36_310 = A36_310/10,
                    A36_320 = A36_320/10,
                    A36_330 = A36_330/10,
                    A36_340 = A36_340/10,
                    A36_350 = A36_350/10)]

        rsf[, ':=' (V36_000 = A36_000*gamma(1+1/k36_000),
                    V36_010 = A36_010*gamma(1+1/k36_010),
                    V36_020 = A36_020*gamma(1+1/k36_020),
                    V36_030 = A36_030*gamma(1+1/k36_030),
                    V36_040 = A36_040*gamma(1+1/k36_040),
                    V36_050 = A36_050*gamma(1+1/k36_050),
                    V36_060 = A36_060*gamma(1+1/k36_060),
                    V36_070 = A36_070*gamma(1+1/k36_070),
                    V36_080 = A36_080*gamma(1+1/k36_080),
                    V36_090 = A36_090*gamma(1+1/k36_090),
                    V36_100 = A36_100*gamma(1+1/k36_100),
                    V36_110 = A36_110*gamma(1+1/k36_110),
                    V36_120 = A36_120*gamma(1+1/k36_120),
                    V36_130 = A36_130*gamma(1+1/k36_130),
                    V36_140 = A36_140*gamma(1+1/k36_140),
                    V36_150 = A36_150*gamma(1+1/k36_150),
                    V36_160 = A36_160*gamma(1+1/k36_160),
                    V36_170 = A36_170*gamma(1+1/k36_170),
                    V36_180 = A36_180*gamma(1+1/k36_180),
                    V36_190 = A36_190*gamma(1+1/k36_190),
                    V36_200 = A36_200*gamma(1+1/k36_200),
                    V36_210 = A36_210*gamma(1+1/k36_210),
                    V36_220 = A36_220*gamma(1+1/k36_220),
                    V36_230 = A36_230*gamma(1+1/k36_230),
                    V36_240 = A36_240*gamma(1+1/k36_240),
                    V36_250 = A36_250*gamma(1+1/k36_250),
                    V36_260 = A36_260*gamma(1+1/k36_260),
                    V36_270 = A36_270*gamma(1+1/k36_270),
                    V36_280 = A36_280*gamma(1+1/k36_280),
                    V36_290 = A36_290*gamma(1+1/k36_290),
                    V36_300 = A36_300*gamma(1+1/k36_300),
                    V36_310 = A36_310*gamma(1+1/k36_310),
                    V36_320 = A36_320*gamma(1+1/k36_320),
                    V36_330 = A36_330*gamma(1+1/k36_330),
                    V36_340 = A36_340*gamma(1+1/k36_340),
                    V36_350 = A36_350*gamma(1+1/k36_350),
                    Vave = Aave*gamma(1+1/k))]
    }


# picker list for interactive plot *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    names_col <- sort(names(rsf))
    picker_list <- names_col[names_col != "Label"]
    picker_list <- picker_list[picker_list != "X"]
    picker_list <- picker_list[picker_list != "Y"]
    #picker_list <- picker_list[picker_list != "Z"]
    picker_list <- picker_list[picker_list != "Height"]
    picker_list <- picker_list[picker_list != "Blank"]
    picker_list <- picker_list[picker_list != "Sector"]
    picker_list <- picker_list[picker_list != "Vave"]
    picker_list <- picker_list[picker_list != "Aave"]
    picker_list <- picker_list[picker_list != "k"]
    picker_list <- c("Vave","Aave","k", picker_list)

    picker_list <- as.list(picker_list)
    picker_list2 <- as.list(unique(rsf$Height))


# Interactive plot using manipulate *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    manipulate(
        factor= picker(picker_list),
        level = picker(picker_list2),

        if (!is.null(layout_csv)) {

            p <- ggplot(subset(rsf, Height == level),aes(X,Y, z = Z)) +
                geom_raster(aes(fill= subset(rsf,Height == level, select = factor)[[1]]))+  #[[1]] --> i guess because of input format condition of geom_raster
                stat_contour(aes(colour = ..level..), binwidth = 10) +
                geom_point(data = layout, aes(x = X, y = Y, z = NULL), shape = 1, size = 2.5, color = "black")+
                guides(colour = FALSE)+
                labs(x = NULL, y = NULL, fill = factor)


            if(substr(factor,1,4) == "Vave") {   #fix scale for comparaison, limits depend on the parameters.

                print("Processing. Please wait.")
                p + scale_fill_distiller(palette = "Spectral", limits = c(5.5, 10.5))

            }else if (substr(factor,1,1) == "V") {

                print("Processing. Please wait.")
                p + scale_fill_distiller(palette = "Spectral")

            }else if (substr(factor,1,1) == "k") {

                print("Processing. Please wait.")
                p + scale_fill_distiller(palette = "PRGn")

            }else if (substr(factor,1,1) == "A") {

                print("Processing. Please wait.")
                p + scale_fill_distiller(palette = "RdYlBu")

            } else {

                print("Processing. Please wait.")
                p + scale_fill_gradientn(colours = terrain.colors(10))
                # scale_fill_distiller(palette = palette(terrain.colors(12)))
            }

        } else {

            p <- ggplot(subset(rsf, Height == level),aes(X,Y, z = Z)) +
                geom_raster(aes(fill= subset(rsf,Height == level, select = factor)[[1]]))+  #[[1]] --> i guess because of input format condition of geom_raster
                stat_contour(aes(colour = ..level..), binwidth = 10) +
                # geom_point(data = layout, aes(x = X, y = Y, z = NULL), shape = 1, size = 2.5, color = "black")+ #removed if no layout input
                guides(colour = FALSE)+
                labs(x = NULL, y = NULL, fill = factor)


            if(substr(factor,1,4) == "Vave") {   #fix scale for comparaison, limits depend on the parameters.

                print("Processing. Please wait.")
                p + scale_fill_distiller(palette = "Spectral", limits = c(5.5, 10.5))

            }else if (substr(factor,1,1) == "V") {

                print("Processing. Please wait.")
                p + scale_fill_distiller(palette = "Spectral")

            }else if (substr(factor,1,1) == "k") {

                print("Processing. Please wait.")
                p + scale_fill_distiller(palette = "PRGn")

            }else if (substr(factor,1,1) == "A") {

                print("Processing. Please wait.")
                p + scale_fill_distiller(palette = "RdYlBu")

            } else {

                print("Processing. Please wait.")
                p + scale_fill_gradientn(colours = terrain.colors(10))


            }
        }
    )
}





#' Shear_plot Function
#'
#' Interactive plot (raster/map) of a shear table file (*.csv) (X,Y,direction versus shear)
#' @param rsf_file: Input *.rsf file name (name should be surrounded by "").
#' @param layout_csv: Input *.csv file name (name should be surrounded by "") containing the WTG positions (X,Y).Default value = NULL. if provided, turbine positions would be plotted on the map.
#' @keywords shear, plot, vizualisation
#' @export
#' @examples
#' Shear_plot("Aldermyrberget 166m.rsf",layout_csv = "layout V14A.csv")
Shear_plot <- function(shear_file,layout_csv = NULL){

    library(ggplot2)
    library(manipulate)


# Reading Inputs *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

    shear <- read.csv2(shear_file)

    if (length(shear) == 16) {     #12 sectors format
        names(shear)[-c(1,2,3)] <- paste(c("Shear_"), c("000","030","060","090","120","150","180","210","240","270","300","330","Ave"), sep = "")
    }else if (length(shear) == 40) { #36 sectors format
        names(shear)[-c(1,2,3)] <- paste(c("Shear_"), c("000","010","020","030","040","050","060","070","080","090","100","110","120","130","140","150","160","170","180","190","200","210","220","230","240","250","260","270","280","290","300","310","320","330","340","350","Ave"), sep = "")
    }

    if (!is.null(layout_csv)) {
        layout <- read.csv(layout_csv)
        layout <- layout[complete.cases(layout), ]
        names(layout) <- c("X", "Y")
    }


# picker list for interactive plot *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    names_col <- sort(names(shear))
    picker_list <- names_col[names_col != "X"]
    picker_list <- picker_list[picker_list != "Y"]
    picker_list <- picker_list[picker_list != "Shear_Ave"]
    picker_list <- c("Shear_Ave", picker_list)   # just put Shear_Ave in the first row
    picker_list <- as.list(picker_list)


# Interactive plot using manipulate *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    manipulate(
        factor= picker(picker_list),

        if (!is.null(layout_csv)) {
            print("Processing. Please wait.")
            ggplot(shear,aes(X,Y, z = Z)) +
                geom_raster(aes(fill= subset(shear, select = factor)[[1]]))+  #[[1]] --> i guess because of input format condition of geom_raster
                scale_fill_distiller(palette = "Spectral", limits = c(-0.1, 0.7))+
                stat_contour(aes(colour = ..level..), binwidth = 10)+
                geom_point(data = layout, aes(x = X, y = Y, z = NULL), shape = 1, size = 2.5, color = "black") +
                guides(colour = FALSE)+
                labs(x = NULL, y = NULL, fill = factor)


        } else {
            print("Processing. Please wait.")
            ggplot(shear,aes(X,Y, z = Z)) +
                geom_raster(aes(fill= subset(shear, select = factor)[[1]]))+  #[[1]] --> i guess because of input format condition of geom_raster
                scale_fill_distiller(palette = "Spectral", limits = c(-0.1, 0.7))+
                stat_contour(aes(colour = ..level..), binwidth = 10)+
                guides(colour = FALSE)+
                labs(x = NULL, y = NULL, fill = factor)

        }
    )
}

Extrapol_RSF <- function(rsf_H1, H2, shear = 0.2, shearmap) {
    #interpol a map from H1 to H2 using a shear value (or shear map if provided)

    N_Sector <- unique(rsf_H1[, Sector])
    H1 <- unique(rsf_H1[, Height])

    rsf_H2 <- data.table(rsf_H1)
    rsf_H2[, ':='(
        Aave = round(rsf_H1[, Aave]  * (H2 / H1) ^ shear,2),
        Height = H2)]


    for(i in 1:N_Sector){
        j <- i * 3 + 7 # column number in the rsf data.tables for the concerned sector


        rsf_H2[, j + 1 := round(rsf_H1[, j + 1, with = FALSE] * (H2 / H1) ^ shear), with = FALSE]     #Aave

    }

    rsf_H2[1, 1] <-"Extrapol"
    return(rsf_H2)

}


#' RSF_scale
#'
#' This function allows to add scale a rsf with a factor x.
#' The new rsf file will be created in the same folder, with the scale factor tagged in the file name.
#' Note that the scale is purely apply to the A parameter of the weibull distribution (k stay unchanged)
#' @param rsf_file: Input *.rsf file name (name should be surrounded by "").
#' @param scale_factor: the A value of the rsf will be multiply by this value.
#' @param output_name: Output *.rsf file name " (name should be surrounded by "").
#' @keywords rsf, scale
#' @export
#' @examples
#' Scale the rsf file "RSF-CFD_Douglas_West.rsf" by 1.11 (add 11 percents to the wind speed)  :
#' RSF_scale(rsf_file_file = "RSF-CFD_Douglas_West.rsf",scale = 1.11)
RSF_scale <- function(rsf_file, scale_factor, output_name = rsf_file){

    # read inputs -------------------------------------------------------------------------
    rsf <- Read_RSF(rsf_file)
    Num_Sectors <- rsf[, unique(Sector)]


    # creation of new levels --------------------------------------------------------------

    # scaling the A columns*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


    if(Num_Sectors == 12) {


        rsf[, ':='(Aave = round(Aave*scale_factor,2),
                   A12_000 = as.integer(round( A12_000*scale_factor)),
                   A12_030 = as.integer(round( A12_030*scale_factor)),
                   A12_060 = as.integer(round( A12_060*scale_factor)),
                   A12_090 = as.integer(round( A12_090*scale_factor)),
                   A12_120 = as.integer(round( A12_120*scale_factor)),
                   A12_150 = as.integer(round( A12_150*scale_factor)),
                   A12_180 = as.integer(round( A12_180*scale_factor)),
                   A12_210 = as.integer(round( A12_210*scale_factor)),
                   A12_240 = as.integer(round( A12_240*scale_factor)),
                   A12_270 = as.integer(round( A12_270*scale_factor)),
                   A12_300 = as.integer(round( A12_300*scale_factor)),
                   A12_330 = as.integer(round( A12_330*scale_factor)))]

    }

    if(Num_Sectors == 36) {


        rsf[, ':=' (Aave = round( Aave*scale_factor,2),
                    A36_000 = as.integer(round( A36_000*scale_factor)),
                    A36_010 = as.integer(round( A36_010*scale_factor)),
                    A36_020 = as.integer(round( A36_020*scale_factor)),
                    A36_030 = as.integer(round( A36_030*scale_factor)),
                    A36_040 = as.integer(round( A36_040*scale_factor)),
                    A36_050 = as.integer(round( A36_050*scale_factor)),
                    A36_060 = as.integer(round( A36_060*scale_factor)),
                    A36_070 = as.integer(round( A36_070*scale_factor)),
                    A36_080 = as.integer(round( A36_080*scale_factor)),
                    A36_090 = as.integer(round( A36_090*scale_factor)),
                    A36_100 = as.integer(round( A36_100*scale_factor)),
                    A36_110 = as.integer(round( A36_110*scale_factor)),
                    A36_120 = as.integer(round( A36_120*scale_factor)),
                    A36_130 = as.integer(round( A36_130*scale_factor)),
                    A36_140 = as.integer(round( A36_140*scale_factor)),
                    A36_150 = as.integer(round( A36_150*scale_factor)),
                    A36_160 = as.integer(round( A36_160*scale_factor)),
                    A36_170 = as.integer(round( A36_170*scale_factor)),
                    A36_180 = as.integer(round( A36_180*scale_factor)),
                    A36_190 = as.integer(round( A36_190*scale_factor)),
                    A36_200 = as.integer(round( A36_200*scale_factor)),
                    A36_210 = as.integer(round( A36_210*scale_factor)),
                    A36_220 = as.integer(round( A36_220*scale_factor)),
                    A36_230 = as.integer(round( A36_230*scale_factor)),
                    A36_240 = as.integer(round( A36_240*scale_factor)),
                    A36_250 = as.integer(round( A36_250*scale_factor)),
                    A36_260 = as.integer(round( A36_260*scale_factor)),
                    A36_270 = as.integer(round( A36_270*scale_factor)),
                    A36_280 = as.integer(round( A36_280*scale_factor)),
                    A36_290 = as.integer(round( A36_290*scale_factor)),
                    A36_300 = as.integer(round( A36_300*scale_factor)),
                    A36_310 = as.integer(round( A36_310*scale_factor)),
                    A36_320 = as.integer(round( A36_320*scale_factor)),
                    A36_330 = as.integer(round( A36_330*scale_factor)),
                    A36_340 = as.integer(round( A36_340*scale_factor)),
                    A36_350 = as.integer(round( A36_350*scale_factor)))]

    }


    # Name of the output file: tag tthe scaling factor ------------------------------------

    tag <- as.character(scale_factor)


    output_name <- gsub(pattern = ".rsf|.RSF",replacement = "",x = output_name) # remove the .rsf extension

    output_name <- paste("[scl_",tag,"]", output_name,sep = "")


    # Write the results into a new rsf file ------------------------------------------------

    Write_RSF(rsf, output_file_name = output_name)
}
