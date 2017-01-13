
## Function list
Function are split into two groups:

+ **Master functions**: Work directly with files (as inputs and outputs). They are the "common"" functions directly call by the user.
+ **Root functions**: Basic functions (work with object already loaded into memory). They are called (indirectly) via the master functions.



### Master functions:
  
  * **RSF_convert <- function(rsf_file1,
                         output_name = rsf_file1,
                         crop = FALSE,
                         layout_csv = NULL,
                         crop_buffer = 700,
                         twelve_S = FALSE,
                         rsf_file2 = NULL,
                         shear_out = FALSE,
                         three_D = FALSE,
                         layer_H = NULL)**    
Depending on the provided inputs, do basic transformation on *.rsf file (rsf = wind ReSource File): Crop, convert in 12 sector, create new height levels, compute and export shear table.
    + rsf_file1:	Input *.rsf file name (name should be surounded by "").
    + output_name:	Output *.rsf file name " (name should be surounded by "").
    + crop:	 Should the rsf be croped ? TRUE/FALSE. Default value = FALSE.
    + layout_csv:	 Input *.csv file name (name should be surounded by "") containing the WTG positions (X,Y).Default value = NULL. This input is used to crop and benchmark the production (12 Sectors versus 36 sectors).
    + crop_buffer:	 Integer value: Distance in meters around the WTG postions that should be crop (default value = 700m)
    + twelve_S:	 Should the rsf be converted into 12 sectors ? TRUE/FALSE. Default value = FALSE.
    + rsf_file2:	Input *.rsf file name (name should be surounded by ""). Second rsf file used for shear and 3 dimensionals rsf creation. Default value = NULL.
    + shear_out:	Should the shear should be extracted and exported (csv file) ? TRUE/FALSE. Default value = FALSE.
    + three_D:	Should a three dimensional rsf created ? TRUE/FALSE. Default value = FALSE.
    + layer_H:	vector containing the level of the 3 dimensional rsf (used if three_D = TRUE). if not provided and three_D = TRUE, the levels multiple of 10 would be created. Default value = NULL.
```r
            #Example 1: 
            # Transform one rsf file (crop and reduce the number of sectors)
            RSF_convert(rsf_file1 = "RSF-windresource-CFD_Hennoid [80m].rsf",
            crop = TRUE,
            layout_csv = "layout V2 Hennoid.csv",
            twelve_S = TRUE,
            output_name = "Hennoid")
            
            #Example 2: 
            # Create a 3 dimensional rsf based on two rsf files.
            RSF_convert(rsf_file1 = "RSF-windresource-CFD_Hennoid [80m].rsf",
            rsf_file2 = "RSF-windresource-CFD_Hennoid_[117].rsf",
            crop = TRUE,
            layout_csv = "layout V2 Hennoid.csv",
            shear_out = TRUE,
            three_D = TRUE,
            layer_H = c(94,105,112),
            output_name = "Hennoid")
```
  
  * **RSF_plot <- function(rsf_file,layout_csv = NULL)**  
Raster plot of A, k and wind speed (Average or Sectorwise). If provided it will also plot the layout.
    + rsf_file:	Input *.rsf file name (name should be surounded by "").
    + layout_csv:	Input *.csv file name (name should be surounded by "") containing the WTG positions (X,Y).Default value = NULL. if provided, turbine postions would be plotted on the map.
```r            
            Example: 
            RSF_plot(rsf_file1 = "[3D][crop_700]Hennoid.rsf") 
```            
            
   * **Shear_plot <- function(shear_file,layout_csv = NULL)**  
    Raster plot of wind shear (Average or Sectorwise). If provided it will also plot the layout.
    + rsf_file:	Input *.rsf file name (name should be surounded by "").
    + layout_csv:	Input *.csv file name (name should be surounded by "") containing the WTG positions (X,Y).Default value = NULL. if provided, turbine postions would be plotted on the map.
```r
            Example: 
            Shear_plot(rsf_file1 = "[shear_table][crop_700]Hennoid.csv", layout_csv = "layout V2 Hennoid.csv") 
```      
      
  * **ShearExtract <- function(shear_csv,layout_csv, output_name = "shear_VSC_input.shr")**  
      Extract shear at turbine locations, and create a *.shr file for VSC (use as CFD input)  
    + shear_csv: shear file in CSV format (output of "ShearTable" function)  
    + layout_csv: layout file in csv format (X,Y)  
    + output_name: name of the output file (*.shr)  
```r        
            Example: 
            ShearExtract("sheartable.csv","layout_Saint_Martin.csv")
```          


### Root functions:
"Root" functions work with object already loaded on the memory. They are core modules of this tool.


 - **Read_RSF <- function("Input_file")**  
     Read an *.rsf ("Input_file"). it recognizes if it is 36 or 12 sectors, load it accordingly into memory, with explicit columns names
```r     
            #Example: 
            rsf87 <- Read_RSF("RSF-windresource-CFD_Tonstaad [87m].rsf")
```
 - **read.csv**  
      Standard R function to read and load into memory a csv file. use here to load the layout
```r      
            #Example: 
            layout_T <- read.csv("[V10A] Tonstaad layout.csv")
```
 - **Write_RSF <- function(rsf, "output_file_name")**  
     Write an rsf table (rsf) in memory to a *.rsf file ("output_file_name")
```r     
            #Example: 
            Write_RSF(rsf87_12S,"[12S][87]Tonstaad.rsf")
```
- **Crop_RSF <- function(rsf, Xmax, Xmin, Ymax, Ymin)**  
    Crop rsf with new boundaries.
```r    
            #Example: 
            rsf87C2 <- Crop_RSF(rsf87,645000,635000,598800,597800)
```
- **S36_to_S12_RSF <- function(rsf36)**  
     Take as argument a 36 sectors rsf, and transform it into 12 sectors rsf.
     it combines the neighboring sectors 3 by 3 by doing a ponderate average.
```r     
            #Example: 
            rsf87_12S <- S36_to_S12_RSF(rsf87)
```
- **Gen_power_curve <- function(turbine = "V117-3.45")**  
      Generate a power curve data frame, that can be used in other function. to match the dim of other function, the power curve need to go from 0ms to 40ms by 0.5ms. (to dev: add other turbine types). For the moment only the V117-3.45 is available (or benchmark purpose it should be enough)
```r
            #Example: 
            V117_PC <- Gen_power_curve(V117)
```
- **Prod_RSF <- function(rsf, Xp, Yp,power_curve)**  
     Compute the production [MWh] for one turbine type (power_curve) on one position (Xp, Yp) based on the rsf (rsf).
```r     
            #Example: 
            rod_RSF(rsf87, 644354,544568, V117_PC)
```
- **Bensh_RSF <- function(rsf36, rsf12, layout, power_curve)**  
    Compute the ratio of AEP from 36S (rsf36) and 12S (rsf12) rsf. it returns a table of production for each turbine positions (layout,power_curve).
    To do so it call the Prod_RSF.
```r    
            #Example: 
            Bensh_table <- Bensh_RSF(rsf87, rsf87_12S, layout_T, V117_PC)
```  
- **Interpol_RSF <- function(rsf_H1,rsf_H2,layer_H)**  
    create a combined rsf file for the height(s) provided in the vector "layer_H", based on the interpolation of rsf at H1 and H2. it assumes a shear profile of the wind speed.  The function can actually extrapolate if needed. It is recommended that the difference H1-H2 is at least 20m.
```r    
            #Example: 
            rsf3D <- Interpol_RSF(rsf87_12S, rsf137_12S, c(91.5,94,117,125,147))
```  
- **Shear_RSF <- function(rsf_H1,rsf_H2)**  
    Create a Shear map (in rsf format) based on the two input rsf files (rsf_H1,rsf_H2). It is recommended to have abs(H1-H2)>20m. The output format is also in a rsf format.
```r    
            #Example:
            shear_ <- Shear_RSF(rsf87_12S,rsf137_12S)
```    
- **Shear2png <- function(RSF_shear)**
    Save the "RSF_shear" rsf to a serie (average and for all sectors) of georefenced images (*.png and *.pgw). 
    for more information about the *.pgw (word file): https://en.wikipedia.org/wiki/World_file). those images can be      directly used in windpro as background maps.
```r    
            #Example:
            Shear2png(shear_) 
```    
- **ShearTable <- function(RSF_shear, output_name = "sheartable.csv")**  
    Re-format the shear table (result from Shear_RSF) from (x,y,sector) table (3 x n) --> ((x,y) ~ sector) table (2+n_sector x n).it then writes the result in a csv file (can be used afterward to extrapolate to new HH or extract the shear value on some positions)
```r
            #Example:
            ShearTable(shear_) 
```            
Note: the last tree functions on the shear (Shear_RSF, Shear2png, ShearTable) would be merge in the future.
    
