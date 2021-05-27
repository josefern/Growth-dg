                          
                                   U S E R  M A N U A L 

                                    "G R O W T H - dg"                                                           
          
            A 3-D inversion tool for modelling gravity changes on a small network.   
     
    Reference:  

         Antonio G. Camacho, Peter Vajda, Craig A. Miller, José Fernández 
         A free-geometry geodynamic modelling of surface gravity changes using Growth-dg software
         S.D.  2021.           

         Corresponding author: José Fernández (jft@mat.ucm.es)




1. INTRODUCTION


   This computational tool, Growth-dg, consist on a Fortran 77 code.
   It calculates in an almost automatic way the possible sources of density change responsible for the observed gravity variations. 
   The process is based on a non-linear adjustment of anomalous structures defined by aggregation of small filled cells. 
   It works on a step by step growth process that fills prismatic cells of a previous 3D grid, and incorporates them to the growing anomalous structure.

   Some characteristics of the approach are the following: 

    (1) solution for 3-D model space, 
    (2) acceptance of non-gridded non-planar low-precision imprecise data, 
    (3) optional determination of an offset gravity component in the data, 
    (4) managing of outlier values by means of weigthing accross the iterative process,
    (5) simultaneous inversion for both positive and negative density contrasts, 
    (6) inversion for irregularly-shaped structures composed of several individual bodies, 
    (7) non homogeneous density distribution, highligting the core anomalous structures
    (8) semi-automated data inversion routine, 
    (9) several graphical outputs (process and results) helping to take modelling decisions,
   (10) only two main modelling values: model size and smoothing


   The main input file, "grad.txt", contains the coordinates (m) of the gravity stations and the observed gravity changes (uGal).


        
           
2. INPUTS 
                                                            
2.1. Input files:

  -> file "grad.txt": The basic input ASCII file, named "grad.txt" (as default), contains a line corresponding to each gravity benchmark with the following values: 
     x-coordinate (m), y-coordinate (m), altitude (m) (above sea level), gravity change (uGal) (1 uGal=10-8 m/s2)

      example:

    212876.1 3079402.0  560.30  -4 
    212903.9 3079303.6  571.02   1 
    213149.2 3079854.5  533.59 -36 
    213128.0 3080707.0  524.20 -59 
    212667.0 3081849.6  403.01 -18 
    ...............................
    201634.6 3073315.1  283.73  -9
     0 0 0 0 

    The end can be denoted by a line of zeros.


  -> file "map.bln" (default name): ASCII file containing a base map (coast, roads, etc.) to be included on the graphic presentation. 
   The format of this file is:
 
       j1
     xx(1),yy(1)
     xx(2),yy(2)
     ..........
     xx(j1),yy(j1)
        j2
     xx(1),yy(1) 
     xx(2),yy(2)
     ...........
     xx(j2),yy(j2) 
       ...

   where j1,j2,... are the number of points for each line and xx,yy are plane coordinates for the sucessive points of the lines.
       
      example

       10 
    491925 4280999
    491984 4280863
    .............
    493287 4279100 




                 
2.2. Values for several parameters in a dialog form through the Graphical User Interface (GUI): 


  First step (top left hand pane): Input and output files

  - Names of the input and output files (default values are provided). 
    Paths can be used. It is necessary to press “Run/OK” at the end of each step. 
    The program’s response to successful data loading is displaying the number count of benchmarks.



  Second step (bottom left hand pane): 3D grid of cells

 - Mean length (m) for the side of the parallelepiped cells (default value is provided). 
   Shallow cells will have smaller length and deep cells will have larger length. 
   The program’s response to successful cell generation is displaying the number count of resulting cells.



  Third step (right hand pane): values for the inversion parameters.
    (default value are provided)

 a) Gravity offset and regional trend.
     - By pressing the corresponding “g0 Adj” button (default option), 
       an additional offset value g0 is adjusted as part of the inversion. 
       The offset parameter has a particular bearing on the deeper parts of the model domain. 
       We suggest trying with and without offset option.
     - By pressing the corresponding “gx,gy Adj” buttont, specially for very large data set, 
       an additional linear trend values gx, gy (µGal/km) are adjusted in the inversion approach. 
       This option has been included to allow the use of this software for CBA anomaly data.

 b) Coefficient r for random search. This is an integer value (>=1) that allows for an exploratory process with a more or less random character. 
    For a value of 1, the program carries out a systematic exploration. 
    For a value r>1, the program tests 100/r % cells, randomly selected, thus allowing for a faster (but less precise) solution.

 c) Balance factor. This positive value corresponds to the parameter lambda for balance between “fitness” and “smoothness” of the model. 
    The value is dependent on the quality of the data and on the complexity of the subsurface structure. 
    For low values, the inversion results in the generation of an excessively complex model of the subsurface structure with very good fit to the data. 
    For high values, the inversion produces a simple model but with poor fit to the data. 
    The operator may try different values of it close to the default value and draw resulting models and residuals.

 d) Ending values.
   There are two alternative parameters to define the normal end of the inversion process.
    - A minimum avarage value for density increase (second square). 
      It is active with a non-null value.
    - (Default) A value for relative size (%) of the resulting anomalous model volume with respect 
      to the available subsoil volume in the 3D grid partition.

 e) Number of possible density levels for any cells. It is the number of possible revisiting on the cells. 
   For value 1 the model becomes homogeneous.

 f) Blunder value for detecting outlier values and for further weighting in the growth process.

 g) Upward weighting factor. It is an optional parameter to get shallower models. 
   It could be useful for models characterized by a shallow structures (sedimentary structures). 
   Value cero for u corresponds to neutral weighting. Higher value as 1,2, ... produce shallower models. 




3. RUNNING

   Once all the required parameters are introduced, the program initiates the inversion. 
   The evolution of the numbers of filled cells, density contrast values, standard deviation, 
   offset value, etc. are shown on-screen. 
   Several graphical presentations (plan and elevation views) shows the model growth and the evolution of modelled and residual values.

   The process ends when the ending value (mode size or density contrast) is reached 
   or when it is not possible to aggregate a new cell within the model constrains. 




4. OUTPUTS

   The program shows several screens with numerical and graphical information about the inversion process and about the results. 
   The final model is shwon by means of plan and elevation views and by means of several horizontal and vertical sections with colour code.


   In addition, two output files with results are offered: "mod.txt" and "fil.txt"                                                            
                                                               
 -> file "mod.txt": inversion model as 3-D partition of the subsurface volume into filled cells of anomalous density. 
    The content and general format of this ASCII file is: 
 
    nc, side, lambda, model size , dens 
    x(1),y(1),z(1),sx(1),sy(1),sz(1),den(1),s(1),f(1)   
    x(2),y(2),z(2),sx(2),sy(2),sz(2),den(2),s(2),f(2)          
    ..........................................
    x(nc),y(nc),z(nc),sx(nc),sy(nc),sz(nc),den(nc),s(nc),f(nc)

    where:                                                         

    nc:  total number of filled cells (parallelepiped) of the model,          
    side: adopted mean side (in metres) for the cells in the 3D partition,    
    lambda: adopted lambda balance value, 
    model size: % of the available space volume ocupied by model cells,
    dens: assumed a priori density contrast,    
    x(j),y(j),z(j): coordinates (in metres) UTM and depth for the center of the j-th cell of the model,
    sx(j),sy(j),sz(j): horizontal sides along x,y,and z of the j-th cell of the model,
    den(j): anomalous density contrast (kg/m3) of the j-th cell,
    s(j): relative sensitivity for each cell from the survey points,
    f(j): density value level at the moment of filling that cell

    j=1,...,nc


      example:
 
    1775     cell side=  178.    lamda=  30.    %-model= 3.0    dens=  0.0
     X(m)   Y(m)   Z(m asl)   sx(m)  sy(m)   sz(m)  den(kg/m3)   sen    ord
   368625 6009136    2229      213     96    213      3.6        9.7    3.8
   368838 6009125    2229      213     98    213      3.6        9.4    3.7
   369050 6009183    2229      213     97    213      3.6        9.5    4.0
   369263 6010017    2229      213     84    213      3.6        8.3    3.6
   360039 6011656    2069      108    138    108     -7.1        3.5    4.0
    . . . . . .. . . . . .


 -> file "fil.txt": ASCII file containing, for each data point (in planar coordinates), the following values (uGal) resulting from the inversion process:
              offset (regional) gravity anomaly
              observed residual gravity (= observed gravity change - offset value)
              modelled gravity change
              final residual values
              standardized residuals and test for outlier

     example:

       X       Y         Off      dG    Mod-dg    Res-dg    Test
    11000.  362000.      0.0    193.0    136.8     56.2     1.03        
    13023.  359471.      0.0     30.3     28.7      1.6     0.03        
     8586.  362467.      0.0    453.1    366.1     87.0     1.59        
     7741.  363047.      0.0    334.1    259.6     74.5     1.36        
   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .




5. DIMENSIONS                                                         

    The memory requirements of the program execution depend on the following parameters:
                                                         
       ms:  max. number of gravity stations. 
       mc:  max. number of whole prismatic cells for the 3-D grid.   

   The values of these size parameters can be changed in the code according to the applications and the memory abailability.





6. COMPILING SPECIFICATIONS. 

   Compiled with Intel Visual Fortran, by using a Quickwin Application Project. 
   The complete compilation requires reveral additional files (that are available in the zip file of codes and):
     
   Resource files:  Growth-dg.rc    (resource script)
                    Growth-dg.ico   (icono)
   
   Fortran source files:  Growth-dg.fd
                          Growth-dg.for 

   C/C++Header:  Growth-dg.h   (dialog window)    






     