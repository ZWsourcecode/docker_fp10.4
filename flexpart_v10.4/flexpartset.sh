#!/bin/bash 
### create pathnames, RELEASES, COMMAND file for flexpart
### script should be in parent folder of options
### run : bash flexpartset.sh 
### zhendong.wu@nateko.lu.se

HOME_PATH=$(pwd)
OPTION_PATH=$HOME_PATH/options

SPINUPTIME=10 # in days

echo ------ make pathnames, RELEASES, COMMAND file ------
echo

# ======================================
# arguments of flexpartset.sh
# ======================================
input=$1 # meteo input path
output=$2 # flexpart output path
id=$3 # species id, NOTE one id for now 
start=$4 # start date, e.g. YYYYMMDD
end=$5 # end date, e.g. YYYYMMDD
step=$6 # time step in hour
lon=$7 # longitude of release box -180 < LON1 <180
lat=$8 # latitude of release box, -90 < LAT1 < 90
z=$9 # height of release
particles=${10} # Total number of particles to be released

# ======================================
# backup orignal files
# ======================================
cp -n ${OPTION_PATH}/RELEASES ${OPTION_PATH}/RELEASES_backup
cp -n ${OPTION_PATH}/COMMAND ${OPTION_PATH}/COMMAND_backup

# ======================================
# make pathnames file
# ======================================
cat <<EOM >$HOME_PATH/pathnames
./options/
${output}
${input}
${input}AVAILABLE
EOM
echo pathnames is done, under $HOME_PATH/
echo

# ======================================
# make RELEASES file
# ======================================
cat <<EOM >$OPTION_PATH/RELEASES
***************************************************************************************************************
*                                                                                                             *
*                                                                                                             *
*                                                                                                             *
*   Input file for the Lagrangian particle dispersion model FLEXPART                                          *
*                        Please select your options                                                           *
*                                                                                                             *
*                                                                                                             *
*                                                                                                             *
***************************************************************************************************************
&RELEASES_CTRL
 NSPEC      =           1, ! Total number of species
 SPECNUM_REL=          ${id}, ! Species numbers in directory SPECIES
 /
EOM

start_release=$(date -u -d "$start -$step hours")
end_release=$(date -u -d "$end -$step hours")

# i=-$step
# while [ "$thedate" != "$end_release" ]
# do
#  	thedate=$(date -u -d "$start +$i hours")
# # 	echo $thedate
# 	IDATE1=$(date -u -d "$start +$i hours" +%Y%m%d)
# 	ITIME1=$(date -u -d "$start +$i hours" +%H%M%S)
# # 	echo $IDATE1
# # 	echo $ITIME1
# 	
# 	i=$(( i + $step ))
# 	
# 	IDATE2=$(date -u -d "$start +$i hours" +%Y%m%d)
# 	ITIME2=$(date -u -d "$start +$i hours" +%H%M%S)
# # 	echo $IDATE2
# # 	echo $ITIME2

i=0
while [ "$thedate" != "$start_release" ]
do
 	thedate=$(date -u -d "$end -$i hours")
# 	echo $thedate
	IDATE2=$(date -u -d "$end -$i hours" +%Y%m%d)
	ITIME2=$(date -u -d "$end -$i hours" +%H%M%S)
	
# 	echo $IDATE1
# 	echo $ITIME1
	
	i=$(( i - $step ))
	
	IDATE1=$(date -u -d "$end -$i hours" +%Y%m%d)
	ITIME1=$(date -u -d "$end -$i hours" +%H%M%S)
	
# 	echo $IDATE2
# 	echo $ITIME2
	
cat <<EOM >>$OPTION_PATH/RELEASES
&RELEASE                   ! For each release 
IDATE1  =       $IDATE1, ! Release start date, YYYYMMDD: YYYY=year, MM=month, DD=day
ITIME1  =         $ITIME1, ! Release start time in UTC HHMISS: HH hours, MI=minutes, SS=seconds
IDATE2  =       $IDATE2, ! Release end date, same as IDATE1
ITIME2  =         $ITIME2, ! Release end time, same as ITIME1
LON1    =          $lon, ! Left longitude of release box -180 < LON1 <180
LON2    =          $lon, ! Right longitude of release box, same as LON1
LAT1    =         $lat, ! Lower latitude of release box, -90 < LAT1 < 90
LAT2    =         $lat, ! Upper latitude of release box same format as LAT1 
Z1      =         $z, ! Lower height of release box meters/hPa above reference level
Z2      =         $z, ! Upper height of release box meters/hPa above reference level
ZKIND   =              1, ! Reference level 1=above ground, 2=above sea level, 3 for pressure in hPa
MASS    =       1.0000E0, ! Total mass emitted, only relevant for fwd simulations
PARTS   =          $particles, ! Total number of particles to be released
COMMENT =    "RELEASE 1 $IDATE1:$ITIME1", ! Comment, written in the outputfile
/
EOM
	
done
echo RELEASES file is done, under $OPTION_PATH/. i.e. $(( $i / $step + 2 )) releases from $start to $end. 
echo



# ======================================
# make COMMAND file
# ======================================
IBDATE=$(date -u -d "$start -$SPINUPTIME days" +%Y%m%d)
IBTIME=$(date -u -d "$start -$SPINUPTIME days" +%H%M%S)
IEDATE=$(date -u -d "$end" +%Y%m%d)
IETIME=$(date -u -d "$end" +%H%M%S)
OUTSTEP=$(( $step * 60 * 60 )) # in seconds

cat <<EOM >$OPTION_PATH/COMMAND
***************************************************************************************************************
*                                                                                                             *
*      Input file for the Lagrangian particle dispersion model FLEXPART                                       *
*                           Please select your options                                                        *
*                                                                                                             *
***************************************************************************************************************
&COMMAND
 LDIRECT=               -1, ! Simulation direction in time   ; 1 (forward) or -1 (backward)
 IBDATE=         $IBDATE, ! Start date of the simulation   ; YYYYMMDD: YYYY=year, MM=month, DD=day  
 IBTIME=           $IBTIME, ! Start time of the simulation   ; HHMISS: HH=hours, MI=min, SS=sec; UTC
 IEDATE=         $IEDATE, ! End date of the simulation     ; same format as IBDATE 
 IETIME=           $IETIME, ! End  time of the simulation    ; same format as IBTIME
 LOUTSTEP=           $OUTSTEP, ! Interval of model output; average concentrations calculated every LOUTSTEP (s)  
 LOUTAVER=           $OUTSTEP, ! Interval of output averaging (s)
 LOUTSAMPLE=          900, ! Interval of output sampling  (s), higher stat. accuracy with shorter intervals
 ITSPLIT=        99999999, ! Interval of particle splitting (s) 
 LSYNCTIME=           900, ! All processes are synchronized to this time interval (s)
 CTL=          -5.0000000, ! CTL>1, ABL time step = (Lagrangian timescale (TL))/CTL, uses LSYNCTIME if CTL<0
 IFINE=                 4, ! Reduction for time step in vertical transport, used only if CTL>1 
 IOUT=                  13, ! Output type: [1]mass 2]pptv 3]1&2 4]plume 5]1&4, +8 for NetCDF output     
 IPOUT=                 1, ! Particle position output: 0]no 1]every output 2]only at end 3]time averaged
 LSUBGRID=              1, ! Increase of ABL heights due to sub-grid scale orographic variations;[0]off 1]on 
 LCONVECTION=           1, ! Switch for convection parameterization;0]off [1]on
 LAGESPECTRA=           1, ! Switch for calculation of age spectra (needs AGECLASSES);[0]off 1]on
 IPIN=                  0, ! Warm start from particle dump (needs previous partposit_end file); [0]no 1]yes  
 IOUTPUTFOREACHRELEASE= 1, ! Separate output fields for each location in the RELEASE file; [0]no 1]yes 
 IFLUX=                 0, ! Output of mass fluxes through output grid box boundaries
 MDOMAINFILL=           0, ! Switch for domain-filling, if limited-area particles generated at boundary
 IND_SOURCE=            1, ! Unit to be used at the source   ;  [1]mass 2]mass mixing ratio 
 IND_RECEPTOR=          1, ! Unit to be used at the receptor; [1]mass 2]mass mixing ratio 3]wet depo. 4]dry depo.
 MQUASILAG=             0, ! Quasi-Lagrangian mode to track individual numbered particles 
 NESTED_OUTPUT=         0, ! Output also for a nested domain 
 LINIT_COND=            0, ! Output sensitivity to initial conditions (bkw mode only) [0]off 1]conc 2]mmr 
 SURF_ONLY=             0, ! Output only for the lowest model layer, used w/ LINIT_COND=1 or 2
 CBLFLAG=               0, ! Skewed, not Gaussian turbulence in the convective ABL, need large CTL and IFINE
 OHFIELDS_PATH= "../../flexin/", ! Default path for OH file
 /
EOM

echo COMMAND file is done, under $OPTION_PATH/.  
echo
