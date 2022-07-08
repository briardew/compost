#!/bin/csh

# COMPOST  COMParison of Observations to Simulated Tracers
#
# Author(s):	Brad Weir (brad.weir@nasa.gov)
#
# Changelog:
# 2019/03/27	New version, based on compare_obspack_loop.csh
# 2021/12/17	Major updates
#
# Todo:
# * Run comparison jobs in parallel? (needs output change)
# * Utilities to write output in read format
# * Rewrite C shell parts in bash or python
# * Seek forgiveness for sins: using C shell, etc.
#===============================================================================

#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --job-name=compost
#SBATCH --account=s1460
#SBATCH --qos=long

##SBATCH --ntasks=1
##SBATCH --time=12:00:00
##SBATCH --job-name=compost
##SBATCH --account=s1460
##SBATCH --qos=covid19_s1460
##SBATCH --partition=compute

#SBATCH --output=compost.log.o%j
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=briardew@gmail.com

# 1. SETUP ENVIRONMENT
#===============================================================================
alias mcons /discover/vis/mathworks/matlab_r2020a/bin/matlab -nosplash -nodesktop

# Grab a license
set nolic = 1
while ("$nolic" == "1")
  mcons -r "exit"
  set nolic = $?
  sleep 10
end
## Don't think this does what I want it to (check)
#mcons -r "while 1, pause(60); end; exit;" &

# Load run settings
# -----------------
if (${#argv} == 1 && -f ${argv[1]}) then
  set frun = ${argv[1]}

  set drun = `dirname  ${argv[1]} .rc`
  set brun = `basename ${argv[1]} .rc`
  set orig = `pwd`

# Default values, can override in input settings
  set DIROUT = "${drun}/${brun}"

  source ${frun}
else
  echo "You must provide a single rc file as a command line argument ..."
  echo "Try looking in the etc directory"
  exit 1
endif

# Create work directory
# ---------------------
set dbit = `date +"${DIROUT}/%Y%m%d"`

foreach tag ( "a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m" \
              "n" "o" "p" "q" "r" "s" "t" "u" "v" "w" "x" "y" "z" )
  set dout = "${dbit}${tag}"
  if (! (-d ${dout})) break
end

mkdir -p ${dout}
mkdir    ${dout}/molefracs
cp ${frun} ${dout}/compost.rc
cp -r +compost/    ${dout}
cp -r +fit/        ${dout}
cp -r +atmosutils/ ${dout}

echo "#==============================================================================="
echo "#"
echo "#   COMParison of Observations to Simulated Tracers (COMPOST)"
echo "#"  
echo "#==============================================================================="
echo ""

echo "Running in directory ${dout} ..."
echo ""

cd ${dout}


# 2. GENERATE COMPOST SCRIPTS
#===============================================================================
@ NFIELDS = 9
@ NCOMPS  = ${#OBSINFO} / ${NFIELDS}
if (${#OBSINFO} != ${NFIELDS} * ${NCOMPS}) then
   echo ""
   echo "Incorrect specification of OBSINFO variable, exiting ..."
   exit 1
endif

# Generate compare scripts
# ------------------------
set join_typeid = ()
set join_dataid = ()
foreach nn (`seq 1 $NCOMPS`)
#  Read next line of OBSINFO
   @ mm = $NFIELDS * ($nn - 1)

   @ mm++; set typeid  = ${OBSINFO[$mm]}
   @ mm++; set obstype = ${OBSINFO[$mm]}
   @ mm++; set dataid  = ${OBSINFO[$mm]}
   @ mm++; set dirobs  = ${OBSINFO[$mm]}
   @ mm++; set varobs  = ${OBSINFO[$mm]}
   @ mm++; set sclobs  = ${OBSINFO[$mm]}
   @ mm++; set varmod  = ${OBSINFO[$mm]}
   @ mm++; set sclmod  = ${OBSINFO[$mm]}
   @ mm++; set isdry   = ${OBSINFO[$mm]}

#  Hack for portability
   set compfn = compare
   if ("$obstype" != "-") then
     set compfn = ${obstype}_compare
   endif

#  Fill template with user settings
   cat +compost/${typeid}_${compfn}.tmpl.m		| \
       sed -e "s?>>>ISDRY<<<?${isdry}?"			| \
       sed -e "s?>>>SCLOBS<<<?${sclobs}?"		| \
       sed -e "s?>>>SCLMOD<<<?${sclmod}?"		| \
       sed -e "s?>>>VAROBS<<<?'${varobs}'?"		| \
       sed -e "s?>>>VARMOD<<<?'${varmod}'?"		| \
       sed -e "s?>>>VARPHIS<<<?'${VARPHIS}'?"		| \
       sed -e "s?>>>VARPS<<<?'${VARPS}'?"		| \
       sed -e "s?>>>VARZL<<<?'${VARZL}'?"		| \
       sed -e "s?>>>VARQW<<<?'${VARQW}'?"		| \
       sed -e "s?>>>TTMET<<<?'${TTMET}'?"		| \
       sed -e "s?>>>TTGAS<<<?'${TTGAS}'?"		| \
       sed -e "s?>>>DIROBS<<<?'${dirobs}'?"		| \
       sed -e "s?>>>DIRMOD<<<?'molefracs/'?"		| \
       sed -e "s?>>>HDMET<<<?'${EXPID}.${COLL}.'?"	| \
       sed -e "s?>>>HDGAS<<<?'${EXPID}.${COLL}.'?"	> \
       +compost/${typeid}_${dataid}_${compfn}.m

#  Compile list of unique comparison ids
   if (" $join_dataid " !~ *" $dataid "*) then
      set join_typeid = ($join_typeid $typeid)
      set join_dataid = ($join_dataid $dataid)
   endif
end
set NJOIN = ${#join_dataid}

# Generate join scripts
# ---------------------
foreach kk (`seq 1 $NJOIN`)
   set typeid = ${join_typeid[$kk]}
   set dataid = ${join_dataid[$kk]}
   set headid = ${EXPID}__${typeid}_${dataid}

   cat +compost/${typeid}_join.tmpl.m			| \
       sed -e "s?>>>HEADID<<<?${headid}?"		| \
       sed -e "s?>>>YEAR0<<<?${YEAR0}?"			| \
       sed -e "s?>>>YEARF<<<?${YEARF}?"			> \
       +compost/${typeid}_${dataid}_join.m
end


foreach nyear (`seq ${YEAR0} ${YEARF}`)
# 3. COLLECT MODEL DATA
#===============================================================================
   echo "#"  
   echo "#        Collecting model files for ${nyear}"
   echo "#==============================================================================="
   echo ""

   cd molefracs

   @ nprev = $nyear - 1
   @ nnext = $nyear + 1

   if ($UNTAR == 0) then
#     Find the files and create symbolic links
      find ${DIRAR} -name "${EXPID}.${COLL}.${nprev}1231_*z.nc4" \
           -exec ln -s {} . \;
      find ${DIRAR} -name "${EXPID}.${COLL}.${nyear}????_*z.nc4" \
           -exec ln -s {} . \;
      find ${DIRAR} -name "${EXPID}.${COLL}.${nnext}0101_*z.nc4" \
           -exec ln -s {} . \;
   else
#     Archive directories (must conform with run settings)
      set DPRV = ${DIRAR}/${COLL}/Y${nprev}
      set DNOW = ${DIRAR}/${COLL}/Y${nyear}
      set DNXT = ${DIRAR}/${COLL}/Y${nnext}

#     Untar archived data
      tar xf ${DPRV}/${EXPID}.${COLL}.daily.${nprev}12.nc4.tar \
          --wildcards "${EXPID}.${COLL}.${nprev}1231_*z.nc4" &
      set tarids = $!
      foreach nmon (`seq -w 01 12`)
         tar xf ${DNOW}/${EXPID}.${COLL}.daily.${nyear}${nmon}.nc4.tar &
         set tarids = ( $tarids $! )
      end
      tar xf ${DNXT}/${EXPID}.${COLL}.daily.${nnext}01.nc4.tar \
          --wildcards "${EXPID}.${COLL}.${nnext}0101_*z.nc4" &
      set tarids = ( $tarids $! )

#     Some logic to kill tar jobs in case of interrupts
      onintr stoptar
      while (`ps -p "$tarids" | wc -l` > 1)
         sleep 10
      end
      set tarids = ()
   endif

   cd ..


# 4. RUN COMPARISONS IN MATLAB
#===============================================================================
   echo ""
   echo "#"  

   set t0 = `date +"%s"`
   foreach nn (`seq 1 $NCOMPS`)
      @ mm = $NFIELDS * ($nn - 1) + 1; set typeid  = ${OBSINFO[$mm]}
      @ mm = $NFIELDS * ($nn - 1) + 2; set obstype = ${OBSINFO[$mm]}
      @ mm = $NFIELDS * ($nn - 1) + 3; set dataid  = ${OBSINFO[$mm]}

      set compid = ${typeid}_${dataid}
      if ("$obstype" != "-") then
        set compid = ${compid}_${obstype}
      endif

      set frun = ${EXPID}__${compid}__${nyear}.m		# temporary
      set fout = ${EXPID}__${compid}__${nyear}.mat
      set flog = ${EXPID}__${compid}__${nyear}.log

      echo "#        Running ${compid}_compare for ${nyear}: ${flog}"

#     Extra work to get the logging, etc. right
      echo "compost.${compid}_compare; save('${fout}'); exit" > ${frun}
      mcons < ${frun} >& ${flog} &
      rm ${frun}
   end

   wait
   set tF = `date +"%s"`
   @ dmin = ($tF - $t0) / 60

   echo ""
   echo "#        ${nyear} completed (${dmin} minutes elapsed)"
   echo "#==============================================================================="

#  Clean up
   rm -f molefracs/${EXPID}.${COLL}.*.nc4
end
rmdir molefracs


# 5. JOIN YEARLY COMPARISONS
#===============================================================================
foreach kk (`seq 1 $NJOIN`)
   set typeid = ${join_typeid[$kk]}
   set dataid = ${join_dataid[$kk]}
   set joinid = ${typeid}_${dataid}

   echo ""
   echo "#"  
   echo "#        Joining ${joinid} comparison"
   echo "#==============================================================================="
   echo ""

   set fout = ${EXPID}__${joinid}.mat

   mcons -r "compost.${joinid}_join; save('${fout}'); exit"

   echo ""
   echo "#        Created comparison file ${fout}"
   echo "#==============================================================================="
   echo ""
end

cd ${orig}

echo ""
echo "#"  
echo "#        THE END: Compost output in ${dout}"
echo "#==============================================================================="
echo ""

exit 0

stoptar:
   echo "Killing untarring processes and exiting ..."
   kill $tarids
   exit 1
