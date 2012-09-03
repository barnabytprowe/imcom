PRO example_write_config, USERXY=userxy

;
;    EXAMPLE_WRITE_CONFIG.PRO - writes a set of example config files for the
;                               prototype IMCOM package
;    Copyright (C) 2011 Barnaby Rowe
;
;    This program is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see <http://www.gnu.org/licenses/>.
;
;    DEPENDENCIES - IMCOM_WRITE_CONFIG.PRO
;                   IMCOM_WRITE_INCONFIG.PRO
;                   IMCOM_WRITE_OUTCONFIG.PRO
;
;    DESCRIPTION
;    Calls imcom_write_config.pro to generate example config files
;    ...Shows default use setting USERXY=true, all not-required keywords
;    are commented out
;

!PATH="../idl:"+!PATH

estring = "example"
if not keyword_set(userxy) then begin
  dithers = [[0.d0, 0.5d0, 0.d0, 0.5d0], [0.d0, 0.d0, 0.5d0, 0.5d0]] ; 2x2 dither pattern
  imcom_write_config, CONFIGFILE="config_"+estring,                                     $
                      NEXP=4,                                                           $
                      USERXY=userxy,                                                    $
                      GIMFILE=["test1.", "test2.", "test3.", "test4."]+estring+".fits", $
                      GIMXSCALE=[1.0d0, 1.0d0, 1.0d0, 1.0d0],                           $
                      GIMYSCALE=[1.0d0, 1.0d0, 1.0d0, 1.0d0],                           $
                      PSFFILE=replicate("psf."+estring+".fits", 4),                     $
                      PSFXSCALE=0.10d0,                                                 $
                      PSFYSCALE=0.10d0,                                                 $
                      ROTANGDEG=[0.d0, 0.d0, 0.d0, 0.d0],                               $
                      NOISE=[1.d0, 1.d0, 1.d0, 1.d0],                                   $
                      DITHERS=dithers,                                                  $
                      OUTFILE="H."+estring+".fits",                                     $
                      AFILE="A."+estring+".fits",                                       $
                      BFILE="B."+estring+".fits",                                       $
                      QFILE="Q."+estring+".fits",                                       $
                      LFILE="L."+estring+".fits",                                       $
                      PFILE="P."+estring+".fits",                                       $
                      KFILE="K."+estring+".fits",                                       $
                      TFILE="T."+estring+".fits",                                       $
                      SFILE="S."+estring+".fits",                                       $
                      UFILE="U."+estring+".fits",                                       $
                      GAMFILE="psf."+estring+".fits",                                   $
                      OUTXSCALE=0.5d0,                                                  $
                      OUTYSCALE=0.5d0,                                                  $
                      OUTPOS=[0.d0, 0.d0],                                              $
                      NOUT=[65L, 65L]
endif else begin
  imcom_write_config, CONFIGFILE="config_"+estring,                                     $
                      NEXP=4,                                                           $
                      USERXY=userxy,                                                    $
                      GIMFILE=["test1.", "test2.", "test3.", "test4."]+estring+".fits", $
                      GIMXFILE=["x1.", "x2.", "x3.", "x4."]+estring+".fits",            $
                      GIMYFILE=["y1.", "y2.", "y3.", "y4."]+estring+".fits",            $
                      PSFFILE=replicate("psf."+estring+".fits", 4),                     $
                      PSFXSCALE=0.10d0,                                                 $
                      PSFYSCALE=0.10d0,                                                 $
                      ROTANGDEG=[0.d0, 0.d0, 0.d0, 0.d0],                               $
                      NOISE=[1.d0, 1.d0, 1.d0, 1.d0],                                   $
                      OUTFILE="H."+estring+".fits",                                     $
                      AFILE="A."+estring+".fits",                                       $
                      BFILE="B."+estring+".fits",                                       $
                      QFILE="Q."+estring+".fits",                                       $
                      LFILE="L."+estring+".fits",                                       $
                      PFILE="P."+estring+".fits",                                       $
                      KFILE="K."+estring+".fits",                                       $
                      TFILE="T."+estring+".fits",                                       $
                      SFILE="S."+estring+".fits",                                       $
                      UFILE="U."+estring+".fits",                                       $
                      GAMFILE="psf."+estring+".fits",                                   $
                      OUTXFILE="xout."+estring+".fits",                                 $
                      OUTYFILE="yout."+estring+".fits"
endelse

END
