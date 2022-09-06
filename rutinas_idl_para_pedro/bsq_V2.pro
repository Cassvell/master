;
;Name:
;	bsq_V2.pro
;purpose:
;	derivación de las líneas base de variación diurna o Bsq
;author:
;	Carlos Isaac Castellanos Velazco
;	Estudiante de Maestría en Ciencias de la Tierra
;	Instituto de Geofísica, Unidad Michoacan
;	UNAM
;	ccastellanos@igeofisica.unam.mx
;
;category:
;   data analysis
;
;calling sequence:
;   .r bsq_V2
;   bsq_V2, teo_mag, idate[yyyy,mm,dd], fdate[yyyy,mm,dd]
;parameters:
;
;
;dependencies:
;
;
;input files
;   teoloyucan magnetic field measurements
;
;output files:
;   figure of diurnal base line and SQ files in 24 h resolution



FUNCTION teo, date

	On_error, 2
	COMPILE_OPT idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;###############################################################################
;reading data files
;###############################################################################
        date = STRING(year, month, day, FORMAT = '(I4, I02, I02)')
		
		file_name = '../rutidl/teoloyucan/teoloyucan/'+'teo_'+date+'.clean.dat'

		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)

		data = STRARR(number_of_lines)

		OPENR, LUN, file, /GET_LUN, ERROR=err
		READF, LUN, data, FORMAT = '(A)'
		CLOSE, LUN
		FREE_LUN, LUN
;###############################################################################
;extracting data and denfining an structure data
;###############################################################################
        DStruct = {year:0, month:0, day:0, hour:0, minuntes:0, DOY:0, $ 
                     TEOD:0., TEOH:0., TEOZ:0., TEOF:0.}

		teo_mag = REPLICATE(DStruct, number_of_lines)	
        header = 0             ; Defining number of lines of the header 

		READS, data[header:number_of_lines-1], teo_mag, $
		FORMAT='(I4,X,I02,X,I02,X,I02,X,I02,8X,I03,F12,F10,F10,F10)'

		RETURN, teo_mag		
END


PRO qd_selec, date_i, date_f
	On_error, 2
	COMPILE_OPT idl2, HIDDEN
	
	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]

        file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1
        data_file_name = strarr(file_number)
        string_date     = strarr(file_number)
       
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                data_file_name[i] = '../rutidl/teoloyucan/teoloyucan/'+'teo_'+string_date[i]+'.clean.dat'
        ENDFOR

        exist_data_file   = FILE_TEST(data_file_name)
        capable_to_plot   = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_plot NE N_ELEMENTS(data_file_name) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        H    = FLTARR(file_number*1440)                               
        FOR i = 0, N_ELEMENTS(exist_data_file)-1 DO BEGIN
                IF exist_data_file[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        dat = teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*1440:(i+1)*1440-1] = dat.TEOH[*]                                                
                ENDIF ELSE BEGIN
                         H[i*1440:(i+1)*1440-1] = 999999.0
                ENDELSE                
        ENDFOR

        i_nan1 = WHERE(H EQ 999999.00, ncount)
        i_nan2 = WHERE(H EQ 99999.0, n2count) 
       
        prcent_nan = FLOAT(ncount+n2count)*100.0     
        PRINT,'porcentaje de valores NaN:', prcent_nan/N_ELEMENTS(H),'%'       
        PRINT, 'cantidad de valores NaN: ', N_ELEMENTS(i_nan1)+N_ELEMENTS(i_nan2)
        PRINT, 'Posición de valores NaN: '
        ;PRINT, i_nan1, i_nan2
PRINT, '#######################################################################'
        FOR i=0, N_ELEMENTS(H)-1 DO BEGIN
            IF H[i] EQ 999999.0 THEN BEGIN
                H[WHERE(H[*] EQ 999999.0)] = !Values.F_NAN          
            ENDIF
        ENDFOR

        FOR i=0, N_ELEMENTS(H)-1 DO BEGIN
            IF H[i] EQ 99999.0 THEN BEGIN
                H[WHERE(H[*] EQ 99999.0)] = !Values.F_NAN          
            ENDIF
        ENDFOR        
;Cálculo de las desviaciones estándar en formato horario.
H_STDESV = FINDGEN(N_ELEMENTS(H)/60)  
    	
FOR i=0, N_ELEMENTS(H_STDESV)-1 DO BEGIN 
    ;PRINT, STDDEV(H[i*60:(i+1)*60-1], /NAN) 
    H_STDESV[i] = STDDEV(H[i*60:(i+1)*60-1], /NAN)
ENDFOR   	
;Determinación de las máximas desviaciones estándar por cada día
H_std_max = FINDGEN(N_ELEMENTS(H_STDESV)/24)
FOR i=0, file_number-1 DO BEGIN
    tmp_year    = 0
    tmp_month   = 0
    tmp_day     = 0
    tmp_julday  = JULDAY(mh_i, dy_i, yr_i)
    CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
    string_date[i]    = STRING(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')        
 ;  PRINT, string_date[i], max(H_hr[i*24:(i+1)*24-1]), $
  ; FORMAT='(8A, 4X, I03)'   
    ENDFOR 	    
PRINT, '#######################################################################'
PRINT, '#######################################################################'       
tw = FINDGEN(file_number*24)/24.    
;###############################################################################   
;print figures
    time = FINDGEN(file_number *24)/24.0

    path = '../rutidl/output/eventos_tgm/'

        Device_bak = !D.Name 
        SET_PLOT, 'Z'

        
        Xsize=fix(1600)
        Ysize=250
        DEVICE, SET_RESOLUTION = [Xsize,Ysize]
        DEVICE, z_buffering=O
        DEVICE, set_character_size = [10, 12]
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;definición de color
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-        
        chr_size1 = 0.9
        chr_thick1= 1.0
        space     = 0.015
        rojo      = 248
        naranja   = 220
        amarillo  = 198
        verde     = 160
        negro     = 0
        gris_o    = 100
        blanco    = 255
        gris      = 130
        morado    = 248
                     
    TVLCT, R_bak, G_bak, B_bak, /GET        
    LOADCT, 39, /SILENT

     X_label   = STRARR(file_number+1)+' '
        months    = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        old_month = mh_i
        ;print, old_month
        FOR i =0,  N_ELEMENTS(X_label)-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                
                IF i LT N_ELEMENTS(X_label)-1 THEN $
                        X_label[i]  = (tmp_month EQ old_month) ? string(tmp_day, FORMAT='(I02)') : string(tmp_day, FORMAT='(I02)')
                old_month = tmp_month
        ENDFOR        
;##############################################################################
;Declaración de fechas
;##############################################################################
;mes del primer día quieto
mh_i2 = mh_i
CASE mh_i2 of
    1: mh_i2 = 'Enero'
    2: mh_i2 ='Febrero'
    3: mh_i2 ='Marzo'
    4: mh_i2 ='Abril'
    5: mh_i2 ='Mayo'
    6: mh_i2 ='Junio'
    7: mh_i2 ='Julio'
    8: mh_i2 ='Agosto'
    9: mh_i2 ='Septiembre'
    10:mh_i2 ='Octubre'
    11:mh_i2 ='Noviembre'
    12:mh_i2 ='Diciembre'
    ELSE: PRINT, 'fuera de rango'
ENDCASE
;###############################################################################
;mes del segundo día quieto
mh_f2 = mh_f
CASE mh_f2 of
    1: mh_f2 = 'Enero'
    2: mh_f2 ='Febrero'
    3: mh_f2 ='Marzo'
    4: mh_f2 ='Abril'
    5: mh_f2 ='Mayo'
    6: mh_f2 ='Junio'
    7: mh_f2 ='Julio'
    8: mh_f2 ='Agosto'
    9: mh_f2 ='Septiembre'
    10:mh_f2 ='Octubre'
    11:mh_f2 ='Noviembre'
    12:mh_f2 ='Diciembre'
    ELSE: PRINT, 'fuera de rango'
ENDCASE   
;###############################################################################
;fecha del primer dia quieto seleccionado
dy_i2 = dy_i
idate1 = string(yr_i, mh_i, format='(I4,I02)')
dq1 = idate1
case dq1 of
    '200310' : dq1 = 6
    '200411' : dq1 = 5
    '200505' : dq1 = 4
    '201503' : dq1 = 6
    '201705' : dq1 = 6
    '201708' : dq1 = 18
    else: print, 'fuera de rango'
endcase  

IF idate1 EQ '200505' THEN BEGIN
        for i=0, n_elements(H_STDESV)-1 do begin
            if H_STDESV[i] gt 1000.0 then begin
                H_STDESV[where(H_STDESV[*] gt 1000.0)] = !Values.F_NAN          
            endif
        endfor        
ENDIF

IF idate1 EQ '201705' THEN BEGIN
        for i=0, n_elements(H_STDESV)-1 do begin
            if H_STDESV[i] gt 40.0 then begin
       ;         H_STDESV[where(H_STDESV[*] gt 40.0)] = !Values.F_NAN          
            endif
        endfor        
ENDIF

IF idate1 EQ '201708' THEN BEGIN
        for i=0, n_elements(H_STDESV)-1 do begin
            if H_STDESV[i] GT 100.0 then begin
                H_STDESV[where(H_STDESV[*] GT 100.0)] = !Values.F_NAN          
            endif
        endfor        
ENDIF
;###############################################################################
;fecha del n1 dia quieto 
dy_i01 = dy_i
fdate1 = string(yr_i, mh_f, format='(I4,I02)')
dq01 = idate1
case dq01 of
    '200310' : dq01 = 5
    '200411' : dq01 = 1
    '200505' : dq01 = 99999
    '201503' : dq01 = file_number-5    
    '201705' : dq01 = 5
    '201708' : dq01 = 5
    else: print, 'fuera de rango'
endcase                   

;###############################################################################
;fecha del n2 dia quieto 
dy_i02 = dy_i
fdate1 = string(yr_i, mh_f, format='(I4,I02)')
dq02 = idate1
case dq02 of
    '200310' : dq02 = 7
    '200411' : dq02 = 99999
    '200505' : dq02 = 99999
    '201503' : dq02 = 22
    '201705' : dq02 = 4
    '201708' : dq02 = 16
    else: print, 'fuera de rango'
endcase    
;###############################################################################
;fecha del n3 dia quieto 
dy_i02 = dy_i
fdate1 = string(yr_i, mh_f, format='(I4,I02)')
dq03 = idate1
case dq03 of
    '200310' : dq03 = 18
    '200411' : dq03 = 99999
    '200505' : dq03 = 99999
    '201503' : dq03 = file_number-6
    '201706' : dq03 = 19
    '201708' : dq03 = 18
    else: print, 'fuera de rango'
endcase                 
;###############################################################################
;día del segundo día quieto seleccionado
dy_f2 = dy_f
fdate1 = string(yr_i, mh_f, format='(I4,I02)')
dq2 = fdate1
case dq2 of
    '200311' : dq2 = file_number-3
    '200411' : dq2 = file_number-7
    '200505' : dq2 = file_number-7
    '201504' : dq2 = file_number-5
    '201706' : dq2 = 15
    '201709' : dq2 = file_number-9
    else: print, 'fuera de rango'
endcase 
;###############################################################################
;día del segundo día n1
dy_f2 = dy_f
fdate1 = string(yr_i, mh_f, format='(I4,I02)')
dq04 = fdate1
case dq04 of
    '200311' : dq04 = file_number-4
    '200411' : dq04 = file_number-4
    '200505' : dq04 = file_number-8
    '201504' : dq04 = 31
    '201706' : dq04 = 15
    '201709' : dq04 = file_number-5
    else: print, 'fuera de rango'
endcase   
;###############################################################################
;día del segundo día n2
dy_f2 = dy_f
fdate1 = string(yr_i, mh_f, format='(I4,I02)')
dq05 = fdate1
case dq05 of
    '200311' : dq05 = file_number-2
    '200411' : dq05 = 99999
    '200505' : dq05 = file_number-4
    '201504' : dq05 = 35
    '201706' : dq05 = file_number-11
    '201709' : dq05 = file_number-8
    else: print, 'fuera de rango'
endcase 
;###############################################################################
;día del segundo día n3
dy_f2 = dy_f
fdate1 = string(yr_i, mh_f, format='(I4,I02)')
dq06 = fdate1
case dq06 of
    '200311' : dq06 = 99999
    '200411' : dq06 = 99999
    '200505' : dq06 = file_number-6
    '201504' : dq06 = file_number-2
    '201706' : dq06 = file_number-1
    '201709' : dq06 = file_number-6
    else: print, 'fuera de rango'
endcase  
;###############################################################################
;fecha en que inició la respectiva TGM 
idate0 = string(yr_i, mh_i, format='(I4,I02)')
TGM_i = idate0
case TGM_i of
    '200310' : TGM_i = 23
    '200411' : TGM_i = 6
    '200505' : TGM_i = 14
    '201503' : TGM_i = 8
    '201705' : TGM_i = 11
    '201708' : TGM_i = file_number-25
    else: print, 'fuera de rango'
endcase  
;###############################################################################
;fecha en que terminó la respectiva TGM 
fdate0 = string(yr_i, mh_f, format='(I4,I02)')
TGM_f = fdate0
case TGM_f of
    '200311' : TGM_f = 27
    '200411' : TGM_f = 12
    '200505' : TGM_f = 16
    '201504' : TGM_f = 16
    '201706' : TGM_f = 7
    '201709' : TGM_f = file_number-20
    else: print, 'fuera de rango'
endcase                   
;###############################################################################
       days = intarr(file_number+1)
       for i=0, n_elements(days)-1 do begin
            days[i] = dy_i+i
       endfor
       tot_days = days*24/24. 
       day_time = findgen(24)                       
;###############################################################################                 
IF mh_i EQ mh_f THEN BEGIN
    time_title = STRING(mh_i2, yr_i, FORMAT='((A, X, I4))')
ENDIF ELSE BEGIN
    time_title = STRING(mh_i2, FORMAT='((A))')+' y '+$
    STRING(mh_f2, ',', yr_i, FORMAT='((A,A,X,I4))')   
ENDELSE
;###############################################################################
FOR i=0, file_number-1 DO BEGIN
    tmp_year    = 0
    tmp_month   = 0
    tmp_day     = 0
    tmp_julday  = JULDAY(mh_i, dy_i, yr_i)
    CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
    string_date[i]    = STRING(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')        

   ; OPENW, LUN, outfile[i], /GET_LUN        
  ; PRINT, string_date[i], max(H_STDESV[i*24:(i+1)*24-1]), $
   FORMAT='(8A, 4X, I03)' 
    
   ; CLOSE, LUN
   ; FREE_LUN, LUN    
    ENDFOR   
;###############################################################################
    sigma = TeXtoIDL('\sigmaH')
    dH    = TeXtoIDL('\DeltaH')
    time_title = 'Desviacion estandar de '+dH+' para '+ time_title
;###############################################################################    
     down0 = min(H_STDESV)   
     up0   = max(H_STDESV)  

    plot, time, H_STDESV, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=rojo,$
     CHARSIZE = 1.0, CHARTHICK=chr_thick1, POSITION=[0.05,0.15,0.95,0.90], $
     XSTYLE = 5, XTICKNAME=REPLICATE(' ', file_number+1), XRANGE=[0, file_number], $
     ySTYLE = 6, YRANGE=[down0,up0], THICK=4         
;###############################################################################    
    POLYFILL, [dq01, dq01+1, dq01+1, dq01],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3, SPACING=0.3             

    POLYFILL, [dq01, dq01+1, dq01+1, dq01],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=-45, LINESTYLE=0, THICK=3, SPACING=0.3   
                            
    POLYFILL, [dq02, dq02+1, dq02+1, dq02],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3, SPACING=0.3   

    POLYFILL, [dq02, dq02+1, dq02+1, dq02],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=-45, LINESTYLE=0, THICK=3, SPACING=0.3  
              
    POLYFILL, [dq03, dq03+1, dq03+1, dq03],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3, SPACING=0.3   

    POLYFILL, [dq03, dq03+1, dq03+1, dq03],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=-45, LINESTYLE=0, THICK=3, SPACING=0.3       
              
    POLYFILL, [dq04, dq04+1, dq04+1, dq04],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3, SPACING=0.3   

    POLYFILL, [dq04, dq04+1, dq04+1, dq04],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=-45, LINESTYLE=0, THICK=3, SPACING=0.3      

    POLYFILL, [dq05, dq05+1, dq05+1, dq05],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3, SPACING=0.3   

    POLYFILL, [dq05, dq05+1, dq05+1, dq05],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=-45, LINESTYLE=0, THICK=3, SPACING=0.3      
              
    POLYFILL, [dq06, dq06+1, dq06+1, dq06],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3, SPACING=0.3   

    POLYFILL, [dq06, dq06+1, dq06+1, dq06],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=-45, LINESTYLE=0, THICK=3, SPACING=0.3                                                                   
;###############################################################################    
    ;días quietos seleccionados
    POLYFILL, [TGM_i, TGM_f, TGM_f, TGM_i], $
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=amarillo
              
    POLYFILL, [dq1, dq1+1, dq1+1, dq1],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde;, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3           
              
    POLYFILL, [dq2, dq2+1, dq2+1, dq2],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde;, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3                              

IF idate0 EQ '200310' THEN BEGIN
    POLYFILL, [TGM_i+22, TGM_f+22, TGM_f+22, TGM_i+22], $
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=amarillo   
    OPLOT, [tgm_i+22, tgm_i+22], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5,$
 color=negro, THICK=3 
    OPLOT, [tgm_f+22, tgm_f+22], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5, $
    color=negro, THICK=3                 
ENDIF             
;###############################################################################           
;días quietos    
    OPLOT, [dq01,dq01], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2    
    OPLOT, [dq01+1,dq01+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2  

    OPLOT, [dq02,dq02], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2    
    OPLOT, [dq02+1,dq02+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2 
    
    OPLOT, [dq03,dq03], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2    
    OPLOT, [dq03+1,dq03+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2 
    
    OPLOT, [dq04,dq04], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2    
    OPLOT, [dq04+1,dq04+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2 
    
    OPLOT, [dq05,dq05], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2    
    OPLOT, [dq05+1,dq05+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2         
 
    OPLOT, [dq06,dq06], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2    
    OPLOT, [dq06+1,dq06+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2        
;###############################################################################       
;días quietos seleccionados              
    OPLOT, [dq1,dq1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=3    
    OPLOT, [dq1+1,dq1+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=3   
        
    OPLOT, [dq2,dq2], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=3
    OPLOT, [dq2+1,dq2+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=3   
    
    OPLOT, [tgm_i, tgm_i], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5, color=negro, THICK=3 
    OPLOT, [tgm_f, tgm_f], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5, color=negro, THICK=3  
    
    OPLOT, time, H_STDESV, LINESTYLE=0, THICK=4 , COLOR=rojo                                                              
;###############################################################################       
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 1.0, $
                         xtitle='',$
                         CHARTHICK=2
                         
        AXIS, XAXIS = 1, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', file_number+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down0,up0],$
                         Ystyle=2, $
                         COLOR=negro, $
                         CHARSIZE = 1.0, $
                         CHARTHICK=1.5;, $

        AXIS, YAXIS = 1, YRANGE=[down0,up0],$
                         Ystyle=2, $
                         COLOR=negro, $
                         CHARSIZE = 1.0, $
                         CHARTHICK=1.5;, $
;###############################################################################                         
   y = (!Y.Window[1] - !Y.Window[0]) / 2. + !Y.Window[0] 
   XYOUTS, 0.02, y, sigma+' [nT]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.2, ORIENTATION=90, CHARTHICK=1.5                             
;###############################################################################
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.92   
   XYOUTS, X, y, time_title, /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.65, CHARTHICK=2 
   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.02   
   XYOUTS, X, y, 'Tiempo universal [dias]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.0, CHARTHICK=1.0
;###############################################################################    
IF idate0 EQ '200310' THEN BEGIN
;first panel legend                   

       
        POLYFILL, [0.07,0.1,0.1,0.07], [0.832,0.832,0.872,0.872], color = verde, $
        /NORMAL, /LINE_FILL, ORIENTATION=45, THICK=3, LINESTYLE=0, SPACING=0.3
        
        POLYFILL, [0.07,0.1,0.1,0.07], [0.832,0.832,0.862,0.862], color = verde, $
        /NORMAL, /LINE_FILL, ORIENTATION=-45, THICK=3, LINESTYLE=0, SPACING=0.3 

        POLYFILL, [0.07,0.1,0.1,0.07], [0.752,0.752,0.802,0.802], color = verde, $
        /NORMAL         
        
        
        POLYFILL, [0.07,0.1,0.1,0.07], [0.672,0.672,0.722,0.722], color = amarillo, $
        /NORMAL
        

                             
        XYOUTS, 0.105, 0.83 , /NORMAL, $
                'DQ', COLOR=negro, $
                CHARSIZE = 1.1, $
                CHARTHICK=2
                
        XYOUTS, 0.105, 0.75, /NORMAL, $
                'DQS', COLOR=negro, $
                CHARSIZE = 1.1, $
                CHARTHICK=2                
                
        XYOUTS,  0.105, 0.67 , /NORMAL, $
                'TGM', COLOR=negro, $
                CHARSIZE = 1.1, $
                CHARTHICK=2                                
ENDIF
;###############################################################################                     
    Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!
    
    TVLCT, R_bak, G_bak, B_bak
        
    SET_PLOT, Device_bak
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; open the post stript device
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')
        IF keyword_set(jpeg) THEN BEGIN
                info = size(Image)
                nx = info[1]
                ny = info[2]
                true_image = bytarr(3,nx,ny)
                true_image[0,*,*] = reds[image]
                true_image[1,*,*] = greens[image]
                true_image[2,*,*] = blues[image]
                write_jpeg, path+'DH_stdesv_V2_'+fecha+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'DH_stdesv_V2_'+fecha+'.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'DH_stdesv_V2_'+fecha
                print, ''
        ENDIF
        RETURN   
	
END	

PRO bsq_v2, teo_mag, date_i, date_f
	On_error, 2
	COMPILE_OPT idl2, HIDDEN
	
	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]
;###############################################################################	
	dat1    = teo([yr_i, mh_i, dy_i])	
	H1      = dat1.TEOH

    dat2    = teo([yr_f, mh_f, dy_f])
    H2      = dat2.TEOH
    
    ndays   = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1
    td      = FINDGEN(ndays*1440)/1440.0
    td_h    = FINDGEN(ndays*24)/24.0
    datetime= TIMEGEN(N_ELEMENTS(td), FINAL=JULDAY(mh_f, dy_f, yr_f, 23), $
                START=JULDAY(mh_i, dy_i, yr_i, 0), UNITS='H')
    CALDAT, datetime, mh, dy, yr, hr
;###############################################################################                        
        i_nan1 = WHERE(H1 EQ 999999.00, ncount)
        
        prcent_nan = FLOAT(ncount)*100.0
        PRINT,'porcentaje de valores NaN:', prcent_nan/N_ELEMENTS(H1),'%'
        
        i_nan2 = WHERE(H2 EQ 999999.00, ncount2)
        
        prcent_nan = FLOAT(ncount2)*100.0
        PRINT,'porcentaje de valores NaN:', prcent_nan/N_ELEMENTS(H2),'%'
        
        for i=0, n_elements(H1)-1 do begin
            if H1[i] eq 999999.0 then begin
                H1[where(H1[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor
        
        for i=0, n_elements(H2)-1 do begin
            if H2[i] eq 999999.0 then begin
                H2[where(H2[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor  
        
        FOR i=0, N_ELEMENTS(H1)-1 DO BEGIN
            IF H1[i] EQ 99999.0 THEN BEGIN
                H1[WHERE(H1[*] EQ 99999.0)] = !Values.F_NAN          
            ENDIF
        ENDFOR

        FOR i=0, N_ELEMENTS(H2)-1 DO BEGIN
            IF H2[i] EQ 99999.0 THEN BEGIN
                H2[WHERE(H2[*] EQ 99999.0)] = !Values.F_NAN          
            ENDIF
        ENDFOR           
;###############################################################################        
    ;implementar una función de interpolación en caso de que el porcentaje de 
    ;nan sea muy bajo       
    H1_tmp   = H1
    H1_exist = WHERE(finite(H1_tmp), ngooddata1, complement=baddata1, $
    ncomplement=nbaddata1)

    H2_tmp   = H2
    H2_exist = WHERE(finite(H2_tmp), ngooddata2, complement=baddata2, $
    ncomplement=nbaddata2)
       
    ; interpolate at the locations of the bad data using the good data    
    IF nbaddata1 GT 0 THEN H1_tmp[baddata1] = INTERPOL(H1_tmp[H1_exist], $
    H1_exist, baddata1, /QUADRATIC)
    H1 = H1_tmp  
    
    IF nbaddata1 GT 0 THEN H2_tmp[baddata2] = INTERPOL(H2_tmp[H2_exist], $
    H2_exist, baddata2, /QUADRATIC)
    H2 = H2_tmp 
;###############################################################################
;Extend QDS data ndays for a quadratic interpolation
QDS1 = REFORM(REBIN(H1, 1440, ndays), N_ELEMENTS(td))
QDS2 = REFORM(REBIN(H2, 1440, ndays), N_ELEMENTS(td))
;###############################################################################
;Interpolate between the QDS
slope1   = (td - td[719])
slope2   = (td[N_ELEMENTS(td)-721]-td[719])
slope    = slope1/slope2

Bsq1     = (QDS2-QDS1)*slope
Bsq     = QDS1+Bsq1                           
;###############################################################################
    ;Generate a QD file in days
    outfile = STRARR(ndays) 
    
    string_date        = STRARR(ndays)                       
    data_file_name_h  = STRARR(ndays)

    Bsq_H = FINDGEN(N_ELEMENTS(Bsq)/60)
    FOR i=0, N_ELEMENTS(Bsq_H)-1 DO BEGIN
        ;print, Bsq[i*60:(i+1)*60-1]
        Bsq_H[i] = MEDIAN(Bsq[i*60:(i+1)*60-1])        
    
    ENDFOR

;Generación de archivo en muestreo de horas     
FOR i=0, ndays-1 DO BEGIN
    tmp_year    = 0
    tmp_month   = 0
    tmp_day     = 0
    tmp_julday  = JULDAY(mh_i, dy_i, yr_i)
    CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
    string_date[i]    = STRING(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')        

    outfile[i] = '../rutidl/output/Bsq_baselines/Bsq_'+string_date[i]+'h.dat'    
    OPENW, LUN, outfile[i], /GET_LUN        
    PRINTF, LUN, Bsq_H[i*24:(i+1)*24-1], format='(F10.4)'
    CLOSE, LUN
    FREE_LUN, LUN    
ENDFOR 

;archivo en muestreo de minutos
FOR i=0, ndays-1 DO BEGIN
    tmp_year    = 0
    tmp_month   = 0
    tmp_day     = 0
    tmp_julday  = JULDAY(mh_i, dy_i, yr_i)
    CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
    string_date[i]    = STRING(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')        

    outfile[i] = '../rutidl/output/Bsq_baselines/Bsq_'+string_date[i]+'m.dat'    
    OPENW, LUN, outfile[i], /GET_LUN        
    PRINTF, LUN, Bsq[i*1440:(i+1)*1440-1], format='(F10.4)'
    CLOSE, LUN
    FREE_LUN, LUN    
ENDFOR    
    
;###############################################################################
;Generación de la figura de Bsq
path = '../rutidl/output/eventos_tgm/'
Date = string(yr_i, mh_i, FORMAT='(I4, "-", I02)')

        Device_bak = !D.Name 
        SET_PLOT, 'Z'    
        
        Xsize=FIX(1600)
        Ysize=300
        ;DEVICE, decompose=0
        DEVICE, SET_RESOLUTION = [Xsize,Ysize]
        DEVICE, z_buffering=O
        DEVICE, set_character_size = [10, 12]

    TVLCT, R_bak, G_bak, B_bak, /GET        
    LOADCT, 39, /SILENT
             
        chr_size1 = 0.9
        chr_thick1= 1.0
        space     = 0.015
        rojo      = 248
        amarillo  = 190
        verde     = 170
        negro     = 0
        azul      = 70
        blanco    = 255
        gris      = 130
        morado    = 16
        
    TVLCT, R_bak, G_bak, B_bak, /GET   
    LOADCT, 39, /SILENT
;###############################################################################    
     X_label   = STRARR(ndays+1)+' '
        months    = ['Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', $
                     'Sep', 'Oct', 'Nov', 'Dic']
        old_month = mh_i
        FOR i =0,  N_ELEMENTS(X_label)-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                
                IF i LT N_ELEMENTS(X_label)-1 THEN $
                        X_label[i]  = (tmp_month EQ old_month) ? STRING(tmp_day, $
                        FORMAT='(I02)') : STRING(tmp_day, FORMAT='(I02)')
                        old_month = tmp_month
        ENDFOR 
;##############################################################################
;Declaración de fechas
;##############################################################################
;mes del primer día quieto
mh_i2 = mh_i
CASE mh_i2 of
    1: mh_i2 = 'Enero'
    2: mh_i2 ='Febrero'
    3: mh_i2 ='Marzo'
    4: mh_i2 ='Abril'
    5: mh_i2 ='Mayo'
    6: mh_i2 ='Junio'
    7: mh_i2 ='Julio'
    8: mh_i2 ='Agosto'
    9: mh_i2 ='Septiembre'
    10:mh_i2 ='Octubre'
    11:mh_i2 ='Noviembre'
    12:mh_i2 ='Diciembre'
    ELSE: PRINT, 'fuera de rango'
ENDCASE
;###############################################################################
;mes del segundo día quieto
mh_f2 = mh_f
CASE mh_f2 of
    1: mh_f2 = 'Enero'
    2: mh_f2 ='Febrero'
    3: mh_f2 ='Marzo'
    4: mh_f2 ='Abril'
    5: mh_f2 ='Mayo'
    6: mh_f2 ='Junio'
    7: mh_f2 ='Julio'
    8: mh_f2 ='Agosto'
    9: mh_f2 ='Septiembre'
    10:mh_f2 ='Octubre'
    11:mh_f2 ='Noviembre'
    12:mh_f2 ='Diciembre'
    ELSE: PRINT, 'fuera de rango'
ENDCASE   
;###############################################################################
;fecha en que inició la respectiva TGM 
idate0 = string(yr_i, mh_i, format='(I4,I02)')
TGM_i = idate0
case TGM_i of
    '200310' : TGM_i = td_h[360]
    '200411' : TGM_i = td_h[24]
    '200505' : TGM_i = td_h[240]
    '201503' : TGM_i = td_h[144]
    '201705' : TGM_i = td_h[24]
    '201708' : TGM_i = td_h[216]
    else: print, 'fuera de rango'
endcase  
;###############################################################################
;fecha en que terminó la respectiva TGM 
fdate0 = string(yr_i, mh_f, format='(I4,I02)')
TGM_f = fdate0
case TGM_f of
    '200311' : TGM_f = td_h[456]
    '200411' : TGM_f = td_h[144]
    '200505' : TGM_f = td_h[288]
    '201504' : TGM_f = td_h[216]
    '201706' : TGM_f = td_h[72]
    '201709' : TGM_f = td_h[288]
    else: print, 'fuera de rango'
endcase                   
;###############################################################################
       days = intarr(ndays+1)
       for i=0, n_elements(days)-1 do begin
            days[i] = dy_i+i
       endfor
       tot_days = days*24/24. 
       day_time = findgen(24)                       
;###############################################################################                 
IF mh_i EQ mh_f THEN BEGIN
    time_title = STRING(mh_i2, yr_i, FORMAT='((A, X, I4))')
ENDIF ELSE BEGIN
    time_title = STRING(mh_i2, FORMAT='((A))')+' y '+$
    STRING(mh_f2, ',', yr_i, FORMAT='((A,A,X,I4))')   
ENDELSE
;###############################################################################
;ylimits
IF MAX(Bsq_H) GT MAX(QDS1) THEN up = MAX(Bsq_H) ELSE up = MAX(QDS1)
IF MIN(Bsq_H) LT MIN(QDS1) THEN down=MIN(Bsq_H) ELSE down=MIN(QDS1)
;###############################################################################
;xlimits
ini = td[0]
fin = td[N_ELEMENTS(td)-1]
;###############################################################################                                 
    PLOT, td_h, Bsq_H, position = [0.09, 0.14, 0.91, 0.9], XRANGE=[ini,fin], $
    YRANGE=[down,up], YSTYLE=6, XSTYLE=5, BACKGROUND = blanco, COLOR=negro,$
    XTICKNAME=REPLICATE(' ', ndays+1)

    POLYFILL, [TGM_i, TGM_f, TGM_f, TGM_i], $
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=amarillo
              
    POLYFILL, [td_h[0], td_h[23], td_h[23], td_h[0]],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3, $
              SPACING=0.3  
              
    POLYFILL, [td_h[0], td_h[23], td_h[23], td_h[0]],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=-45, LINESTYLE=0, THICK=3, $
              SPACING=0.3                                  

    POLYFILL, [td_h[N_ELEMENTS(td_h)-1], td_h[N_ELEMENTS(td_h)-24], $
               td_h[N_ELEMENTS(td_h)-25], td_h[N_ELEMENTS(td_h)-1]],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=45, LINESTYLE=0, THICK=3, $
              SPACING=0.3     
              
    POLYFILL, [td_h[N_ELEMENTS(td_h)-1], td_h[N_ELEMENTS(td_h)-24], $
               td_h[N_ELEMENTS(td_h)-25], td_h[N_ELEMENTS(td_h)-1]],$
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=verde, /LINE_FILL, ORIENTATION=-45, LINESTYLE=0, THICK=3, $
              SPACING=0.3        
              
              
    OPLOT, [td_h[23], td_h[23]], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, $
    color=negro, THICK=2                
    OPLOT, [td_h[N_ELEMENTS(td_h)-25], td_h[N_ELEMENTS(td_h)-25]], $
    [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=2
              
    OPLOT, [TGM_i, TGM_i], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5, $
    color=negro, THICK=2                   
    OPLOT, [TGM_f, TGM_f], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5, $
    color=negro, THICK=2     
              
IF idate0 EQ '200310' THEN BEGIN
    POLYFILL, [TGM_i+22, TGM_f+22, TGM_f+22, TGM_i+22], $
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=amarillo   
    OPLOT, [tgm_i+22, tgm_i+22], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5,$
 color=negro, THICK=2 
    OPLOT, [tgm_f+22, tgm_f+22], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5, $
    color=negro, THICK=2                
ENDIF                          
    OPLOT, td_h, Bsq_H, COLOR=negro, THICK=4   
;###############################################################################
        AXIS, XAXIS = 0, XRANGE=[0,ndays], $
                         XTICKS=ndays, $                     
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 1.1 , $
                         CHARTHICK=1.5,$
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,ndays], $     
                         XTICKFORMAT="(A1)", $                                            
                         XMINOR=8, $ 
                         CHARSIZE = 1.1 , $ 
                         CHARTHICK=1.5,$                      
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down,up], $                     
                         COLOR=negro, $
                         CHARTHICK=1.5,$
                         YSTYLE=2, $
                         CHARSIZE = 1.1 ;, $
                        
        AXIS, YAXIS = 1, YRANGE=[down,up], $
                         COLOR=negro, $
                         CHARTHICK=1.5,$                         
                         YSTYLE=2, $
                         CHARSIZE = 1.1 ;, $  
;###############################################################################                         
   y = (!Y.Window[1] - !Y.Window[0]) / 2. + !Y.Window[0] 
   XYOUTS, 0.02, y, '[nT]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.2, ORIENTATION=90, CHARTHICK=1.2                            
;###############################################################################
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.925   
   XYOUTS, X, y, time_title, /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.75, CHARTHICK=1.2 
;###############################################################################  
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.01, 'Tiempo Universal en dias', /Normal, $
   color=negro, Alignment=0.5, Charsize=1, CHARTHICK=1.2   
;###############################################################################
IF idate0 EQ '200310' THEN BEGIN
;first panel legend                          
        POLYFILL, [0.80,0.82,0.82,0.80], [0.832,0.832,0.862,0.862], color = verde, $
        /NORMAL, /LINE_FILL, ORIENTATION=45, THICK=3, LINESTYLE=0, SPACING=0.3   
        
        POLYFILL, [0.80,0.82,0.82,0.80], [0.832,0.832,0.862,0.862], color = verde, $
        /NORMAL, /LINE_FILL, ORIENTATION=-45, THICK=3, LINESTYLE=0, SPACING=0.3  

        POLYFILL, [0.80,0.82,0.82,0.80], [0.764,0.764,0.794,0.794], color = amarillo, /NORMAL         
                             
        XYOUTS, 0.83, 0.832 , /NORMAL, $
                'DGS', COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=2
                
        XYOUTS, 0.83, 0.76 , /NORMAL, $
                'TGM', COLOR=negro, $
                CHARSIZE = 1.2, $
                CHARTHICK=2                                
ENDIF                   
;###############################################################################                
     Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!
    
    TVLCT, R_bak, G_bak, B_bak
        
    ;DEVICE, /CLOSE
    SET_PLOT, Device_bak   
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; open the post stript device
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        IF keyword_set(jpeg) THEN BEGIN
                info = size(Image)
                nx = info[1]
                ny = info[2]
                true_image = BYTARR(3,nx,ny)
                true_image[0,*,*] = reds[image]
                true_image[1,*,*] = greens[image]
                true_image[2,*,*] = blues[image]
                write_jpeg, path+'Bsq_'+Date+'_V3.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'Bsq_'+Date+'_V3.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'Bsq_'+Date+'_V3.png'
                print, ''
        ENDIF
        RETURN         
    
    
END

PRO H_hour, date_i, date_f
	On_error, 2
	COMPILE_OPT idl2, HIDDEN
	
	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]

        file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1
        data_file_name = strarr(file_number)
        string_date     = strarr(file_number)
       
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                data_file_name[i] = '../rutidl/teoloyucan/teoloyucan/'+'teo_'+string_date[i]+'.clean.dat'
        ENDFOR

        exist_data_file   = FILE_TEST(data_file_name)
        capable_to_plot   = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_plot NE N_ELEMENTS(data_file_name) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        H    = FLTARR(file_number*1440)                               
        FOR i = 0, N_ELEMENTS(exist_data_file)-1 DO BEGIN
                IF exist_data_file[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        dat = teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*1440:(i+1)*1440-1] = dat.TEOH[*]                                                
                ENDIF ELSE BEGIN
                         H[i*1440:(i+1)*1440-1] = 999999.0
                         H_STDESV[i*1440:(i+1)*1440-1] = 999999.0                        
                ENDELSE                
        ENDFOR

PRINT, '#######################################################################'
        FOR i=0, N_ELEMENTS(H)-1 DO BEGIN
            IF H[i] EQ 999999.0 THEN BEGIN
                H[WHERE(H[*] EQ 999999.0)] = !Values.F_NAN          
            ENDIF
        ENDFOR

        FOR i=0, N_ELEMENTS(H)-1 DO BEGIN
            IF H[i] EQ 99999.0 THEN BEGIN
                H[WHERE(H[*] EQ 99999.0)] = !Values.F_NAN          
            ENDIF
        ENDFOR             

    H_hr = FINDGEN(N_ELEMENTS(H)/60)
    FOR i=0, N_ELEMENTS(H_hr)-1 DO BEGIN
        H_hr[i] = MEAN(H[i*60:(i+1)*60-1])        
    
    ENDFOR

    outfile = STRARR(file_number)
;Generación de archivo en muestreo de horas     
FOR i=0, file_number-1 DO BEGIN
    tmp_year    = 0
    tmp_month   = 0
    tmp_day     = 0
    tmp_julday  = JULDAY(mh_i, dy_i, yr_i)
    CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
    string_date[i]    = STRING(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')        

    outfile[i] = '../rutidl/teoloyucan/hourly/teo_'+string_date[i]+'h.dat'    
    OPENW, LUN, outfile[i], /GET_LUN        
    PRINTF, LUN, H_hr[i*24:(i+1)*24-1], format='(F10.4)'
    CLOSE, LUN
    FREE_LUN, LUN    
ENDFOR    

END








