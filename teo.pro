;
;Name:
;	get_data_date
;purpose:
;	this routine will call read a certain number of files containing 
;   interplanetary data measurements
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
;
;
;parameters:
;
;
;dependencies:
;
;
;input files
;
;
;output files:
;
;

function DH_teo, date

	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        date = string(year, month, day, format = '(I4, I02, I02)')		
		file_name = '../rutidl/dH_teo/'+'teo_'+date+'.dst.early'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        DStruct = {hora : 0, D_stdesv : 0., D : 0., H_stdesv : 0., H : 0., $
        Z_stdesv : 0., Z : 0., N_stdesv : 0., N : 0., F_stdesv : 0., F : 0.}

		teo_mag = REPLICATE(DStruct, number_of_lines)	
  
		READS, data[0:number_of_lines-1], teo_mag, $
		FORMAT='(I2, F10, F8, F10, F10, F10, F10, F10, F10, F10, F10)'
		
		return, teo_mag		
end

function kmex, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        date = string(year, month, day, format = '(I4, I02, I02)')		
		name = 'teo_'+date+'.index.'
		
		file_name = '../rutidl/Kmex/'+name+'final'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
	    IF opened_files NE N_ELEMENTS(file) THEN begin
	        file_name = '../rutidl/Kmex/'+name+'early'
	        file = FILE_SEARCH(file_name, COUNT=opened_files) 
    	    IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+'not found'  
	    
	    endif

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		IF err NE 0 THEN MESSAGE, 'Error opening '+file_name[0]
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;extracting data and denfining an structure data
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-        
        idx_kmex      = {k_mex      : intarr(8), $
                         k_mex_sum  : 0, $
                         a_mex      : intarr(8), $
                         a_med      : 0}
        
        struct = {x : [0, 0, 0, 0, 0, 0, 0, 0], y : 0}
        
        tmp_var = replicate(struct, 2)  
		READS, data, tmp_var, FORMAT='(I3, I4, I4, I4, I4, I4, I4, I4, I4)'
		
		idx_kmex.k_mex[*]   = tmp_var[0].x
		idx_kmex.a_mex[*]   = tmp_var[1].x
        idx_kmex.k_mex_sum  = tmp_var[0].y
        idx_kmex.a_med      = tmp_var[1].y		
		
		return, idx_kmex	
end

pro teo, date_i, date_f
	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]	
;##############################################################################
; reading data files
;##############################################################################
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
                data_file_name[i] = '../rutidl/dH_teo/'+'teo_'+string_date[i]+'.dst.early'
        ENDFOR

        exist_data_file   = FILE_TEST(data_file_name)
        capable_to_plot   = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_plot NE N_ELEMENTS(data_file_name) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

        H    = FLTARR(file_number*24)        
        H_STDESV    = FLTARR(file_number*24)
                       
        FOR i = 0, N_ELEMENTS(exist_data_file)-1 DO BEGIN
                IF exist_data_file[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        dat = DH_teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*24:(i+1)*24-1] = dat.H[*]                                                
                        H_STDESV[i*24:(i+1)*24-1] = dat.H_stdesv[*]
                                                                                              
                ENDIF ELSE BEGIN
                         H[i*24:(i+1)*24-1] = 999999.0
                         H_STDESV[i*24:(i+1)*24-1] = 999999.0                        
                ENDELSE                
        ENDFOR

        for i=0, n_elements(H)-1 do begin
            if H[i] eq 999999.0 then begin
                H[where(H[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor
        
        for i=0, n_elements(H)-1 do begin
            if H[i] ge 100.0 then begin
                H[where(H[*] ge 100.0)] = !Values.F_NAN          
            endif
        endfor

        for i=0, n_elements(H)-1 do begin
            if H_STDESV[i] eq 999999.0 then begin
                H_STDESV[where(H_STDESV[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor
   ;     for i=0, n_elements(H)-1 do begin
    ;        if H_STDESV[i] gt 40.0 then begin
     ;           H_STDESV[where(H_STDESV[*] gt 40.0)] = !Values.F_NAN          
      ;      endif
       ; endfor        
        ;print, H_STDESV
    time = findgen(file_number *24)/24.0

    path = '../rutidl/output/eventos_tgm/'

        Device_bak = !D.Name 
        SET_PLOT, 'Z'
        
        tmp_spam = 1
        IF file_number GT 7 THEN tmp_spam = 1.5
        IF file_number GT 15 THEN tmp_spam = 2.
        
        Xsize=fix(800*tmp_spam)
        Ysize=400
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
;fecha del primer dia quieto
dy_i2 = dy_i
idate1 = string(yr_i, mh_i, format='(I4,I02)')
dq1 = idate1
case dq1 of
    '200310' : dq1 = 7
    '200411' : dq1 = 5
    '200505' : dq1 = 4
    '201503' : dq1 = 2
    '201705' : dq1 = 6
    '201708' : dq1 = 6
    else: print, 'fuera de rango'
endcase  
;###############################################################################
;fecha del segundo dia quieto
dy_f2 = dy_f
fdate1 = string(yr_i, mh_f, format='(I4,I02)')
dq2 = fdate1
case dq2 of
    '200311' : dq2 = file_number-4
    '200411' : dq2 = file_number-6
    '200505' : dq2 = file_number-6
    '201504' : dq2 = file_number-4
    '201706' : dq2 = 15
    '201709' : dq2 = file_number-5
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

    OPLOT, [dq1,dq1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=3    
    OPLOT, [dq1+1,dq1+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=3   
        
    OPLOT, [dq2,dq2], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=3
    OPLOT, [dq2+1,dq2+1], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=3, color=negro, THICK=3   
    
    OPLOT, [tgm_i, tgm_i], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5, color=negro, THICK=3 
    OPLOT, [tgm_f, tgm_f], [!Y.CRANGE[0], !Y.CRANGE[1]], linestyle=5, color=negro, THICK=3     
;###############################################################################       
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         xtitle='Tiempo universal [dias]'
                         CHARTHICK=chr_thick1
                         
        AXIS, XAXIS = 1, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', file_number+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down0,up0],$
                         Ystyle=2, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1;, $

        AXIS, YAXIS = 1, YRANGE=[down0,up0],$
                         Ystyle=2, $
                         COLOR=negro, $
                         CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1;, $
;###############################################################################                         
   y = (!Y.Window[1] - !Y.Window[0]) / 2. + !Y.Window[0] 
   XYOUTS, 0.02, y, sigma+' [nT]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.2, ORIENTATION=90                           
;###############################################################################
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.92   
   XYOUTS, X, y, time_title, /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.85 
;###############################################################################                         
    Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!
    
    TVLCT, R_bak, G_bak, B_bak
        
    SET_PLOT, Device_bak
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; open the post stript device
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
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
end

pro kmx_plot, date_i, date_f
	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]    
;##############################################################################
; reading data files
;##############################################################################
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
                data_file_name[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.final'
		        
		        file = FILE_SEARCH(data_file_name[i], COUNT=opened_files)
	            IF opened_files NE N_ELEMENTS(file) THEN begin
	                data_file_name[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.early'    
	            ENDIF		        
	        
        ENDFOR
	
        exist_data_file   = FILE_TEST(data_file_name)
        capable_to_plot   = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_plot NE N_ELEMENTS(data_file_name) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

        k_mex    = INTARR(file_number*8)
        a_mex    = INTARR(file_number*8)

        FOR i = 0, N_ELEMENTS(exist_data_file)-1 DO BEGIN
                IF exist_data_file[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        dat = kmex([tmp_year, tmp_month, tmp_day])
                        
                        k_mex[i*8:(i+1)*8-1] = dat.k_mex[*]/10
                        a_mex[i*8:(i+1)*8-1] = dat.a_mex[*]
                ENDIF 
        ENDFOR
    
    time = findgen(file_number *8)/8.0     
    path = '../rutidl/output/teofig/'

        Device_bak = !D.Name 
        SET_PLOT, 'Z'
               
        Xsize=fix(1200)
        Ysize=400
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
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; Write a post Script
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-    
     X_label   = STRARR(file_number+1)+' '
        months    = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        old_month = mh_i
        FOR i =0,  N_ELEMENTS(X_label)-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                
                IF i LT N_ELEMENTS(X_label)-1 THEN $
                        X_label[i]  = (tmp_month EQ old_month) ? string(tmp_day, FORMAT='(I02)') : months[tmp_month-1]+' '+string(tmp_day, FORMAT='(I02)')
                old_month = tmp_month
        ENDFOR        
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; PLOTTING
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
time_title = 'Time window begin: '+string(yr_i, mh_i, dy_i, $
                FORMAT='(I4, "/", I02, "/", I02)')+' 00:00 UTC'
      str = TeXtoIDL('\DeltaH, K_{mex}')          
MAG_source = 'Source: International Service of Geomagnetic Indices'                
    
    plot, time, k_mex, psym=6, /NoDATA, MAX_VALUE=9., XTICKS=file_number, xminor=8, $
                    TITLE = str, SUBTITLE = time_title, XTITLE = 's', YTITLE = 's', $
                    BACKGROUND = blanco, COLOR=negro, YRANGE=[0,9], YTICKS=9, $
                    YMINOR=0, CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
                    POSITION=[0.07,0.35,0.95,0.9], XSTYLE = 5, ystyle = 5,$
                    XTICKNAME=REPLICATE(' ', file_number+1), XRANGE=[0, file_number]

        j = N_ELEMENTS(k_mex)
        FOR i = 0, j-1 DO BEGIN
                IF k_mex[i] LE 9 THEN BEGIN
                                        color = 0
                                        step  = (k_mex[i] EQ 0) ? 0.1 : 0.
                                        CASE 1 OF
                                                k_mex[i] EQ 4 : color = amarillo
                                                k_mex[i] GE 4 : color = rojo
                                                ELSE       : color = verde
                                        ENDCASE
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+time[i], [0,0,k_mex[i]$
                                        +step,k_mex[i]+step], color=color
                                      ENDIF $
                                      ELSE BEGIN
                                        LOADCT, 0, /SILENT
                                        
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+time[i], $
                                        [0,0,9.,9.], color=morado, /LINE_FILL, $
                                        ORIENTATION=45., linestyle = 0
                                         
                                         
                                        POLYFILL, [0.+space,0.125-space,0.125-$
                                        space,0.+space]+time[i], $
                                        [0,0,9.,9.],color=morado, /LINE_FILL, $
                                         ORIENTATION=-45., linestyle = 0
                                      ENDELSE
        ENDFOR
 
 FOR i = 0, file_number-1 DO BEGIN
                OPLOT, [i,i], [0.,9.], linestyle=1, COLOR=negro
        ENDFOR

        FOR i=5, 8, 1 DO OPLOT, [0,file_number], [i,i], linestyle=1, COLOR=negro

        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.7, $
                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.04
                         
        AXIS, XAXIS = 1, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', file_number+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTITLE = 'Kmex index', $
                         COLOR=negro, $
                         CHARSIZE = 0.7, $
                         CHARTHICK=chr_thick1;, $

        AXIS, YAXIS = 1, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTICKNAME=[' ', ' ', ' ', ' ', ' ', 'G1', 'G2', 'G3', 'G4', 'G5'], $
                         COLOR=negro, $
                         CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1;, $
 
        XYOUTS, 0.01, 0.070 , /NORMAL, $
                MAG_source, COLOR=negro, $
                CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1

        XYOUTS, 0.01, 0.165 , /NORMAL, $
                'Color Code:       quiet,       disturbed,         storm,        intense storm,        data not available.', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1
        POLYFILL, [0.09,0.12,0.12,0.09], [0.165,0.165,0.195,0.195], color = verde, /NORMAL
        POLYFILL, [0.17,0.20,0.20,0.17], [0.165,0.165,0.195,0.195], color = amarillo, /NORMAL
        POLYFILL, [0.29,0.32,0.32,0.29], [0.165,0.165,0.195,0.195], color = naranja, /NORMAL
        POLYFILL, [0.38,0.41,0.41,0.38], [0.165,0.165,0.195,0.195], color = rojo, /NORMAL
        POLYFILL, [0.52,0.55,0.55,0.52], [0.165,0.165,0.195,0.195], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=45., linestyle = 0
        POLYFILL, [0.52,0.55,0.55,0.52], [0.165,0.165,0.195,0.195], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=-45., linestyle = 0

    Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!
    
    TVLCT, R_bak, G_bak, B_bak
        
    SET_PLOT, Device_bak
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
; open the post stript device
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        IF keyword_set(jpeg) THEN BEGIN
                info = size(Image)
                nx = info[1]
                ny = info[2]
                true_image = bytarr(3,nx,ny)
                true_image[0,*,*] = reds[image]
                true_image[1,*,*] = greens[image]
                true_image[2,*,*] = blues[image]
                write_jpeg, path+'test_Km_'+fecha+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'test_Km_'+fecha+'.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'test_Km_'+fecha
                print, ''
        ENDIF
        RETURN 
end

pro list, date_i, date_f
	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]	
;##############################################################################
; reading data files
;##############################################################################
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
                data_file_name[i] = '../rutidl/dH_teo/'+'teo_'+string_date[i]+'.dst.early'
        ENDFOR

        exist_data_file   = FILE_TEST(data_file_name)
        capable_to_plot   = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_plot NE N_ELEMENTS(data_file_name) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

        H    = FLTARR(file_number*24)
                      
        FOR i = 0, N_ELEMENTS(exist_data_file)-1 DO BEGIN
                IF exist_data_file[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        dat = DH_teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*24:(i+1)*24-1] = dat.H[*]
                ENDIF ELSE BEGIN
                         H[i*24:(i+1)*24-1] = 999999.0
                ENDELSE                
        ENDFOR
	idx = WHERE(H LE -100, count)
	IF count LE 0 THEN MESSAGE, 'ERROR'

    time_w = timegen(n_elements(file_number), final=julday(mh_f, dy_f, yr_f, 23), $
                start=julday(mh_i, dy_i, yr_i, 0) , units='H')

    caldat, time_w, m, d, y, hr
    y_idx    = y[idx]   
    m_idx    = m[idx]
    d_idx    = d[idx]
    H_idx    = H[idx]
    hr_idx   = hr[idx]
    
    y_idx_tmp    = y_idx*0
    m_idx_tmp    = m_idx*0
    d_idx_tmp    = d_idx*0
    H_idx_tmp    = H_idx*0
    hr_idx_tmp   = hr_idx*0
    
    j=0l
    y_idx_tmp[j]    = y_idx[0]
    m_idx_tmp[j]    = m_idx[0]
    d_idx_tmp[j]    = d_idx[0]
    H_idx_tmp[j]    = H_idx[0]
    hr_idx_tmp[j]   = hr_idx[0]
    
    FOR i=1l, N_ELEMENTS(H_idx)-1 DO BEGIN
	IF d_idx_tmp[j] EQ d_idx[i] THEN BEGIN
	        IF H_idx[i] GT H_idx_tmp[j] THEN BEGIN
        	        H_idx_tmp[j]   = H_idx[i]
		            y_idx_tmp[j]   = y_idx[i]
	                m_idx_tmp[j]   = m_idx[i]
        	        d_idx_tmp[j]   = d_idx[i]
	                hr_idx_tmp[j]  = hr_idx[i]	               
	        ENDIF
	ENDIF ELSE BEGIN
                j+=1
                
        	    H_idx_tmp[j]   = H_idx[i]
		        y_idx_tmp[j]   = y_idx[i]
	            m_idx_tmp[j]   = m_idx[i]
        	    d_idx_tmp[j]   = d_idx[i]
	            hr_idx_tmp[j]  = hr_idx[i]	
	ENDELSE
    ENDFOR
    
    outfile='../rutidl/output/'+'DH_list_TGM.txt'
        OPENW, lun, outfile, /GET_LUN, /append
        
    print, '###################################################################'
    print, 'lista de eventos de tormenta geomagnética '
    print, '###################################################################'
    print, '                                                                    '
    print, 'Descripción: Identificación de la fecha y hora en formato '
    print, 'fecha yy/mm/dd hh, día del año de cuando el índice DH descendió por'
    print, 'debajo de -100 nT, significando un evento de tormenta geomagnética '
    print, 'intensa.'
    print, '                                                                    '    
    print, format='(2X, "Fecha", 4X, "Hora", 2X, "Indice DH")'
    print, '---------------------------------'
    print, '                                                                    '
         
    for i=0, N_ELEMENTS(idx)-1 do begin
           if d_idx_tmp[i] ne 0 then begin

                print, y_idx_tmp[i], m_idx_tmp[i], d_idx_tmp[i], hr_idx_tmp[i], $
                H_idx_tmp[i], $
                FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 4X, I04)' 
  
                printf, lun, y_idx_tmp[i], m_idx_tmp[i], d_idx_tmp[i], $
                hr_idx_tmp[i], H_idx_tmp[i], $ 
                FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 4X, I04)'                         
        endif        
    endfor
    
    close,lun
    FREE_LUN, lun
    print, '                                                                   '
    print, '###################################################################'
    print, 'A continuación, ejecutar el procedimiento plotting eligiendo dos'
    print, 'límites de ventana del tiempo en formato DOY entorno a cada evento'
 
end

pro list_kmex, date_i, date_f

	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]	
;##############################################################################
; reading data files
;##############################################################################
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
                data_file_name[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.final'
		        
		        file = FILE_SEARCH(data_file_name[i], COUNT=opened_files)
	            IF opened_files NE N_ELEMENTS(file) THEN begin
	                data_file_name[i] = '../rutidl/Kmex/'+'teo_'+string_date[i]+'.index.early'    
	            ENDIF		        	        
        ENDFOR
	
        exist_data_file   = FILE_TEST(data_file_name)
        capable_to_plot   = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_plot NE N_ELEMENTS(data_file_name) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')"              
        ENDIF

        fecha = string(yr_i, mh_i, dy_i, yr_f, mh_f, dy_f, format = '(I4,I02,I02,"_",I4,I02,I02)')

        k_mex    = INTARR(file_number*8)
        a_mex    = INTARR(file_number*8)

        FOR i = 0, N_ELEMENTS(exist_data_file)-1 DO BEGIN
                IF exist_data_file[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        dat = kmex([tmp_year, tmp_month, tmp_day])
                        
                        k_mex[i*8:(i+1)*8-1] = dat.k_mex[*]/10
                        a_mex[i*8:(i+1)*8-1] = dat.a_mex[*]
                ENDIF 
        ENDFOR

	idx = WHERE(k_mex GE 7, count)
	IF count LE 0 THEN MESSAGE, 'ERROR'

    time_w = timegen(n_elements(file_number), final=julday(mh_f, dy_f, yr_f, 23), $
                start=julday(mh_i, dy_i, yr_i, 0) , units='H')

    caldat, time_w, m, d, y, hr

    y_idx        = y[idx]   
    m_idx        = m[idx]
    d_idx        = d[idx]
    k_mex_idx    = k_mex[idx]
    hr_idx       = hr[idx]
    
    y_idx_tmp        = y_idx*0
    m_idx_tmp        = m_idx*0
    d_idx_tmp        = d_idx*0
    k_mex_idx_tmp    = k_mex_idx*0
    hr_idx_tmp       = hr_idx*0
;print, H[idx]
    j=0l
    y_idx_tmp[j]        = y_idx[0]
    m_idx_tmp[j]        = m_idx[0]
    d_idx_tmp[j]        = d_idx[0]
    k_mex_idx_tmp[j]    = k_mex_idx[0]
    hr_idx_tmp[j]       = hr_idx[0]
    
    FOR i=1l, N_ELEMENTS(k_mex_idx)-1 DO BEGIN
	IF d_idx_tmp[j] EQ d_idx[i] THEN BEGIN
	        IF k_mex_idx[i] GT k_mex_idx_tmp[j] THEN BEGIN
        	        k_mex_idx_tmp[j]   = k_mex_idx[i]
		            y_idx_tmp[j]       = y_idx[i]
	                m_idx_tmp[j]       = m_idx[i]
        	        d_idx_tmp[j]       = d_idx[i]
	                hr_idx_tmp[j]      = hr_idx[i]	               
	        ENDIF
	ENDIF ELSE BEGIN
                j+=1
                
        	    k_mex_idx_tmp[j]   = k_mex_idx[i]
		        y_idx_tmp[j]       = y_idx[i]
	            m_idx_tmp[j]       = m_idx[i]
        	    d_idx_tmp[j]       = d_idx[i]
	            hr_idx_tmp[j]      = hr_idx[i]	
	ENDELSE
    ENDFOR
    
    ;outfile='../rutidl/output/'+'Kmex_list_TGM.txt'
     ;   OPENW, lun, outfile, /GET_LUN, /append
        str = TeXtoIDL('\rho^2 + 2\Gamma_{ij}')
    print, '###################################################################'
    print, 'lista de eventos de tormenta geomagnética '
    print, '###################################################################'
    print, '                                                                    '
    print, 'Descripción: Identificación de la fecha y hora en formato '
    print, 'fecha yy/mm/dd hh, día del año de cuando el índice DH ascendió por'
    print, 'encima de 6, significando un evento de tormenta geomagnética '
    print, 'intensa.'
    print, '                                                                    '    
    print, format='(2X, "Fecha", 4X, "Hora", 2X, "Indice Kmex")'
    print, '---------------------------------'
    ;print, srt
         
    for i=0, N_ELEMENTS(idx)-1 do begin
           if d_idx_tmp[i] ne 0 then begin

                print, y_idx_tmp[i], m_idx_tmp[i], d_idx_tmp[i], hr_idx_tmp[i], $
                k_mex_idx_tmp[i], $
                FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 4X, I03)' 
  
               ; printf, lun, y_idx_tmp[i], m_idx_tmp[i], d_idx_tmp[i], $
               ; hr_idx_tmp[i], k_mex_idx_tmp[i], $ 
               ; FORMAT = '(I4, "-", I02, "-", I02, 2X, I02, 4X, I03)'                         
        endif        
    endfor
    ;close,lun
    ;FREE_LUN, lun
   ;PRINT, DATETIME[i], Dst[i], format = "(A29)"
   print, '                                                                    '
    print, '###################################################################'
    print, 'A continuación, ejecutar el procedimiento plotting eligiendo dos'
    print, 'límites de ventana del tiempo en formato DOY entorno a cada evento'
 
end

function get_list_dh, dh_flist
	On_error, 2
	compile_opt idl2, HIDDEN
	
		dh_flist = '../rutidl/output/DH_list_TGM.txt'
	
		file = FILE_SEARCH(dh_flist, COUNT=opened_dh_files)
		n_lines = FILE_LINES(file)
		data = STRARR(n_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun

        DStruct_dst = {year : 0, month : 0, day : 0, hour : 0, DH: 0}   

		r_DH = REPLICATE(DStruct_dst, n_lines)	
        ;print, nlines-header_dst-1, nlines-header_dst+1
        
		READS, data[0:n_lines-1], r_DH, $
		FORMAT='(I4,X,I02,X,I02,2X,I02,4X,I4)'
			
		return, r_DH; funciones solo retornan un valor, buscar juntar ambos archivos u otra salida
end

function get_list_km, km_flist
	On_error, 2
	compile_opt idl2, HIDDEN
	
		km_flist = '../rutidl/output/Kmex_list_TGM.txt'
	
		file = FILE_SEARCH(km_flist, COUNT=opened_km_files)
		n_lines = FILE_LINES(file)
		data = STRARR(n_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun

        DStruct_km = {year : 0, month : 0, day : 0, hour : 0, km : 0}   

		r_km = REPLICATE(DStruct_km, n_lines)	
        ;print, nlines-header_dst-1, nlines-header_dst+1
        
		READS, data[0:n_lines-1], r_km, $
		FORMAT='(I4,X,I02,X,I02,2X,I02,4X,I3)'
			
		return, r_km; funciones solo retornan un valor, buscar juntar ambos archivos u otra salida
end

pro list_ev, r_DH, r_km, dh_flist, km_flist
	On_error, 2
	compile_opt idl2, HIDDEN

		dh_flist = '../rutidl/output/DH_list_TGM.txt'
		km_flist = '../rutidl/output/Kmex_list_TGM.txt'
			
	    dh_dat = get_list_dh(dh_flist)
        km_dat = get_list_km(km_flist)
    
    outfile = '../rutidl/output/'+'TGM_Intensas_list_reg.txt'
    fecha_dh = string(dh_dat.year, dh_dat.month, dh_dat.day)
    fecha_km = string(km_dat.year, km_dat.month, km_dat.day)    
    
    OPENW, lun, Outfile, /GET_LUN, /append

print, '#######################################################################'
print, 'Lista de fechas donde tanto el índice Kmex como el DH cumplieron con el'
print, 'criterio que indica la ocurrencia de una tormenta geomagnética intensa '
print, '#######################################################################'
print, '                                                                       '
print, format = '(8X, "año", 11X, "mes", 9X, "día", 8X, "DH", X, "Kmex")'
print, '                                                                       '
DH = dh_dat.DH
km = km_dat.km

yr_h = dh_dat.year
mh_h = dh_dat.month
dy_h = dh_dat.day

yr_k = km_dat.year
mh_k = km_dat.month
dy_k = km_dat.day
fecha_k = fltarr(n_elements(km))
fecha_h = fltarr(n_elements(DH))

;idx = where(DH le -120)
;print, n_elements(idx)


    for i=0, n_elements(km)-1 do begin
        ;fecha_k[i] = string(yr_k[i], mh_k[i], dy_k[i], format = '(I4, X, I02, X, I02, X, I4)')
        print, fecha_k[i];, format = '(A16)'
        ;print, yr_h[i], mh_h[i], dy_h[i], DH[i], format = '(I4, X, I02, X, I02, X, I4)'
        for j = 0, n_elements(DH)-1 do begin
            ;fecha_h[j] = string(yr_h[j]+mh_h[j]+dy_h[j], format = '(I4, X, I02, X, I02, X, I4)')
            ;print, yr_k[j], mh_k[j], dy_k[j], km[j], format = '(I4, X, I02, X, I02, X, I3)'
           
            if yr_k[i] eq yr_h[j] && mh_k[i] eq mh_h[j] && dy_k[i] eq dy_h[j] then begin
            ;if fecha_k[i] eq fecha_h[j] then begin
                ;yr_k[i]    = yr_h[j]
                ;mh_k[i]    = mh_h[j]
                ;dy_k[i]    = dy_h[j]
                 print, yr_k[i], mh_k[i], dy_k[i], format = '(I4, X, I02, X, I02)'             
               ; print, yr_k[i], mh_k[i], dy_k[i], km[i],$
                ;format='(I4, X, I02, X, I02, 2X, I3)'
                ;km_dat.km[j], format = '(I4, X, I02, X, I02, 2X, I4, 2X, I02)'
                
                
                ;printf, lun, km_dat.year[j], km_dat.month[j], km_dat.day[j], $
                ;dh_dat.DH[i], km_dat.km[j],$
                ;format = '(4I, X, I02, X, I02, 2X, I4, 2X, I02)'
            endif        
        endfor    
    endfor

    close,lun
    FREE_LUN, lun

end














		
		
		
		
		
		
		

