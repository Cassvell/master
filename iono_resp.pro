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
FUNCTION tec_data, idate
	On_error, 2
	compile_opt idl2, HIDDEN

	iyear	= idate[0]
	imonth	= idate[1]
	iday 	= idate[2]		

        header = 1      ; Defining number of lines of the header 
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
        idate = string(iyear, imonth, iday, format = '(I4, "-", I02, "-", I02)')
		file_name = path+'/tec/'+'tec_'+idate+'.txt'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun

        DStruct = {doy : 0, tec : 0., med : 0.}                                   
		r_tec = REPLICATE(DStruct, number_of_lines-header)	      
		READS, data[header:number_of_lines-1], r_tec, $
	format='(F3, X, F5, X, F5)'		
		RETURN, r_tec
END

FUNCTION dst_data, date

	On_error, 2
	COMPILE_OPT idl2, HIDDEN
    
	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
        date = string(year, month, day, format = '(I4,"-",I02,"-",I02)')
		file_name = path+'/dst/daily/dst_'+date+'.txt'
		;print, date
        header = 1             ; Defining number of lines of the header 
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-	
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'		

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun

        DataStruct = {year : 0, month : 0, day : 0, hour : 0, minute: 0, $
        second : 0, DOY : 0, Dst: 0}

		r_dst = REPLICATE(DataStruct, number_of_lines-header)	        
        
		READS, data[header:number_of_lines-1], r_dst, $
		FORMAT='(I4,X,I2,X,I2,X,I2,X,I2,X,I2,I04,I5)'
		RETURN, r_dst
END

FUNCTION teo, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
;reading data files
;-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        date = string(year, month, day, format = '(I4, I02, I02)')
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'		
		file_name = path+'/teoloyucan/hourly/teo_'+date+'h.dat'
		
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
        DStruct = {H : 0.}

		teo_mag = REPLICATE(DStruct, number_of_lines)	
  
		READS, data[0:number_of_lines-1], teo_mag, $
		FORMAT='(F10)'		
		RETURN, teo_mag		
END



FUNCTION baseline_sq, date
	On_error, 2
	compile_opt idl2, HIDDEN
	
	yr	= date[0]
	mh	= date[1]
	dy  = date[2]

        date = string(yr, mh, dy, format = '(I4,I02,I02)')
        header=0

		file_name = '../rutidl/output/Bsq_baselines/'+'Bsq_'+date+'h.dat'
	
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
		
        DStruct = {Bsq : 0.}
		B_sq = REPLICATE(DStruct, number_of_lines-header)	
		READS, data[header:number_of_lines-1], B_sq, $
	 format='(F08.4)' 		
		RETURN, B_sq
END

	
PRO iono_resp, date_i, date_f, JPEG = jpeg 

	On_error, 2
	compile_opt idl2, HIDDEN

	yr_i	= date_i[0]
	mh_i	= date_i[1]
	dy_i 	= date_i[2]	

	yr_f	= date_f[0]
	mh_f	= date_f[1]
	dy_f 	= date_f[2]
    file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1 	
;###############################################################################
idate0 = string(yr_i, mh_i, format='(I4,I02)')
TGM_n = idate0
case TGM_n of
    '200311' : TGM_n = 1
    '200411' : TGM_n = 2
    '200505' : TGM_n = 3
    '201503' : TGM_n = 4
    '201705' : TGM_n = 5
    '201709' : TGM_n = 6
    else: print, 'fuera de rango'
endcase  
;###############################################################################   
    dst_days= findgen(file_number*24)/24.0    
    Date    = string(yr_i, mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
;###############################################################################
; define DH variables
;###############################################################################
        string_date        = strarr(file_number)
               
        data_file_name_h  = strarr(file_number)                        
        data_file_name_bsq  = strarr(file_number)        
        data_file_name_dst = strarr(file_number)                
        data_file_name_tec  = strarr(file_number)        
        string_date_2        = strarr(file_number)
        
     data_path = '/home/c-isaac/Escritorio/proyecto/master_thesis/datos'                         
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                string_date_2[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4, "-", I02, "-", I02)')
        
                data_file_name_h[i]  = data_path+'/teoloyucan/hourly/'+'teo_'+string_date[i]+'h.dat'
                data_file_name_tec[i] = data_path+'/tec/'+'tec_'+string_date_2[i]+'.txt'
                data_file_name_bsq[i] = '../rutidl/output/Bsq_baselines/Bsq_'+string_date[i]+'h.dat'
                data_file_name_dst[i]= data_path+'/dst/daily/dst_'+string_date_2[i]+'.txt'                
		        
		        file_h  = FILE_SEARCH(data_file_name_h[i], COUNT=opened_files)
		        file_tec = FILE_SEARCH(data_file_name_tec[i], COUNT=opened_files)
		        file_bsq = FILE_SEARCH(data_file_name_bsq[i], COUNT=opened_files)
	            IF opened_files NE N_ELEMENTS(file_dh) THEN begin
	                data_file_name_h[i] = data_path+'/dH_teo/'+'teo_'+string_date[i]+'.dst.early'    
	            ENDIF	            	                            
        ENDFOR

        exist_data_file_h   = FILE_TEST(data_file_name_h)
        capable_to_plot_h   = N_ELEMENTS(where(exist_data_file_h EQ 1))

        exist_data_file_tec   = FILE_TEST(data_file_name_tec)
        capable_to_plot_tec   = N_ELEMENTS(where(exist_data_file_tec EQ 1))
        
        exist_data_file_bsq   = FILE_TEST(data_file_name_bsq)
        capable_to_plot_bsq   = N_ELEMENTS(where(exist_data_file_bsq EQ 1))        
        
        IF capable_to_plot_h NE N_ELEMENTS(data_file_name_h) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dh_index.',A,' impossible to plot all data.')"              
        ENDIF

        IF capable_to_plot_bsq NE N_ELEMENTS(data_file_name_bsq) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.Bsq.txt',A,' impossible to plot all data.')"              
        ENDIF

        IF capable_to_plot_tec NE N_ELEMENTS(data_file_name_tec) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.tec_index.',A,' impossible to plot all data.')"              
        ENDIF
;###############################################################################
; Generate the time variables to plot H time series                     
        H    = FLTARR(file_number*24)                       
        FOR i = 0, N_ELEMENTS(exist_data_file_h)-1 DO BEGIN
                IF exist_data_file_h[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        d_h = teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*24:(i+1)*24-1] = d_h.H[*]
                                                                                                                       
                ENDIF ELSE BEGIN
                        H[i*24:(i+1)*24-1] = 999999.0
                ENDELSE                
        ENDFOR
;###############################################################################
; Generate the time variables to plot TEC time series         
        tec  = fltarr(file_number*12)
        med  = fltarr(file_number*12)
        FOR i = 0, N_ELEMENTS(exist_data_file_tec)-1 DO BEGIN
                IF exist_data_file_tec[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date_2[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,X,I02,X,I02)'                 
                        d_tec = tec_data([tmp_year, tmp_month, tmp_day])
                        
                        tec[i*12:(i+1)*12-1] = d_tec.tec[*]
                        med[i*12:(i+1)*12-1] = d_tec.med[*] 
                ENDIF ELSE BEGIN
                        tec[i*12:(i+1)*12-1] = 999.0
                        med[i*12:(i+1)*12-1] = 999.0
                ENDELSE                
        ENDFOR
;###############################################################################
;Dst Data                       
        exist_data_file_dst   = FILE_TEST(data_file_name_dst)
        capable_to_plot_dst   = N_ELEMENTS(where(exist_data_file_dst EQ 1))

        IF capable_to_plot_dst NE N_ELEMENTS(data_file_name_dst) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dst_index.',A,' impossible to plot all data.')"              
        ENDIF

        dst    = FLTARR(file_number*24)                               
        FOR i = 0, N_ELEMENTS(exist_data_file_dst)-1 DO BEGIN
                IF exist_data_file_dst[i] EQ 1 THEN BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date_2[i] = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,"-",I02,"-",I02)')                
                        dat = dst_data([tmp_year, tmp_month, tmp_day])
                        
                        dst[i*24:(i+1)*24-1] = dat.Dst[*]                                                
                                                                                              
                ENDIF ELSE BEGIN
                         dst[i*24:(i+1)*24-1] = 999999.0                      
                ENDELSE                
        ENDFOR   
;###############################################################################
;Identifying the NAN values        
        for i=0, n_elements(tec)-1 do begin
            if tec[i] eq 999.0 then begin
                tec[where(tec[*] eq 999.0)] = !Values.F_NAN          
            endif
        endfor    

        for i=0, n_elements(med)-1 do begin
            if med[i] eq 999.0 then begin
                med[where(med[*] eq 999.0)] = !Values.F_NAN          
            endif
        endfor
        
        i_nan1 = where(H eq 999999.0, ncount)
        i_nan2 = where(H gt 100.0, n2count)
        
        prcent_nan = FLOAT(ncount+n2count)*100.0
        print,'porcentaje de valores NaN:', prcent_nan/n_elements(H),'%'
        
        for i=0, n_elements(H)-1 do begin
            if H[i] eq 999999.0 then begin
                H[where(H[*] eq 999999.0)] = !Values.F_NAN          
            endif
        endfor
        
        for i=0, n_elements(H)-1 do begin
            if H[i] ge 99999.0 then begin
                H[where(H[*] ge 99999.0)] = !Values.F_NAN          
            endif
        endfor        
    ;implementar una función de interpolación en caso de que el porcentaje de 
    ;nan sea muy bajo
       
    H_tmp   = H
    H_exist = where(finite(H_tmp), ngooddata, complement=baddata, ncomplement=nbaddata)
    ; interpolate at the locations of the bad data using the good data
    
    if nbaddata gt 0 then H_tmp[baddata] = interpol(H_tmp[H_exist], H_exist, baddata)
    H = H_tmp
;###############################################################################      
    tec_days= findgen(file_number*12)/12.0                         
    tec_diff = tec-med    
    new_tecdays = findgen(file_number*1440)/1440.0 ;se genera un arreglo de tiempo con 
;    muestreo cada 15 min. para mejorar la resolución de las gráficas    
    
    new_tecdiff = FLTARR(N_ELEMENTS(new_tecdays))     	    
    tmp_tecdif  = INTERPOL(tec_diff, N_ELEMENTS(new_tecdays))
    new_tecdiff = tmp_tecdif
;###############################################################################
    new_dstdays = findgen(file_number*1440)/1440.0 ;se genera un arreglo de tiempo con 
;    muestreo cada 15 min. para mejorar la resolución de las gráficas    
    
    new_dst = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_dst  = INTERPOL(dst, N_ELEMENTS(new_dstdays))
    new_dst = tmp_dst 
    
    new_H = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_H  = INTERPOL(H, N_ELEMENTS(new_dstdays))
    new_H = tmp_H     
;############################################################################### 
; define diurnal baseline
;###############################################################################  
        Bsq    = FLTARR(file_number*24)                       
        FOR i = 0, N_ELEMENTS(exist_data_file_bsq)-1 DO BEGIN
                IF exist_data_file_bsq[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        sqline = baseline_sq([tmp_year, tmp_month, tmp_day])
                        
                        Bsq[i*24:(i+1)*24-1] = sqline.Bsq[*]
                                                                                                                       
                ENDIF                
        ENDFOR
;###############################################################################
; define frequencies
;###############################################################################  
   mlat         = 28.06*!pi
   ld           = cos(mlat/180)
   p_a          = dst*ld
   baseline     = Bsq + p_a             
        diono   = H-baseline
    n           = n_elements(diono) 
time = 3600.0

fn      = FLOAT(1.0/(2.0*time)) ; frecuencia de Nyquist
y       = FFT(diono)

pws     = abs(y[0:n/2])^2
pws_s   = smooth(pws, 1)

f_k     = (1+findgen(n))/(n*time)
print, 'Nyquist freq: ', fn, 'Hz'
;###############################################################################
    i_diff = diono
    new_idiff = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_idiff  = INTERPOL(i_diff, N_ELEMENTS(new_dstdays))
    new_idiff = tmp_idiff  
;###############################################################################
; define pass band frequencies
;###############################################################################  
passband_l = idate0
CASE passband_l of
    '200311' : passband_l = 1.2e-5  ;23:08 hr
    '200411' : passband_l = 1e-5    ;27:46
    '200505' : passband_l = 1.2e-5  ;23:08
    '201503' : passband_l = 1.14e-5 ;24:21
    '201705' : passband_l = 0.9e-5  ;30:51
    '201709' : passband_l = 0.95e-5  ;23:08
    ELSE: PRINT, 'fuera de rango'
ENDCASE  
;###############################################################################
passband_u = idate0
CASE passband_u of
    '200311' : passband_u = 3.2e-5  ;08:40
    '200411' : passband_u = 2.5e-5  ;11:06
    '200505' : passband_u = 2.8e-5  ;09:55
    '201503' : passband_u = 2.8e-5  ;09:55
    '201705' : passband_u = 3.2e-5  ;08:40
    '201709' : passband_u = 3.4e-5  ;08:10
    ELSE: PRINT, 'fuera de rango'
ENDCASE  
;###############################################################################
; define high band frequencies
highpass_l = idate0
CASE highpass_l of
    '200311' : highpass_l = 7e-5    ;03:58
    '200411' : highpass_l = 7.5e-5  ;03:42
    '200505' : highpass_l = 7.3e-5  ;03:48
    '201503' : highpass_l = 7.5e-5  ;03:42
    '201705' : highpass_l = 7.5e-5  ;03:42
    '201709' : highpass_l = 7e-5    ;03:58
    ELSE: PRINT, 'fuera de rango'
ENDCASE
;###############################################################################
; define filtering
;############################################################################### 
coeff_ddyn    = digital_filter(passband_l/fn, passband_u/fn, 50, 18)
coeff_dp2   = digital_filter(highpass_l/fn, 1.0, 50, 4)
;###############################################################################
; define disturbing effects
;############################################################################### 
ddyn        = convol(diono, coeff_ddyn, /edge_wrap)
dp2         = convol(diono, coeff_dp2, /edge_wrap)  
;###############################################################################      
    new_ddyn = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_ddyn  = INTERPOL(ddyn, N_ELEMENTS(new_dstdays))
    new_ddyn = tmp_ddyn           
;############################################################################### 
    new_dp2 = FLTARR(N_ELEMENTS(new_dstdays))     	    
    tmp_dp2  = INTERPOL(dp2, N_ELEMENTS(new_dstdays))
    new_dp2 = tmp_dp2         
;###############################################################################
; define device and color parameters 
;###############################################################################      
        Device_bak2 = !D.Name         
        SET_PLOT, 'Z'      
        
        Xsize=fix(1600)
        Ysize=1000
        DEVICE, SET_RESOLUTION = [Xsize,Ysize],Set_Pixel_Depth=24, DECOMPOSED=1  
        DEVICE, z_buffer=4
        DEVICE, set_character_size = [10, 12] 
        
        chr_size1 = 0.9
        chr_thick1= 1.5
        space     = 0.015
        rojo      = 248
        amarillo  = 190
        verde     = 150
        negro     = 0
        azul      = 70
        blanco    = 255
        gris      = 110
        morado    = 16
        naranja  = 220
                
    TVLCT, R_bak, G_bak, B_bak, /GET
        
    LOADCT, 39, /SILENT
     X_label   = STRARR(file_number+1)+' '
        months    = ['Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic']
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
        
case old_month of
    1: old_month = 'Enero'
    2: old_month ='Febrero'
    3: old_month ='Marzo'
    4: old_month ='Abril'
    5: old_month ='Mayo'
    6: old_month ='Junio'
    7: old_month ='Julio'
    8: old_month ='Agosto'
    9: old_month ='Septiembre'
    10:old_month ='Octubre'
    11:old_month ='Noviembre'
    12:old_month ='Diciembre'
    else: print, 'fuera de rango'
endcase             
;###############################################################################
    spam_i = idate0
case spam_i of
    '200311' : spam_i = 1000
    '200411' : spam_i = 50
    '200505' : spam_i = 0
    '201503' : spam_i = 0
    '201705' : spam_i = 0
    '201709' : spam_i = 0
    else: print, 'fuera de rango'
endcase 

    spam_f = idate0
case spam_f of
    '200311' : spam_f = 2800
    '200411' : spam_f = 4400
    '200505' : spam_f = 0
    '201503' : spam_f = 0
    '201705' : spam_f = 0
    '201709' : spam_f = 0
    else: print, 'fuera de rango'
endcase    
;############################################################################### 
       days = intarr(file_number+1)
       for i=0, n_elements(days)-1 do begin
            days[i] = dy_i+i
       endfor
       days = days*24/24. 
       day_time = findgen(24)   
;############################################################################### 
    time_title = ' Tiempo Universal ['+textoidl("dias")+' de '+old_month+'].'
    window_title = 'TGM'+ STRING(TGM_n, FORMAT='(I01)')+', '+ $
                    STRING(old_month, yr_i, FORMAT='(A, X, I4)')
    
    periodo = 'Periodo [h]'        
;###############################################################################               
    PLOT, f_k, pws_s, /XLOG, /YLOG, POSITION=[0.07,0.1,0.95,0.9],$
    BACKGROUND = blanco, COLOR=negro, $
    CHARSIZE = chr_size1, XSTYLE=5, YSTYLE=5, SUBTITLE='', THICK=4, /NODATA    

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.95   
   XYOUTS, X, y, window_title, /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=2, CHARTHICK=1.5               
;###############################################################################           
    ysup = MAX(pws)+10
    yinf = MIN(pws)-0.0001
    
    freqs = [1.0/(96.0*3600.0), 1.0/(48.0*3600.0), 1.0/(24.0*3600.0), $
              1.0/(12.0*3600.0), 1.0/(6.0*3600.0), 1.0/(3.0*3600.0)]
               
    periods = [96.0, 48.0, 24.0, 12.0, 6.0, 3.0]
    
    PLOT, f_k, pws_s, /XLOG, /YLOG, XRANGE = [freqs[0], fn], POSITION=[0.07,0.1,0.45,0.9],$
    YRANGE=[yinf, ysup], BACKGROUND = blanco, COLOR=negro, $
    CHARSIZE = chr_size1, XSTYLE=5, YSTYLE=5, SUBTITLE='', THICK=4, /NODATA,$
    /NOERASE


    POLYFILL, [passband_l, passband_u ,passband_u, passband_l], $
              [!Y.CRANGE[0], !Y.CRANGE[0], ysup, ysup], COLOR=amarillo

    POLYFILL, [highpass_l, fn ,fn, highpass_l], $
              [!Y.CRANGE[0], !Y.CRANGE[0], ysup, ysup], COLOR=amarillo
    OPLOT, f_k, pws_s, COLOR=negro, THICK=5    
;###############################################################################    
        AXIS, XAXIS = 0, XRANGE=[freqs[0], fn], $
                         /XLOG,$
                         XSTYLE=1,$
                         xTITLE = 'Frecuencia [Hz]',$
                         COLOR=negro, $
                         CHARSIZE = 1.2, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5
                                           
        AXIS, XAXIS = 1, XRANGE=[freqs[0], fn], $;.0/(!X.CRANGE), $
                         /XLOG,$
                         XTICKS=6,$
                         XMINOR=4,$
                         XTICKV=freqs,$                         
                         XTICKN=STRING(periods, FORMAT='(F4.1)'),$
                         XSTYLE=1,$
                         CHARSIZE = 1.2,$
                         COLOR=negro, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5                     

        AXIS, YAXIS = 0, yrange=[yinf, ysup], $
                         YTITLE = '', $
                         ystyle=1,$                          
                         COLOR=negro, $
                         /ylog,$
                         CHARSIZE = 1.0,$
                         CHARTHICK=1.5
                        
        AXIS, YAXIS = 1, yrange=[yinf, ysup], $
                         COLOR=negro, $
                         /ylog,$
                         ystyle=1, $
                         CHARSIZE = 1.0,$
                         CHARTHICK=1.5

   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOUTS, X, 0.928, periodo, /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.2, CHARTHICK=1.5   
   
   y = (!Y.Window[1] - !Y.Window[0]) / 2. + !Y.Window[0] 
   XYOUTS, 0.02, y, 'Componente espectral [nT]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.2, ORIENTATION=90, CHARTHICK=1.5     
;###############################################################################
    med_ddyn = MEDIAN(new_ddyn)
    std_ddyn = stddev(new_ddyn, /NAN)
    
    ddyn_out = WHERE(new_ddyn GE med_ddyn+std_ddyn OR new_ddyn LE med_ddyn-std_ddyn)
    ddyn_in  = WHERE(new_ddyn LE med_ddyn+std_ddyn AND new_ddyn GE med_ddyn-std_ddyn)
    
    ddyn_diff_out = new_ddyn
    ddyn_diff_out[ddyn_in]=!Values.F_NAN
    
    ddyn_diff_in  = new_ddyn
    ddyn_diff_in[ddyn_out]=!Values.F_NAN
 
     upddyn     = max(ddyn)
     downddyn   = min(ddyn)
     
     updp2     = max(dp2)
     downdp2   = min(dp2)     
;###############################################################################
    ddyn_i = idate0
CASE ddyn_i of
    '200311' : ddyn_i = ddyn_out[0]
    '200411' : ddyn_i = ddyn_out[0]
    '200505' : ddyn_i = ddyn_out[500]
    '201503' : ddyn_i = ddyn_out[400]
    '201705' : ddyn_i = ddyn_out[0]
    '201709' : ddyn_i = ddyn_out[400]
    ELSE: PRINT, 'fuera de rango'
ENDCASE 

    ddyn_si = idate0
CASE ddyn_si of
    '200311' : ddyn_si = -200
    '200411' : ddyn_si = -130
    '200505' : ddyn_si = -160
    '201503' : ddyn_si = -650
    '201705' : ddyn_si = -130
    '201709' : ddyn_si = -350
    ELSE: PRINT, 'fuera de rango'
ENDCASE 

    ddyn_sf = idate0
CASE ddyn_sf of
    '200311' : ddyn_sf = -800
    '200411' : ddyn_sf = 600
    '200505' : ddyn_sf = 5370
    '201503' : ddyn_sf = 4780
    '201705' : ddyn_sf = 1800
    '201709' : ddyn_sf = 4170
    ELSE: PRINT, 'fuera de rango'
ENDCASE       
;###############################################################################   
     dH = TeXtoIDL('\DeltaH') 
;###############################################################################      
    med_tec = MEDIAN(new_tecdiff)
    std_tec = stddev(new_tecdiff)
    
    index_out = WHERE(new_tecdiff GE med_tec+std_tec OR new_tecdiff LE med_tec-std_tec)
    index_in  = WHERE(new_tecdiff LE med_tec+std_tec AND new_tecdiff GE med_tec-std_tec)
    tec_diff_out = new_tecdiff
    tec_diff_out[index_in]=!Values.F_NAN
    
    tec_diff_in  = new_tecdiff
    tec_diff_in[index_out]=!Values.F_NAN

    sup = med_tec+std_tec
    inf = med_tec-std_tec
    
    l_sup = fltarr(n_elements(tec_diff))
    l_sup[*] = sup

    l_inf = fltarr(n_elements(tec_diff))
    l_inf[*] = inf  
     
    up_tecdiff      = MAX(tec_diff)
    down_tecdiff    = MIN(tec_diff) 
    PLOT, tec_days, tec_diff, XTICKS=file_number, xminor=8, BACKGROUND = blanco,$ 
     CHARSIZE = chr_size1, CHARTHICK=chr_thick1, POSITION=[0.55,0.1,0.95,0.27], $
     XSTYLE = 5, XRANGE=[0, file_number], XTICKNAME=REPLICATE(' ', file_number+1), ySTYLE = 6,$
     /NOERASE, YRANGE=[down_tecdiff, up_tecdiff], /NODATA
    
    POLYFILL, [new_dstdays[ddyn_i+ddyn_si], new_dstdays[ddyn_i+spam_f+ddyn_sf],$
              new_dstdays[ddyn_i+spam_f+ddyn_sf], new_dstdays[ddyn_i+ddyn_si]], $
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=amarillo  

     OPLOT, tec_days, tec_diff, LINESTYLE=0, THICK=4, COLOR=morado  
             
     OPLOT, tec_days, l_sup, LINESTYLE=2, THICK=2, COLOR=morado
     OPLOT, tec_days, l_inf, LINESTYLE=2, THICK=2, COLOR=morado
                        
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XTITLE=time_title, $                         
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $
                         XTICKS=file_number, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $                                            
                         XMINOR=8, $ 
                         CHARSIZE = 0.9, $                       
                         COLOR=negro, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5

        AXIS, YAXIS = 0, YRANGE=[down_tecdiff, up_tecdiff], $
                         YTITLE = '', $                          
                         COLOR=negro, $
                         YSTYLE=2, $
                         CHARSIZE = 1.0,$
                         CHARTHICK=1.5
                        
        AXIS, YAXIS = 1, YRANGE=[down_tecdiff, up_tecdiff], $
                         COLOR=negro, $
                         YSTYLE=2, $
                         CHARSIZE = 1.0,$
                         CHARTHICK=1.5                                       
;############################################################################### 
     up_diono=max(diono)
     down_diono=min(diono)          
     PLOT, dst_days, diono, XTICKS=file_number, XMINOR=8, BACKGROUND = blanco, $
     COLOR=negro, CHARSIZE = 0.6, CHARTHICK=chr_thick1, $
     POSITION=[0.55,0.73,0.95,0.9], XSTYLE = 5, XRANGE=[0, file_number], ySTYLE = 6,$
     XTICKNAME=REPLICATE(' ', file_number+1), YRANGE=[down_diono,up_diono], /NOERASE,$
     THICK=4, /NODATA   

    POLYFILL, [new_dstdays[ddyn_i+ddyn_si], new_dstdays[ddyn_i+spam_f+ddyn_sf],$
              new_dstdays[ddyn_i+spam_f+ddyn_sf], new_dstdays[ddyn_i+ddyn_si]], $
              [!Y.CRANGE[0], !Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1]], $
              COLOR=amarillo   
     
     OPLOT, dst_days, diono, THICK=4, LINESTYLE=0, COLOR=negro   
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XTITLE=time_title, $                         
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5
                         
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-0.25), $
                         XTICKS=file_number, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $                                            
                         XMINOR=8, $ 
                         CHARSIZE = 0.9, $                       
                         COLOR=negro, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5
    ppi = TexToIDL('P_{PI}')
        AXIS, YAXIS = 0, YRANGE=[down_diono,up_diono], $
                         YTITLE = '', $                          
                         COLOR=negro, $
                         YSTYLE=2, $
                         CHARSIZE = 1.1,$
                         CHARTHICK=1.5
                        
        AXIS, YAXIS = 1, YRANGE=[down_diono,up_diono], $      
                         COLOR=negro, $
                         YSTYLE=2, $
                         CHARSIZE = 1.1,$
                         CHARTHICK=1.5                   
;###############################################################################                
    if max(ddyn) gt max(dp2) then up = max(ddyn) else up = max(dp2)
    if min(ddyn) lt min(dp2) then down = min(ddyn) else down = min(dp2)
;###############################################################################
    IF upddyn GT updp2 THEN up = upddyn ELSE up=updp2 
    IF downddyn LT downdp2 THEN down = downddyn ELSE down=downdp2 
                               
     PLOT, dst_days, ddyn, XTICKS=file_number, XMINOR=8, BACKGROUND = blanco, $
     COLOR=negro, CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
     POSITION=[0.55,0.34,0.95,0.66], XSTYLE = 5, XRANGE=[0, file_number], YSTYLE = 6,$
     XTICKNAME=REPLICATE(' ', file_number+1), YRANGE=[down,up], /NOERASE, /NODATA
    
        OPLOT, new_dstdays[ddyn_i+ddyn_si:ddyn_i+spam_f+ddyn_sf], $
        ddyn_diff_out[ddyn_i+ddyn_si:ddyn_i+spam_f+ddyn_sf], $
        COLOR=negro, LINESTYLE=0, THICK=5   
                    
        OPLOT, new_dstdays[ddyn_i+ddyn_si:ddyn_i+spam_f+ddyn_sf], $
        ddyn_diff_in[ddyn_i+ddyn_si:ddyn_i+spam_f+ddyn_sf], $
        COLOR=negro, LINESTYLE=0, THICK=5            
;###############################################################################
        OPLOT, new_dstdays, new_ddyn, COLOR=negro, LINESTYLE=3
       ; OPLOT, new_dstdays, ddyn_diff_out, COLOR=negro, LINESTYLE=0, THICK=5
;###############################################################################
    med_dp2 = MEDIAN(new_dp2)
    std_dp2 = stddev(new_dp2, /NAN)
    
    dp2_out = WHERE(new_dp2 GE med_dp2+std_dp2 OR new_dp2 LE med_dp2-std_dp2)
    dp2_in  = WHERE(new_dp2 LE med_dp2+std_dp2 AND new_dp2 GE med_dp2-std_dp2)
    
    dp2_diff_out = new_dp2
    dp2_diff_out[dp2_in]=!Values.F_NAN
    
    dp2_diff_in  = new_dp2
    dp2_diff_in[dp2_out]=!Values.F_NAN     
;###############################################################################
    dp2_i = idate0
case dp2_i of
    '200311' : dp2_i = dp2_out[6]
    '200411' : dp2_i = dp2_out[50]
    '200505' : dp2_i = dp2_out[100]
    '201503' : dp2_i = dp2_out[500]
    '201705' : dp2_i = dp2_out[100]
    '201709' : dp2_i = dp2_out[500]
    else: print, 'fuera de rango'
endcase  

    dp2_si = idate0
case dp2_si of
    '200311' : dp2_si = 0
    '200411' : dp2_si = -100
    '200505' : dp2_si = -100
    '201503' : dp2_si = -730
    '201705' : dp2_si = -50
    '201709' : dp2_si = 0
    else: print, 'fuera de rango'
endcase 

    dp2_sf = idate0
case dp2_sf of
    '200311' : dp2_sf = 20
    '200411' : dp2_sf = 100
    '200505' : dp2_sf = 350
    '201503' : dp2_sf = 550
    '201705' : dp2_sf = 900
    '201709' : dp2_sf = 990
    else: print, 'fuera de rango'
endcase       
;###############################################################################     
     oplot, dst_days, dp2, color=rojo
        oplot, new_dstdays[dp2_i+spam_i+dp2_si:dp2_i+spam_f+dp2_sf], $
        new_dp2[dp2_i+spam_i+dp2_si:dp2_i+spam_f+dp2_sf], color=rojo, $
        linestyle=0, thick=4      
        
        oplot, new_dstdays[dp2_i+spam_i+dp2_si:dp2_i+spam_f+dp2_sf], $
        new_dp2[dp2_i+spam_i+dp2_si:dp2_i+spam_f+dp2_sf], color=rojo, $
        linestyle=0, thick=4          
;############################################################################### 
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XTITLE=time_title, $                         
                         XMINOR=8, $
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5
                         
                         IF TGM_n EQ 3 THEN Ltime = 0.25
                         IF TGM_n EQ 5 THEN Ltime = 0.25
                         IF TGM_n EQ 6 THEN Ltime = 0.25
                         IF TGM_n EQ 1 THEN Ltime = 5./24.
                         IF TGM_n EQ 2 THEN Ltime = 5./24.
                         IF TGM_n EQ 4 THEN Ltime = 5./24.                                                  
        AXIS, XAXIS = 1, XRANGE=(!X.CRANGE+dy_i-Ltime), $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKV=FIX(days), $       
                         XTICKN=STRING(days, FORMAT='(I02)'), $  
                         CHARSIZE = 0.8, $                                                
                         COLOR=negro, $
                         TICKLEN=0.04,$
                         CHARTHICK=1.5

        AXIS, YAXIS = 0, yrange=[down,up], $ 
                         ystyle=2, $  
                         YTITLE = '', $                          
                         COLOR=negro, $
                         CHARSIZE = 1.1,$
                         CHARTHICK=1.5
                        
        AXIS, YAXIS = 1, yrange=[down,up], $ 
                         ystyle=2, $ 
                         COLOR=negro, $
                         CHARSIZE = 1.1,$
                         CHARTHICK=1.5      
;###############################################################################    
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.921, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.9, CHARTHICK=1.5     

   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.679, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.9, CHARTHICK=1.5     
   
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   XYOuts, X, 0.291, 'Tiempo Local', /Normal, $
   color=negro, Alignment=0.5, Charsize=0.9, CHARTHICK=1.5   
   

   XYOuts, 0.93, 0.75, '(a)', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.4, CHARTHICK= 3   
   
   XYOuts, 0.11, 0.14, '(b)', /Normal, $
   color=negro, Alignment=0.5, Charsize=2.4, CHARTHICK= 3  
   
   XYOuts, 0.93, 0.63, '(c)', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.4, CHARTHICK= 3  
   
   XYOuts, 0.93, 0.24, '(d)', /Normal, $
   color=negro, Alignment=0.5, Charsize=1.4, CHARTHICK= 3     
;###############################################################################   
   y = (0.66 - 0.34) / 2. + 0.34 
   XYOUTS, 0.51, y, 'DP2 y Ddyn [nT]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.2, ORIENTATION=90, CHARTHICK=1.5   

   y = (0.9 - 0.73) / 2. + 0.73 
   XYOUTS, 0.51, y, ppi+' [nT]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.2, ORIENTATION=90, CHARTHICK=1.5   

   y = (0.27 - 0.1) / 2. + 0.1 
   XYOUTS, 0.51, y, 'TEC-<TEC> [TECu]', /NORMAL, $
   COLOR=negro, ALIGNMENT=0.5, CHARSIZE=1.2, ORIENTATION=90, CHARTHICK=1.5                
;###############################################################################                      
;second panel legend
        POLYFILL, [0.88,0.91,0.91,0.88], [0.405,0.405,0.407,0.407], color = rojo, /NORMAL
        POLYFILL, [0.88,0.91,0.91,0.88], [0.375,0.375,0.377,0.377], color = negro, /NORMAL        

        XYOUTS, 0.91, 0.4 , /NORMAL, $
                ' DP2', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1 

        XYOUTS, 0.91, 0.37 , /NORMAL, $
                ' Ddyn', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1                                    
;###############################################################################
; saving png
;###############################################################################     
     Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!    
    TVLCT, R_bak, G_bak, B_bak, /get  
        
    ;DEVICE, /CLOSE
    SET_PLOT, Device_bak2  
    path = '../rutidl/output/eventos_tgm/'
        IF keyword_set(jpeg) THEN BEGIN
                info = size(Image)
                nx = info[1]
                ny = info[2]
                true_image = bytarr(3,nx,ny)
                true_image[0,*,*] = R_bak[image]
                true_image[1,*,*] = G_bak[image]
                true_image[2,*,*] = B_bak[image]
                write_jpeg, path+'iono_resp_V8_'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'iono_resp_V8_'+Date+'.png', Image, R_bak, G_bak, B_bak
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'iono_resp_V7_'+Date+'.png'
                print, ''
        ENDIF
        RETURN 	
end	
