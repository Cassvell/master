;
;Name:
;	magdat_iono
;purpose:
;	Este programa está diseñado para retornar una gráfica de comparación entre 
;   la diferencia de los índices geomagnéticos Dst-DH y la diferencia del índice
;   ionosférico TEC - <TEC>. La gráfica resultante servirá para comparar la
;   respuesta ionosférica y la respuesta geomagnética con respecto del tiempo.
;author:
;	Carlos Isaac Castellanos Velazco
;	Estudiante de Maestría en Ciencias de la Tierra
;	Instituto de Geofísica, Unidad Michoacan
;	UNAM
;	ccastellanos@igeofisica.unam.mx
;
;category:
;   Obtención, graficado y comparación de datos magnéticos y ionosféricos. 
;
;calling sequence:
;>   .r mag_iono
;>   mag_iono, r_dst, doy, [idate], [fdate]
;parameters:
;
;
;dependencies:
;
;
;input files
;   archivos de DH
;   archivos de Dst
;   archivos de TEC
;output files:
;   Gráfica Dst-DH y TEC-<TEC>

FUNCTION dst_data, date

	On_error, 2
	COMPILE_OPT idl2, HIDDEN
    
	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;###############################################################################
;reading data files
        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'
        date = string(year, month, day, format = '(I4,"-",I02,"-",I02)')
		file_name = path+'/dst/daily/dst_'+date+'.txt'
        header = 1             ; Defining number of lines of the header 
;###############################################################################
;reading data files
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


FUNCTION DH_teo, date
	On_error, 2
	compile_opt idl2, HIDDEN

	year	= date[0]
	month	= date[1]
	day 	= date[2]	
;###############################################################################
;reading data files
        date = string(year, month, day, format = '(I4, I02, I02)')

        path='/home/c-isaac/Escritorio/proyecto/master_thesis/datos'		
		file_name = path+'/dH_teo/'+'teo_'+date+'.dst.early'
		
		file = FILE_SEARCH(file_name, COUNT=opened_files)
		IF opened_files NE N_ELEMENTS(file) THEN MESSAGE, file_name+' not found'

		number_of_lines = FILE_LINES(file)
	   ; print, number_of_lines
		data = STRARR(number_of_lines)

		openr, lun, file, /GET_LUN, ERROR=err
		readf, lun, data, FORMAT = '(A)'
		CLOSE, lun
		FREE_LUN, lun
;###############################################################################
;extracting data and denfining an structure data
        DStruct = {hora : 0, D_stdesv : 0., D : 0., H_stdesv : 0., H : 0., $
        Z_stdesv : 0., Z : 0., N_stdesv : 0., N : 0., F_stdesv : 0., F : 0.}

		teo_mag = REPLICATE(DStruct, number_of_lines)	
  
		READS, data[0:number_of_lines-1], teo_mag, $
		FORMAT='(I2, F10, F8, F10, F10, F10, F10, F10, F10, F10, F10)'		
		RETURN, teo_mag		
END


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
	   ; print, number_of_lines
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


pro mag_iono, idate, fdate

	On_error, 2
	compile_opt idl2, HIDDEN
	
	yr_i	= idate[0]
	mh_i	= idate[1]
	dy_i 	= idate[2]	

	yr_f	= fdate[0]
	mh_f	= fdate[1]
	dy_f 	= fdate[2]	
;###############################################################################
; Generate the time variables to plot time series of Dst Index
    file_number    = (JULDAY(mh_f, dy_f, yr_f) - JULDAY(mh_i, dy_i, yr_i))+1 
    days_dst= findgen(file_number*24)/24.0    
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
;##############################################################################     
; Generate the time variables to plot time series of DH and TEC Index
        data_path = '/home/c-isaac/Escritorio/proyecto/master_thesis/datos'                           
               
        data_file_name_dh  = strarr(file_number)        
        string_date_dh        = strarr(file_number)

        data_file_name_tec  = strarr(file_number)
        data_file_name_dst = strarr(file_number)                
                        
        string_date_2        = strarr(file_number)        
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(mh_i, dy_i, yr_i)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                string_date_dh[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                string_date_2[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4, "-", I02, "-", I02)')
                data_file_name_dst[i]= data_path+'/dst/daily/dst_'+string_date_2[i]+'.txt'        
                data_file_name_dh[i] = data_path+'/dH_teo/'+'teo_'+string_date_dh[i]+'.dst.early'
                data_file_name_tec[i] = data_path+'/tec/tec_newformat/tec_'+string_date_2[i]+'.txt'
                                
		        file_dh = FILE_SEARCH(data_file_name_dh[i], COUNT=opened_files)
		        file_tec = FILE_SEARCH(data_file_name_tec[i], COUNT=opened_files)
		        
	            IF opened_files NE N_ELEMENTS(file_dh) THEN begin
	                data_file_name_dh[i] = '../rutidl/dH_teo/'+'teo_'+string_date_dh[i]+'.dst.early'    
	            ENDIF
	            
	            IF opened_files NE N_ELEMENTS(file_tec) THEN begin
	                data_file_name_tec[i] = '../rutidl/tec/tec_newformat/'+'tec_'+string_date_tec[i]+'.txt'
	            ENDIF 	                            
        ENDFOR

        exist_data_file_dh   = FILE_TEST(data_file_name_dh)
        capable_to_plot_dh   = N_ELEMENTS(where(exist_data_file_dh EQ 1))

        exist_data_file_tec   = FILE_TEST(data_file_name_tec)
        capable_to_plot_tec   = N_ELEMENTS(where(exist_data_file_tec EQ 1))

        exist_data_file_dst   = FILE_TEST(data_file_name_dst)
        capable_to_plot_dst   = N_ELEMENTS(where(exist_data_file_dst EQ 1))

        IF capable_to_plot_dst NE N_ELEMENTS(data_file_name_dst) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dst_index.',A,' impossible to plot all data.')"              
        ENDIF
        
        IF capable_to_plot_dh NE N_ELEMENTS(data_file_name_dh) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.dh_index.',A,' impossible to plot all data.')"              
        ENDIF

        IF capable_to_plot_tec NE N_ELEMENTS(data_file_name_tec) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.tec_index.',A,' impossible to plot all data.')"              
        ENDIF
;###############################################################################
; Generate the time variables to plot dH time series                           
        H    = FLTARR(file_number*24)                       
        FOR i = 0, N_ELEMENTS(exist_data_file_dh)-1 DO BEGIN
                IF exist_data_file_dh[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date_dh[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        d_dh = DH_teo([tmp_year, tmp_month, tmp_day])
                        
                        H[i*24:(i+1)*24-1] = d_dh.H[*]
                                                                                                                       
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
 
    tec_days= findgen(file_number*12)/12.0   
;############################################################################### 
    tec_diff = tec-med
    
    new_tecdays = findgen(file_number*48)/48.0 ;se genera un arreglo de tiempo con 
;    muestreo cada 15 min. para mejorar la resolución de las gráficas    
    
    new_tecdiff = FLTARR(N_ELEMENTS(new_tecdays))     	
    
;    print, N_ELEMENTS(new_tecdays), N_ELEMENTS(new_tecdiff)
    tmp_tecdif  = INTERPOL(tec_diff, N_ELEMENTS(new_tecdays))
    new_tecdiff = tmp_tecdif
;############################################################################### 
    i_diff = dst-H
    
    new_dstdays = findgen(file_number*96)/96.0 ;se genera un arreglo de tiempo con 
;    muestreo cada 15 min. para mejorar la resolución de las gráficas    
    
    new_idiff = FLTARR(N_ELEMENTS(new_dstdays))     	
    
;    print, N_ELEMENTS(new_tecdays), N_ELEMENTS(new_tecdiff)
    tmp_idiff  = INTERPOL(i_diff, N_ELEMENTS(new_dstdays))
    new_idiff = tmp_idiff    
;###############################################################################
; initiate the figure Device, defining colors and figure dim
        Device_bak = !D.Name 
        SET_PLOT, 'Z'
        
        Xsize=fix(1000)
        Ysize=200
        
    ;    DEVICE, DECOMPOSED = 0
        DEVICE, SET_RESOLUTION = [Xsize,Ysize]
        DEVICE, z_buffering=0
        DEVICE, set_character_size = [10, 12]
             
        chr_size1 = 0.9
        chr_thick1= 1.5
        space     = 0.015
        rojo      = 248
        naranja   = 220
        amarillo  = 198
        verde     = 160
        negro     = 0
        azul      = 100
        blanco    = 255
        gris      = 100
        morado    = 50
        
    TVLCT, R_bak, G_bak, B_bak, /GET
        
    LOADCT, 39, /SILENT
;###############################################################################
; Create a time label based on the DOY initial and DOY final inputs
     X_label   = STRARR(file_number+1)+' '
     ;print, n_elements(x_label)
        months    = ['Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic']
        old_month = mh_i
       ; print, old_month
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
;###############################################################################
; plot the Dst time series for each event
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

window_title = 'TGM'+ string(TGM_n, format='(I01)')+', '+ $
string(old_month, yr_i, format='(A, X, I4)')
;#####################################################################################    
    med_idx = MEDIAN(new_idiff)
    std_idx = stddev(new_idiff, /NAN)
    
    i_out = WHERE(new_idiff GE med_idx+std_idx OR new_idiff LE med_idx-std_idx)
    i_in  = WHERE(new_idiff LE med_idx+std_idx AND new_idiff GE med_idx-std_idx)
    id_diff_out = new_idiff
    id_diff_out[i_in]=!Values.F_NAN
    
    id_diff_in  = new_idiff
    id_diff_in[i_out]=!Values.F_NAN

    sup0 = med_idx+std_idx
    inf0 = med_idx-std_idx
    
    lim_sup = fltarr(n_elements(new_idiff))
    lim_sup[*] = sup0

    lim_inf = fltarr(n_elements(new_idiff))
    lim_inf[*] = inf0   
                    
    up  = max(dst-H)
    down= min(dst-H)
    up2      = max(tec-med)
    down2    = min(tec-med)
    
    if up2 gt up then up0 = up2 else up0 = up
    if down2 lt down then down0 = down2 else down0 = down
       dH = TeXtoIDL('\DeltaH')   
;#####################################################################################  
       days = intarr(file_number+1)
       for i=0, n_elements(days)-1 do begin
            days[i] = dy_i+i
       endfor
       days = days*24/24. 
      day_time = findgen(24)            
      tot_days = fltarr(n_elements(days_dst))      
      for i=0, file_number-1 do begin  
      for j=0, 23 do begin       
      tot_days = days[i]+(day_time[j]/24)                       
      endfor
      endfor
;#####################################################################################         
    plot, new_dstdays, new_idiff, XTICKS=file_number, xminor = 12, POSITION=[0.1,0.18,0.9,0.7],$
    xstyle = 5, ystyle=6, YRANGE=[down,up], XRANGE=[0, file_number],$
    title = '', BACKGROUND = blanco, COLOR=negro,$
	ytitle = 'Indice DST [nT]',  XTICKNAME=REPLICATE(' ', file_number+1), CHARSIZE = 1.2,$
	thick=2, /NODATA

          f_a   = N_ELEMENTS(new_dstdays)
    POLYFILL, [new_dstdays[0], new_dstdays[f_a-1], new_dstdays[f_a-1], new_dstdays[0]],$
    [inf0, inf0, sup0, sup0], COLOR=azul, /LINE_FILL, SPACING=0.001, linestyle=1, THICK=0.1  		
    POLYFILL, [new_dstdays[0], new_dstdays[f_a-1], new_dstdays[f_a-1], new_dstdays[0]],$
    [inf0, inf0, sup0, sup0], COLOR=azul, /LINE_FILL, SPACING=0.001, linestyle=1, ORIENTATION=90, THICK=0.1        

    oplot, new_dstdays, id_diff_in, color=negro, linestyle=3        
    oplot, new_dstdays, id_diff_out, color=negro, linestyle=0, thick=4        
;#####################################################################################        
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
    
    l_sup = fltarr(n_elements(new_tecdiff))
    l_sup[*] = sup

    l_inf = fltarr(n_elements(new_tecdiff))
    l_inf[*] = inf    
;#####################################################################################
        plot, new_tecdays, new_tecdiff, color=rojo, XTICKS=file_number, xminor = 8,$
        xstyle = 5, ystyle=6, YRANGE=[down2,up2], POSITION=[0.1,0.18,0.9,0.7],$
        BACKGROUND = blanco, XRANGE=[0, file_number], /noerase,$
        XTICKNAME=REPLICATE(' ', file_number+1), CHARSIZE = 0.8, thick=3, /NODATA

    POLYFILL, [new_dstdays[0], new_dstdays[f_a-1], new_dstdays[f_a-1], new_dstdays[0]],$
    [inf, inf, sup, sup], COLOR=amarillo, /LINE_FILL, SPACING=0.001, linestyle=1,$
    ORIENTATION=-45  

    POLYFILL, [new_dstdays[0], new_dstdays[f_a-1], new_dstdays[f_a-1], new_dstdays[0]],$
    [inf, inf, sup, sup], COLOR=amarillo, /LINE_FILL, SPACING=0.001, linestyle=1,$
    ORIENTATION=45  	
    
        oplot, new_tecdays, tec_diff_in, color=rojo, linestyle=0
        oplot, new_tecdays, tec_diff_out, color=rojo, linestyle=0, thick=4

    oplot, new_tecdays, l_sup, color=rojo, linestyle=1, thick=1
    oplot, new_tecdays, l_inf, color=rojo, linestyle=1, thick=1


    plot, new_dstdays, new_idiff, XTICKS=file_number, POSITION=[0.1,0.18,0.9,0.7], xminor = 12,$
    xstyle = 5, ystyle=6, YRANGE=[down,up], XRANGE=[0, file_number],$
    title = '', BACKGROUND = blanco, COLOR=negro,$
	ytitle = 'Indice DST [nT]',  XTICKNAME=REPLICATE(' ', file_number+1), CHARSIZE = 1.2,$
	thick=2, /NODATA, /NOERASE

    oplot, new_dstdays, id_diff_in, color=negro, linestyle=0, THICK=1        
    oplot, new_dstdays, id_diff_out, color=negro, linestyle=0, thick=5  
             
        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         xstyle=1, $
                         XTICKS=file_number, $
                         XMINOR=12, $
                         xtitle='Tiempo Universal [Dias]', $
                         XTICKNAME=X_label ,$
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1, $
                         TICKLEN=0.1
                             
        AXIS, XAXIS = 1, XRANGE = (!X.CRANGE+dy_i-0.25),$
                         XTICKS=file_number, $ 
                         XTICKV=days,$
                         xstyle=1, $
                         XMINOR=12, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1, $
                         COLOR=negro, $
                         TICKLEN=0.1

        AXIS, YAXIS = 0, YRANGE=[down0,up0], $
                         YTITLE = 'Dst-'+dH+' [nT]', $                          
                         COLOR=negro, $
                         ystyle=2, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1;, $
                        
        AXIS, YAXIS = 1, YRANGE=[down2,up2], $
                         YTITLE = 'TEC-<TEC> [TECu]', $   
                         COLOR=rojo, $
                         ystyle=2, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1;, $    
                        
   x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
   y = 0.91
   XYOuts, X, y, window_title, /Normal, color=negro, Alignment=0.5, $
    Charsize=1.25, CHARTHICK=1.5

   XYOuts, X, 0.82, 'Tiempo Local [Dias]', /Normal, color=negro, Alignment=0.5, $
    Charsize=0.8, CHARTHICK=1.5   
     Image=TVRD() 
    TVLCT, reds, greens, blues, /get                          ; reads Z buffer !!
    
    TVLCT, R_bak, G_bak, B_bak
        
    ;DEVICE, /CLOSE
    SET_PLOT, Device_bak   
;###############################################################################
; open the post stript device
    path = '../rutidl/output/eventos_tgm/'   
    Date    = string(yr_i, mh_i, dy_i, FORMAT='(I4, "-", I02, "-", I02)')
        IF keyword_set(jpeg) THEN BEGIN
                info = size(Image)
                nx = info[1]
                ny = info[2]
                true_image = bytarr(3,nx,ny)
                true_image[0,*,*] = reds[image]
                true_image[1,*,*] = greens[image]
                true_image[2,*,*] = blues[image]
                write_jpeg, path+'idxV1_'+Date+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'idxV1_'+Date+'.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'mag_idx_'+Date+'_V1.png'
                print, ''
        ENDIF
        RETURN
end    
