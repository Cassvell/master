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
		
		file_name = '../master_thesis/datos/teoloyucan/teoloyucan/'+$
		'teo_'+date+'.clean.dat'

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

pro desv_est_DH, date_i, date_f

	On_error, 2
	compile_opt idl2, HIDDEN

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
                
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, $
                FORMAT='(I4,I02,I02)')
                
                data_file_name[i] = '../master_thesis/datos/teoloyucan/teoloyucan/'$
                +'teo_'+string_date[i]+'.clean.dat'
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
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, $
                        FORMAT='(I4,I02,I02)'
                        
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

  ;  time_w = n_elements(time)
  ;  tot_days= findgen(time_w)/24.0 

    path = '../rutidl/output/globfig_to_reg/'

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
                        ;X_label[i]  = (tmp_month EQ old_month) ? string(tmp_day, FORMAT='(I02)') : months[tmp_month-1]+' '+string(tmp_day, FORMAT='(I02)')
                        X_label[i]  = (tmp_month EQ old_month) ? string(tmp_day, FORMAT='(I02)') : string(tmp_day, FORMAT='(I02)')
                old_month = tmp_month
        ENDFOR        

     down0 = min(H_STDESV)   
     up0   = max(H_STDESV)  
      ;print, max(H_STDESV), min(H_STDESV)   
    sigma = TeXtoIDL('\sigmaH')
    dH    = TeXtoIDL('\DeltaH')
    title = 'Desviación estandard ('+dH+') para Octubre y Noviembre, 2003'                      
    plot, time, H_STDESV, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=rojo,$
     CHARSIZE = 1.0, CHARTHICK=chr_thick1, POSITION=[0.07,0.15,0.95,0.90], $
     XSTYLE = 5, XTICKNAME=REPLICATE(' ', file_number+1), XRANGE=[0, file_number], $
     ySTYLE = 6, YRANGE=[down0,up0]

    plot, time, H_STDESV, /NoDATA, XTICKS=file_number, xminor=8, BACKGROUND = blanco, COLOR=negro,$
     CHARSIZE = 1.0, CHARTHICK=chr_thick1, POSITION=[0.07,0.15,0.95,0.90], $
     XSTYLE = 5, XTICKNAME=REPLICATE(' ', file_number+1), XRANGE=[0, file_number], $
     ySTYLE = 6, YRANGE=[down0,up0], title=title, /noerase    
     ;##########################################################################
     ;Qdays delimitation, TGM1
     ;Qday1
     ;oplot, [7,7], [0,300], linestyle=3, color=negro
     ;oplot, [8,8], [0,300], linestyle=3, color=negro

     ;Qday2
     ;oplot, [file_number-3,file_number-3], [0,300], linestyle=3, color=negro
     ;oplot, [file_number-4,file_number-4], [0,300], linestyle=3, color=negro 
     
     
     ;TGM days delimitation          
     ;oplot, [27,27], [0,300], linestyle=5, color=negro
     ;oplot, [23,23], [0,300], linestyle=5, color=negro
     
     ;oplot, [file_number-7,file_number-7], [0,300], linestyle=5, color=negro
     ;oplot, [file_number-11,file_number-11], [0,300], linestyle=5, color=negro         
     ;##########################################################################
     ;Qdays delimitation, TGM2
     ;Qday1
     ;oplot, [5,5], [0,300], linestyle=3, color=negro
     ;oplot, [5.97,5.97], [0,300], linestyle=3, color=negro

     ;Qday2
     ;oplot, [file_number-5,file_number-5], [0,300], linestyle=3, color=negro
     ;oplot, [file_number-6,file_number-6], [0,300], linestyle=3, color=negro 
     
     
     ;TGM days delimitation          
     ;oplot, [6.03,6.03], [0,300], linestyle=5, color=negro
     ;oplot, [12,12], [0,300], linestyle=5, color=negro               
     ;##########################################################################
     ;Qdays delimitation, TGM3
     ;Qday1
     ;oplot, [5,5], [0,300], linestyle=3, color=negro
     ;oplot, [4,4], [0,300], linestyle=3, color=negro

     ;Qday2
     ;oplot, [file_number-5,file_number-5], [0,300], linestyle=3, color=negro
     ;oplot, [file_number-6,file_number-6], [0,300], linestyle=3, color=negro 
     
     
     ;TGM days delimitation          
     ;oplot, [13.8,13.8], [0,300], linestyle=5, color=negro
     ;oplot, [16,16], [0,300], linestyle=5, color=negro       
     ;########################################################################## 
     ;Qdays delimitation, TGM4
     ;Qday1
     ;oplot, [2,2], [0,300], linestyle=3, color=negro
     ;oplot, [3,3], [0,300], linestyle=3, color=negro

     ;Qday2
     ;oplot, [file_number-3,file_number-3], [0,300], linestyle=3, color=negro
     ;oplot, [file_number-4,file_number-4], [0,300], linestyle=3, color=negro 
     
     
     ;TGM days delimitation          
     ;oplot, [8,8], [0,300], linestyle=5, color=negro
     ;oplot, [16,16], [0,300], linestyle=5, color=negro         
     ;##########################################################################
     ;Qdays delimitation, TGM5
     ;Qday1
     ;oplot, [6,6], [0,400], linestyle=3, color=negro
     ;oplot, [6.97,6.97], [0,400], linestyle=3, color=negro

     ;Qday2
     ;oplot, [16,16], [0,400], linestyle=3, color=negro
     ;oplot, [15,15], [0,400], linestyle=3, color=negro 


     ;TGM days delimitation          
     ;oplot, [7.03,7.03], [0,400], linestyle=5, color=negro
     ;oplot, [11,11], [0,400], linestyle=5, color=negro              
     ;##########################################################################   
     ;Qdays delimitation, TGM6
     ;Qday1
     oplot, [6,6], [0,400], linestyle=3, color=negro
     oplot, [7,7], [0,400], linestyle=3, color=negro

     ;Qday2
     oplot, [file_number-4,file_number-4], [0,400], linestyle=3, color=negro
     oplot, [file_number-5,file_number-5], [0,400], linestyle=3, color=negro 
     

     ;TGM days delimitation          
     oplot, [file_number-25,file_number-25], [0,400], linestyle=5, color=negro
     oplot, [file_number-20,file_number-20], [0,400], linestyle=5, color=negro              
     ;########################################################################## 

        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         ;XTITLE = ' ' , $;Time_title, $; XTICKUNITS = 'Time', XTICKFORMAT='LABEL_DATE',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         xtitle='Time [Days]'
                         CHARTHICK=chr_thick1
                         
        AXIS, XAXIS = 1, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         XTICKNAME=REPLICATE(' ', file_number+1), $
                         COLOR=negro, $
                         TICKLEN=0.04

        AXIS, YAXIS = 0, YRANGE=[down0,up0],$
                         Ystyle=2, $
                         ;YMINOR=1, $
                         YTITLE = sigma, $
                         COLOR=negro, $
                         CHARSIZE = 0.9, $
                         CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00

        AXIS, YAXIS = 1, YRANGE=[down0,up0],$
                         Ystyle=2, $
                         ;YMINOR=1, $
                         ;YTICKNAME=[' ', ' ', ' ', ' ', ' ', 'G1', 'G2', 'G3', 'G4', 'G5'], $
                         COLOR=negro, $
                         CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1;, $
                         
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
                write_jpeg, path+'DH_stdesv'+fecha+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, path+'DH_stdesv'+fecha+'.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+path+'DH_stdesv'+fecha
                print, ''
        ENDIF
        RETURN    
end

