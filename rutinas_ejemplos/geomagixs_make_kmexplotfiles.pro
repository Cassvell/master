;+
; NAME:
;       kmex_update_kplot
;
;
; PURPOSE:
;
;       Plot histogram of K_mex (a_mex) index for a given time period
;
; AUTHOR:
;
;       Pedro Corona Romero
;       Instituto de Geofisica, Unidad Michoacan
;       UNAM, Campus Morelia
;       Antigua Carretera a Patzcuaro
;       Morelia, Michoacan, Mexico
;       piter.cr@gmail.com
;       Mexico, 28.iv.mmxvii
;
; CATEGORY:
;
;       Numerical Data Analize
;
; CALLING SEQUENCE:
;
;       kmex_update_kplot, initial [, final, N_DAYS=n_days]
;
;       Description:
;       Make histrogram-plots of K_mex index given a time period. The time
;       period could be definied by the "initial" and "final" dates or by
;       the "initial" date and the number of days (N_DAYS=#).
;
; PARAMETERS:
;       N_DAYS     : number of days for data to plot
;       STATION    : name or code of the used Geomagnetic Station (GMS)
;
; KEYWORD PARAMETERS:
;
;       /QUIET     : set to supress the messages from the program
;       /REAL_TIME : set to use real-time K_mex data for ploting
;       
;
; DEPENDENCIAS:
;       ?????????? : ????????????
;
; ARCHIVOS ANALIZADOS:
;       GMS_YYYYMMDD.index.final[/early]
;
; ARCHIVOS DE SALIDA:
;
; HISTORIA:
;-
;##############################################################################
;## AUXILARY FUNCTIONS AND PROGRAMS SECTION
;##
;## The fuctions and programs are:
;##     * reading_kmex_data (date [,STATION=station,REAL_TIME=real_time,QUIET=quiet])
;##     * makeing_kmex_plot, initial, final[, STATION = station, REAL_TIME = real_time, QUIET = quiet, PNG = png, JPEG = jpeg, DIR=dir, ERROR=error]
;##     * making_kplot, initial, final, N_DAYS=n_days[, STATION=station, QUIET=quiet, REAL_TIME=real_time, PNG = png, JPEG = jpge, DIR = dir, ERROR=error]
;##############################################################################

;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
FUNCTION reading_kmex_data, date, STATION=station, $
                               REAL_TIME=real_time, $
                               QUIET=quiet

        On_error, 2
        COMPILE_OPT idl2

;##############################################################################
; initialize directories
;##############################################################################
        ;kmex_setup
        @geomagixs_commons

        extention       = keyword_set(real_time) ? '.early' : '.final'
        extention       = STATION EQ 'planetary' ? '' : extention
        string_date     = string(date[0], date[1], date[2], FORMAT='(I4,I02,I02)')

        
        file_name       = system.indexes_dir+gms[system.gms].name+'/'+gms[system.gms].code+'_'+string_date+'.k_index'+extention

        file = FILE_SEARCH(file_name, COUNT=opened_files)
                IF opened_files EQ 0 THEN MESSAGE, 'Error finding K-index data files or directories!'
                
        number_of_lines  = FILE_LINES(file)
        k_index_data = STRARR(number_of_lines)

        OPENR, lun, file, /GET_LUN, ERROR=err
                IF err NE 0 THEN MESSAGE, 'Error opening '+file_name[0]
                READF, lun, k_index_data, FORMAT='(A)'
        CLOSE, lun
        FREE_LUN, lun


        result = { $
                  K_mex        : INTARR(8), $
                  K_SUM        : 0, $
                  a_mex        : INTARR(8), $ 
                  A_median     : 0, $
                  K_mex_max    : INTARR(8), $
                  K_SUM_max    : 0, $
                  a_mex_max    : INTARR(8), $ 
                  A_median_max : 0, $
                  K_mex_min    : INTARR(8), $
                  K_SUM_min    : 0, $
                  a_mex_min    : INTARR(8), $ 
                  A_median_min : 0 $
                 }

        
        dat_str = { z : [0,0,0,0,0,0,0,0], y : 0}
        tmp_var = REPLICATE(dat_str,6)
        READS, k_index_data, tmp_var, FORMAT='(I3,X,I3,X,I3,X,I3,X,I3,X,I3,X,I3,X,I3,X,I3)';'(8(I3,X), I3)'
        result.K_mex[*] = tmp_var[0].z
        result.a_mex[*] = tmp_var[1].z
        result.K_SUM    = tmp_var[0].y
        result.A_median = tmp_var[1].y
        result.K_mex_max[*] = tmp_var[2].z
        result.a_mex_max[*] = tmp_var[3].z
        result.K_SUM_max    = tmp_var[2].y
        result.A_median_max = tmp_var[3].y
        result.K_mex_min[*] = tmp_var[4].z
        result.a_mex_min[*] = tmp_var[5].z
        result.K_SUM_min    = tmp_var[4].y
        result.A_median_min = tmp_var[5].y
        RETURN, result



END

;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
PRO making_kmex_plot, initial, final, STATION = station, $
                                    REAL_TIME = real_time, $
                                    QUIET = quiet, $
                                    PNG = png, $
                                    JPEG = jpeg, $
                                    DIR=dir

        On_error, 2
        COMPILE_OPT idl2

;##############################################################################
; initialize directories
;##############################################################################
        ;kmex_setup
        @geomagixs_commons

;##############################################################################
; initializing dates and hours
;##############################################################################
        initial_year   = initial[0]
        initial_month  = initial[1]
        initial_day    = initial[2]

        final_year     = final[0]
        final_month    = final[1]
        final_day      = final[2]



;##############################################################################
; reading data files
;##############################################################################
        file_number    = (JULDAY(final_month, final_day, final_year) - JULDAY(initial_month, initial_day, initial_year))+1
        data_file_name = strarr(file_number)
        ;years          = intarr(file_number)
        ;months         = intarr(file_number)
        ;days           = intarr(file_number)
        string_date     = strarr(file_number)
        ;exist_file     = intarr(file_number)
        extention       = keyword_set(real_time) ? '.early' : '.final'
        extention       = STATION EQ 'planetary' ? '' : extention
       
        FOR i=0ll, file_number-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(initial_month, initial_day, initial_year)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                ;years[i]  = tmp_year
                ;months[i] = tmp_month
                ;days[i]   = tmp_day
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                ;data_file_name[i] = station_code+'_'+string_date[i]+'.dat'
                data_file_name[i] = gms[system.gms].code+'_'+string_date[i]+'.k_index'+extention
                ;kmex_file_name[i] = station_code+'_'+string_date[i]+'.dat'
                ;print, i+1, '  ', data_file_name[i]
        ENDFOR

        exist_data_file   = FILE_TEST(system.indexes_dir+gms[system.gms].name+'/'+data_file_name)
        capable_to_plot   = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_plot NE N_ELEMENTS(data_file_name) AND NOT keyword_set(quiet) THEN BEGIN 
                PRINT, FORMAT="('CRITICAL ERROR: impossible to read data file(s).')"
                PRINT, FORMAT="('                missing GMS_YYYYMMDD.k_index.',A,' impossible to plot all data.')", extention
                Error.value[1] += 1
                error.log      += 'Impossible to read index k_index-file; file not found or reading permission conflict. '
        ENDIF

        k_mex_data    = INTARR(file_number*8)
        k_SUM_data    = INTARR(file_number)
        a_mex_data    = INTARR(file_number*8)
        a_median_data = INTARR(file_number)


        FOR i = 0, N_ELEMENTS(exist_data_file)-1 DO BEGIN
                IF exist_data_file[i] EQ 1 THEN BEGIN
                        tmp_year    = 0
                        tmp_month   = 0
                        tmp_day     = 0
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        ;print, tmp_year, tmp_month, tmp_day
                        tmp_data    = reading_kmex_data([tmp_year, tmp_month, tmp_day], QUIET=quiet, STATION=station, REAL_TIME=real_time)
                        k_mex_data[i*8:(i+1)*8-1] = tmp_data.K_mex[*]/10
                        a_mex_data[i*8:(i+1)*8-1] = tmp_data.a_mex[*]
                        k_SUM_data[i]             = tmp_data.K_SUM
                        a_median_data[i]          = tmp_data.A_median
                        ;print, string_date[i], k_mex_data[i*8:(i+1)*8-1]
                        ;print, string_date[i], a_mex_data[i*8:(i+1)*8-1]
                        ;print, tmp_data.K_mex[*]
                ENDIF ELSE BEGIN
                        k_mex_data[i*8:(i+1)*8-1] = 99
                        a_mex_data[i*8:(i+1)*8-1] = 999
                        k_SUM_data[i]             = 999
                        a_median_data[i]          = 999
                ENDELSE
        ENDFOR


;##############################################################################
; making plot
;##############################################################################
        
        
        IF gms[system.gms].name EQ 'planetary' THEN BEGIN
                LANCE_banner  = 'https://www.gfz-potsdam.de/en/kp-index/'
                tmp_string    = 'by GFZ German Research Center for Geosciencies'

                MAG_banner    = 'Kp: '+tmp_string
                ;MAG_banner    = MAG_banner+' (LAT '+STRING(kmex_GMS_latitude[kmex_GMS],FORMAT='(F7.2)')+', LON '+STRING(kmex_GMS_longitude[kmex_GMS],FORMAT='(F7.2)')+')'
                Time_title    = 'UTC [Begin: '+string(initial_year, initial_month, initial_day, FORMAT='(I4,"/",I02,"/",I02)')+' 00:00 UTC]'
                K_axis_title  = 'Kp';'Kmex!N index'
                Plot_title    = ' Planetary K index'; +' (GMS: '+kmex_GMS_name[kmex_GMS]+')'
                Plot_title    = keyword_set(real_time) ? 'Nowcast'+Plot_title : 'Definitive'+Plot_title
        ENDIF ELSE BEGIN
                LANCE_banner  = 'LANCE/SCiESMEX - Morelia, Mich., MX'
                tmp_string    = (gms[system.gms].name EQ 'teoloyucan') ? ' Geomagnetic Observatory' : ' Geomagnetic Station'

                MAG_banner    = STRUPCASE(gms[system.gms].code)+': '+STRUPCASE(STRMID(gms[system.gms].name,0,1))+STRMID(gms[system.gms].name,1,STRLEN(gms[system.gms].name))+tmp_string
                MAG_banner    = MAG_banner+' (LAT '+STRING(gms[system.gms].latitude,FORMAT='(F7.2)')+', LON '+STRING(gms[system.gms].longitude,FORMAT='(F7.2)')+')'
                Time_title    = 'UTC [Begin: '+string(initial_year, initial_month, initial_day, FORMAT='(I4,"/",I02,"/",I02)')+' 00:00 UTC]'
                K_axis_title  = STRUPCASE(gms[system.gms].code);'Kmex!N index'
                Plot_title    = ' Local K index'; +' (GMS: '+kmex_GMS_name[kmex_GMS]+')'
                Plot_title    = keyword_set(real_time) ? 'Estimated'+Plot_title : 'Calculated'+Plot_title
        ENDELSE
        
        
        
        Julian = SYSTIME(/JULIAN, /UTC) 
        tmp_year    = 0
        tmp_month   = 0
        tmp_day     = 0
        tmp_hour    = 0
        tmp_minute  = 0
        CALDAT, SYSTIME(/JULIAN, /UTC), tmp_month, tmp_day, tmp_year, tmp_hour, tmp_minute
        UPDATE_banner = 'Updated: '+STRING(tmp_year,tmp_month,tmp_day, tmp_hour,tmp_minute, FORMAT='(I4,"/",I02,"/",I02,"-",I02,":",I02," UTC")')
        
        
;print, time
        ;time          = findgen(file_number*8)/DOUBLE(8)+JULDAY(initial_month,initial_day,initial_year)
        
        Device_bak = !D.Name 
        SET_PLOT, 'Z'
        
        tmp_spam = 1
        IF file_number GT 7 THEN tmp_spam = 1.5
        IF file_number GT 15 THEN tmp_spam = 2.
        
        Xsize=fix(800*tmp_spam)
        Ysize=400
        ;DEVICE, decompose=0
        DEVICE, SET_RESOLUTION = [Xsize,Ysize]
        DEVICE, z_buffering=O
        DEVICE, set_character_size = [10, 12]
        ;window, 1, xsize=Xsize, ysize=Ysize
        
        
        time          = findgen(file_number*8)/8.0;+JULDAY(initial_month,initial_day,initial_year)
        ;print, time

        ;oplot, time+JULDAY(initial_month,initial_day,initial_year), k_mex_data, psym=6

;goto, jump1


        chr_size1 = 0.9
        chr_thick1= 1.0
        space     = 0.015
        rojo      = 248
        amarillo  = 198
        verde     = 160
        negro     = 0
        blanco    = 255
        gris      = 150
        morado    = 16
        
        
        TVLCT, R_bak, G_bak, B_bak, /GET
        
        LOADCT, 39, /SILENT

;!P.FONT = 0

        X_label   = STRARR(file_number+1)+' '
        months    = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        old_month = initial_month
        FOR i =0,  N_ELEMENTS(X_label)-1 DO BEGIN
                tmp_year    = 0
                tmp_month   = 0
                tmp_day     = 0
                tmp_julday  = JULDAY(initial_month, initial_day, initial_year)

                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                
                IF i LT N_ELEMENTS(X_label)-1 THEN $
                        X_label[i]  = (tmp_month EQ old_month) ? string(tmp_day, FORMAT='(I02)') : months[tmp_month-1]+' '+string(tmp_day, FORMAT='(I02)')
                old_month = tmp_month
                ;print, X_label[i]

        ENDFOR

        plot, time, k_mex_data, psym=6, /NoDATA, BACKGROUND = blanco, COLOR=negro, MAX_VALUE=9., $
                                YRANGE=[0,9], YSTYLE=5, YTICKS=9, YMINOR=0, $
                                TITLE=Plot_title, $
                                SUBTITLE= Time_title,$
                                YTITLE = K_axis_title, $
                                XRANGE=[0,file_number], XSTYLE=5, XTICKS=file_number, XMINOR=8, $; XTICKUNITS = 'Days', $
                                XTITLE = Time_title, $; XTICKUNITS = 'Time', XTICKFORMAT='LABEL_DATE',$
                                XTICKNAME=REPLICATE(' ', file_number+1), $
                                CHARSIZE = chr_size1, CHARTHICK=chr_thick1, $
                                POSITION=[0.07,0.35,0.95,0.9];, font=1
                                

        j = N_ELEMENTS(k_mex_data)-1
        ;WHILE k_mex_data[j] GT 9 AND j GT 0 DO j--
        
        
        ;print, j
        
        ;POLYFILL, [0.,0.125,0.125,0.]+time[0], [0,0,k_mex_data[0],k_mex_data[0]], color=255
        ;FOR i = 0, N_ELEMENTS(k_mex_data)-1 DO BEGIN
        FOR i = 0, j DO BEGIN
                IF k_mex_data[i] LE 9 THEN BEGIN
                                        ;LOADCT, 13, /SILENT
                                        color = 0
                                        step  = (k_mex_data[i] EQ 0) ? 0.1 : 0.
                                        CASE 1 OF
                                                k_mex_data[i] EQ 4 : color = amarillo
                                                k_mex_data[i] GT 4 : color = rojo
                                                ELSE               : color = verde
                                        ENDCASE
                                        POLYFILL, [0.+space,0.125-space,0.125-space,0.+space]+time[i], [0,0,k_mex_data[i]+step,k_mex_data[i]+step], color=color
                                      ENDIF $
                                      ELSE BEGIN
                                        ;LOADCT, 0, /SILENT
                                        ;POLYFILL, [0.+space,0.125-space,0.125-space,0.+space]+time[i], [0,0,9.,9.], color=gris, /LINE_FILL, ORIENTATION=45.
                                        ;POLYFILL, [0.+space,0.125-space,0.125-space,0.+space]+time[i], [0,0,9.,9.], color=gris, /LINE_FILL, ORIENTATION=-45.
                                        POLYFILL, [0.+space,0.125-space,0.125-space,0.+space]+time[i], [0,0,9.,9.], color=morado, /LINE_FILL, ORIENTATION=45., linestyle = 0
                                        POLYFILL, [0.+space,0.125-space,0.125-space,0.+space]+time[i], [0,0,9.,9.], color=morado, /LINE_FILL, ORIENTATION=-45., linestyle = 0
                                      ENDELSE
        ENDFOR

        ;LOADCT, 0, /SILENT
        FOR i = 0, file_number-1 DO BEGIN
                OPLOT, [i,i], [0.,9.], linestyle=1, COLOR=negro
                ;print, [i,i], [0.,9.]
        ENDFOR

        FOR i=5, 8 DO OPLOT, [0,file_number], [1.,1.]*i, linestyle=1, COLOR=negro
        ;OPLOT, [0,file_number], [6.,6.], linestyle=1, COLOR=negro

        AXIS, XAXIS = 0, XRANGE=[0,file_number], $
                         XTICKS=file_number, $
                         XMINOR=8, $
                         ;XTITLE = ' ' , $;Time_title, $; XTICKUNITS = 'Time', XTICKFORMAT='LABEL_DATE',$
                         XTICKNAME=X_label, $
                         COLOR=negro, $
                         CHARSIZE = chr_size1, $
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
                         YTITLE = K_axis_title, $
                         COLOR=negro, $
                         CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1;, $
                         ;TICKLEN=0.00

        AXIS, YAXIS = 1, YRANGE=[0,9], $
                         YTICKS=9, $
                         YMINOR=1, $
                         YTICKNAME=[' ', ' ', ' ', ' ', ' ', 'G1', 'G2', 'G3', 'G4', 'G5'], $
                         COLOR=negro, $
                         CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1;, $


        XYOUTS, 0.01, 0.015 , /NORMAL, $
                LANCE_banner, COLOR=negro, $
                CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1

        XYOUTS, 0.01, 0.070 , /NORMAL, $
                MAG_banner, COLOR=negro, $
                CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1

        XYOUTS, 0.01, 0.165 , /NORMAL, $
                'Color Code:        quiet,        disturbed,        storm,          data not available.', COLOR=negro, $
                CHARSIZE = chr_size1, $
                CHARTHICK=chr_thick1
        POLYFILL, [0.15,0.18,0.18,0.15], [0.165,0.165,0.195,0.195], color = verde, /NORMAL
        POLYFILL, [0.27,0.30,0.30,0.27], [0.165,0.165,0.195,0.195], color = amarillo, /NORMAL
        POLYFILL, [0.44,0.47,0.47,0.44], [0.165,0.165,0.195,0.195], color = rojo, /NORMAL
        POLYFILL, [0.58,0.63,0.63,0.58], [0.165,0.165,0.195,0.195], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=45., linestyle = 0
        POLYFILL, [0.58,0.63,0.63,0.58], [0.165,0.165,0.195,0.195], color = morado, /NORMAL, /LINE_FILL, ORIENTATION=-45., linestyle = 0

        XYOUTS, 0.65, 0.015 , /NORMAL, $
                UPDATE_banner, COLOR=negro, $
                CHARSIZE = chr_size1, $
                         CHARTHICK=chr_thick1


        Image=TVRD()                           ; reads Z buffer !!
        TVLCT, reds, greens, blues, /get
        

        TVLCT, R_bak, G_bak, B_bak

		;DEVICE, /CLOSE

        
        SET_PLOT, Device_bak
        
        file_name = gms[system.gms].code+'_'+string(initial_year,initial_month,initial_day,file_number, FORMAT='(I4,I02,I02,"_",I02,"days.k")')+extention
        ;file_dir  = kmex_plot_dir
        file_dir  = system.plots_dir+gms[system.gms].name+'/'
        CASE 1 OF
                file_number EQ 1 : file_dir += '1day/'
                file_number EQ 3 : file_dir += '3days/'
                file_number EQ 7 : file_dir += '7days/'
                ELSE             : file_dir += 'etc/'
        ENDCASE
         
        IF N_ELEMENTS(dir) GT 0 THEN BEGIN
                IF STRLEN(dir) GT 0 THEN BEGIN
                        IF STRPOS(dir, '/', /reverse_search)+1 EQ STRLEN(dir) THEN file_dir = dir $
                                ELSE BEGIN
                                        PRINT, '        It is mandatory a "/" at the end of the given directory!'
                                        IF NOT keyword_set(quiet) THEN PRINT, '        Setting default directory '+file_dir
                                ENDELSE
                ENDIF
        ENDIF
        
        IF NOT FILE_TEST(file_dir, /DIRECTORY) THEN MESSAGE, '        Impossible to save the graph. Directory '+file_dir+' is missing!'
        

        IF keyword_set(jpeg) THEN BEGIN
                info = size(Image)
                nx = info[1]
                ny = info[2]
                true_image = bytarr(3,nx,ny)
                true_image[0,*,*] = reds[image]
                true_image[1,*,*] = greens[image]
                true_image[2,*,*] = blues[image]
                write_jpeg, file_dir+file_name+'.jpg', True_Image, true=1
        ENDIF ELSE BEGIN
                IF NOT (keyword_set(quiet) OR keyword_set(png)) THEN print, '        Setting PNG as default file type.'
                WRITE_PNG, file_dir+file_name+'.png', Image, reds,greens,blues
        ENDELSE

        IF NOT keyword_set(quiet) THEN BEGIN
                print, '        Saving: '+file_name
                print, ''
        ENDIF

;print, file_number
;print, k_mex_data

        RETURN
END

;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
PRO making_kplot, initial, final, N_DAYS=n_days, STATION=station, QUIET=quiet, REAL_TIME=real_time, $
                                PNG = png, $
                                JPEG = jpge, $
                                DIR = dir

        On_error, 2
        COMPILE_OPT idl2, hidden
;##############################################################################
; initialize directories
;##############################################################################
        ;kmex_setup
        @geomagixs_commons
        
;##############################################################################
; depuring inputs
;##############################################################################
        CASE 1 OF
                (N_ELEMENTS(n_days) NE 1) AND $
                (N_PARAMS() EQ 1) $
                                        : BEGIN
                                                MESSAGE, 'It is mandatory a given integer (>= 0) value for N_DAYS!.'
                                                RETURN
                                          END
                (N_ELEMENTS(n_days) EQ 1) AND $
                (N_PARAMS() EQ 1) $
                                        : BEGIN
                                                IF n_days LE 0 THEN BEGIN
                                                        MESSAGE, 'It is mandatory N_DAYS > 0!.'
                                                        RETURN
                                                ENDIF
                                                
                                                IF N_ELEMENTS(initial) LT 3 THEN MESSAGE, 'Error in initial inputs! initial=[YYYY, MM, DD]'
                                                IF N_ELEMENTS(initial) EQ 4 THEN MESSAGE, 'Error in initial inputs! initial=[YYYY, MM, DD, hh, mm]'
                                                IF N_ELEMENTS(initial) GT 5 THEN MESSAGE, 'Error in initial inputs! initial=[YYYY, MM, DD, hh, mm]'

                                                tmp_day   = 0
                                                tmp_month = 0
                                                tmp_year  = 0
                                                CALDAT, JULDAY(initial[1],initial[2],initial[0])-n_days+1, tmp_month, tmp_day, tmp_year
                                                final   = initial
                                                initial = [tmp_year,tmp_month,tmp_day]
                                          END
                ELSE : BEGIN
                        IF N_PARAMS() GT 2 THEN MESSAGE, 'It is mandatory initial and final dates as inputs!'
        
                        IF N_ELEMENTS(initial) LT 3 THEN MESSAGE, 'Error in initial inputs! initial=[YYYY, MM, DD]'
                        IF N_ELEMENTS(final) LT 3 THEN MESSAGE, 'Error in final inputs! final=[YYYY, MM, DD]'

                        IF N_ELEMENTS(initial) EQ 4 THEN MESSAGE, 'Error in initial inputs! initial=[YYYY, MM, DD, hh, mm]'
                        IF N_ELEMENTS(final) EQ 4 THEN MESSAGE, 'Error in final inputs! final=[YYYY, MM, DD, hh, mm]'

                        IF N_ELEMENTS(initial) GT 5 THEN MESSAGE, 'Error in initial inputs! initial=[YYYY, MM, DD, hh, mm]'
                        IF N_ELEMENTS(final) GT 5 THEN MESSAGE, 'Error in final inputs! final=[YYYY, MM, DD, hh, mm]'
                END
        ENDCASE

        initial_year   = initial[0]
        initial_month  = initial[1]
        initial_day    = initial[2]

        final_year     = final[0]
        final_month    = final[1]
        final_day      = final[2]

        initial_JUL = JULDAY(initial_month, initial_day, initial_year)
        final_JUL   = JULDAY(final_month, final_day, final_year)
        
        lower_JUL   = JULDAY(gms[system.gms].dates_index[0,1],gms[system.gms].dates_index[0,2],gms[system.gms].dates_index[0,0])
        upper_JUL   = JULDAY(gms[system.gms].dates_index[1,1],gms[system.gms].dates_index[1,2],gms[system.gms].dates_index[1,0])
        



        
        IF initial_JUL GE lower_JUL AND final_JUL LE upper_JUL THEN $
                making_kmex_plot, [initial_year, initial_month, initial_day], [final_year, final_month, final_day],  $
                                STATION=station, QUIET=quiet, REAL_TIME=real_time, PNG = png, JPEG = jpeg, DIR=dir $
        ELSE BEGIN
                Error.value[3] += 1

                making_kmex_plot, [initial_year, initial_month, initial_day], [final_year, final_month, final_day],  $
                                STATION=station, QUIET=quiet, REAL_TIME=real_time, PNG = png, JPEG = jpeg, DIR=dir

        ENDELSE
        
        RETURN

END




;SECUENCIA PRINCIPAL
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
;##############################################################################
PRO geomagixs_make_kmexplotfiles, initial, final, STATION=station, $
                                       QUIET=quiet, $
                                       FORCE_ALL=force_all, $
                                       REAL_TIME=real_time
        
        On_error, 2
        COMPILE_OPT idl2, HIDDEN
;##############################################################################
; initialize directories
;##############################################################################
        @geomagixs_commons
        geomagixs_setup_commons, QUIET=quiet
        geomagixs_check_system, QUIET=quiet
        geomagixs_setup_dates, STATION=station, QUIET=quiet, /FORCE_ALL

;##############################################################################
; depuring inputs
;##############################################################################
        geomagixs_check_gms, STATION=station, QUIET=quiet

        CASE 1 OF
                N_PARAMS() EQ 0 : BEGIN
                                        initial = system.today_date
                                        final   = initial
                                  END
                N_PARAMS() EQ 1 : final = initial
                N_PARAMS() EQ 2 : geomagixs_check_dates, initial, final, STATION=station, QUIET=quiet
                ELSE            : BEGIN
                                        IF NOT keyword_set(quiet) THEN BEGIN
                                                PRINT, FORMAT="('CRITICAL ERROR: conflict with input date(s), data inconsistent or invalid.')"
                                                PRINT, FORMAT="('                impossible to proceed with an ambiguoss range of dates.')"
                                        ENDIF
                                        error.value[2] += 1
                                        error.log      += 'Ambiguous range of dates, it is required only two input dates to define a time period. '
                                        RETURN
                                  END
        ENDCASE


;##############################################################################
; initializing dates and hours
;##############################################################################
        initial_year   = initial[0]
        initial_month  = initial[1]
        initial_day    = initial[2]

        final_year     = final[0]
        final_month    = final[1]
        final_day      = final[2]




;##############################################################################
; reading data files
;##############################################################################
        plot_number     = (JULDAY(final_month, final_day, final_year) - JULDAY(initial_month, initial_day, initial_year))+1
        ;days_ina_week   = 7
        data_file_name = strarr(plot_number)
        ;years           = intarr(file_number)
        ;months          = intarr(file_number)
        ;days            = intarr(file_number)
        string_date     = strarr(plot_number)
        ;exist_file     = intarr(file_number)
        extention       = keyword_set(real_time) ? '.early' : '.final'
        extention       = STATION EQ 'planetary' ? '' : extention
       
        tmp_year    = 0
        tmp_month   = 0
        tmp_day     = 0
        tmp_julday  = JULDAY(initial_month, initial_day, initial_year)
        FOR i=0ll, plot_number-1 DO BEGIN
                CALDAT, tmp_julday+i, tmp_month, tmp_day, tmp_year
                ;years[i]  = tmp_year
                ;months[i] = tmp_month
                ;days[i]   = tmp_day
                string_date[i]    = string(tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)')
                ;data_file_name[i] = station_code+'_'+string_date[i]+'.dat'
                data_file_name[i] = gms[system.gms].code+'_'+string_date[i]+'.k_index'+extention
                ;kmex_file_name[i] = station_code+'_'+string_date[i]+'.dat'
                ;print, i+1, '  ', data_file_name[i]
        ENDFOR

        exist_data_file   = FILE_TEST(system.indexes_dir+gms[system.gms].name+'/'+data_file_name)

        capable_to_update = N_ELEMENTS(where(exist_data_file EQ 1))

        IF capable_to_update LT 1 THEN BEGIN
                IF not keyword_set(quiet) THEN BEGIN
                        PRINT, '        Data File Error: GMS_YYYYMMDD.k_index.'+extention+' files not found!'
                        PRINT, '        If proceed, plots will be filled with data gaps.'
                ENDIF
                proceed = 'Y'
                REPEAT BEGIN
                        IF not ( keyword_set(quiet) OR keyword_set(force_all) ) THEN READ, proceed, PROMPT='        Continue (Y/N)?: '
                        proceed=STRUPCASE(proceed)
                        IF proceed EQ 'N' OR proceed EQ 'NO' THEN RETURN
                ENDREP UNTIL proceed EQ 'Y' OR proceed EQ 'YES'
        ENDIF

        ;exist_result_file = FILE_TEST(kmex_processed_dir+station_code+'_'+string_date+'.data'+extention) AND not(keyword_set(force_all))
        exist_plot_file = FILE_TEST(system.plots_dir+gms[system.gms].name+'/'+'1day/'+gms[system.gms].code+'_'+string_date+'.01days.k'+extention+'.png') AND not(keyword_set(force_all))

        make_update_file = exist_data_file AND NOT(exist_plot_file)

        IF TOTAL(make_update_file) GT 0 THEN files_to_update = N_ELEMENTS(where(make_update_file EQ 1)) $
                                        ELSE files_to_update = 0
        ;print, where(make_update_file EQ 1)
        
        IF not (keyword_set(quiet) OR keyword_set(force_all) ) THEN BEGIN
                IF TOTAL(files_to_update) GT 0 THEN BEGIN
                        print, ''
                        PRINT, files_to_update, FORMAT='("        A total of ",I," graph(s) need to be updated.")' 
                ENDIF ELSE BEGIN
                        print, ''
                        PRINT, "        No graph requires to be updated."
                        RETURN
                ENDELSE

                IF capable_to_update GT files_to_update THEN BEGIN
                        PRINT, capable_to_update-files_to_update, FORMAT='("        There are still ",I," graph(s) that can be updated.")'
                        PRINT, '        Use the /FORCE_ALL keyword to force the updating of all graph(s).'
                        PRINT, ''
                ENDIF

                IF N_ELEMENTS(exist_data_file) GT capable_to_update THEN BEGIN
                        print, ''
                        PRINT, N_ELEMENTS(exist_data_file) - capable_to_update, FORMAT='("        There are ",I," graph(s) that are unable to be updated.")'
                        PRINT, '        Use the MAKE_INDEXFILES tool for updating K-index data base.'
                        PRINT, ''
               ENDIF
        ENDIF

        proceed = 'Y'
        REPEAT BEGIN
                IF not ( keyword_set(quiet) OR keyword_set(force_all) ) THEN READ, proceed, PROMPT='        Continue (Y/N)?: '
                proceed=STRUPCASE(proceed)
                IF proceed EQ 'N' OR proceed EQ 'NO' THEN RETURN
        ENDREP UNTIL proceed EQ 'Y' OR proceed EQ 'YES'

        tmp_year    = 0
        tmp_month   = 0
        tmp_day     = 0
        FOR i = 0ll, N_ELEMENTS(make_update_file)-1 DO BEGIN
                IF make_update_file[i] EQ 1 THEN BEGIN
                        READS, string_date[i], tmp_year, tmp_month, tmp_day, FORMAT='(I4,I02,I02)'
                        ;print, tmp_year, tmp_month, tmp_day
                        ;kmex_making_magneticdata, [tmp_year, tmp_month, tmp_day], QUIET=quiet, STATION=station, REAL_TIME=real_time
                        making_kplot, [tmp_year, tmp_month, tmp_day], N_DAYS=1, STATION=station, QUIET=quiet, REAL_TIME=real_time
                        making_kplot, [tmp_year, tmp_month, tmp_day], N_DAYS=3, STATION=station, QUIET=quiet, REAL_TIME=real_time
                        making_kplot, [tmp_year, tmp_month, tmp_day], N_DAYS=7, STATION=station, QUIET=quiet, REAL_TIME=real_time
                ENDIF
        ENDFOR

        IF not keyword_set(quiet) THEN BEGIN
                PRINT, ''
                PRINT, '        Kmex-Plot files updated!'
        ENDIF

        RETURN


END
