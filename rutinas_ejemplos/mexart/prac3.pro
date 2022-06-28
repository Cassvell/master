readcol, 'M20190514c48.txt', dt, x ; leer las columnas del archivo .txt
;print, dt	;asegurarse de que los datos leídos son los correctos
;print, x

;plot, x


y = fft(x)	; aplicamos la transformada rápida de fourier a los 			;datos de la segunda columna
;print, y
;plot, y	;revisamos los datos 
;plot, y[9734:*], yrange = [-6e-2, 0.5e-2]	; ajustamos la escala del rango
;help, y
f_k = findgen(n_elements(x))/(n_elements(x)*dt)	;obtenemos las 							;frecuencias
;print, f_k
;plot, f_k			; se grafica para revisar los resultados




;plot, y, abs(y)^2, yrange = [0,0.3e2] ; se obtiene el espectro de 						;potencia
psw = abs(y)^2
plot, f_k[9734:*], psw[9734:*], yrange = [0, 0.3e-2]	; se evalua la frecuencia contra el espectro de potencia y se ajusta el rango de valores de los ejes
