readcol, 'M20190510c48.txt', dt, x ; leer las columnas del archivo .txt
;print, dt	;asegurarse de que los datos leídos son los correctos
;print, x

;plot, x


y = fft(x)	; aplicamos la transformada rápida de fourier a los 			;datos de la segunda columna
print, y
;plot, y	;revisamos los datos 
;plot, y[8209:*], yrange = [-3.8e-2, 0.9e-2]<	; ajustamos la escala del rango y solo consideramos
						; la parte real

f_k = findgen(n_elements(x))/(n_elements(x)*dt)	;obtenemos las 							;frecuencias
;print, f_k
;plot, f_k<			; se grafica para revisar los resultados




;plot, y, abs(y)^2, yrange = [0,0.4e2]< ; se obtiene el espectro de 						;potencia
psw = abs(y)^2 
print, psw
;plot, f_k[8209:*], psw[8209:*], yrange = [0, 0.15e-2]	; se evalua la frecuencia contra 									;el espectro de potencia y se 								;justa el rango de valores de los ejes
