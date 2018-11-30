;
; loads film and makes a tif image from it
;
; hazen 12/98
;

loadct, 5

COMMON colors, R_ORIG, G_ORIG, B_ORIG, R_CURR, G_CURR, B_CURR

; get file to open

run = "asdf"
print, "name of file to analyze"
read, run

film_x = fix(1)
film_y = fix(1)

;;; input film
close, 1				; make sure unit 1 is closed

; open file
;dir=''
dir = "C:\Users\Ha-Myong\Documents\MATLAB\"
openr, 1, dir + run + ".pma"

; figure out size + allocate appropriately

result = FSTAT(1)
readu, 1, film_x
readu, 1, film_y
film_l = long(long(result.SIZE-4)/(long(film_x)*long(film_y)))
;film_l = long(long(result.SIZE-4)/(long(film_x)*long(film_y))) - 50

;*****delete me************
;film_l=11000
;*************

print, "film x,y,l : ", film_x,film_y,film_l
;film_l = 100		; for debugging
;if film_l gt 60 then film_l = 60

window, 0, xsize=film_x, ysize=film_y
frame   = bytarr(film_x,film_y)
ave_arr = fltarr(film_x,film_y)

; compute average image
;skip first 80 frames
;for j = 0, 50 do begin
;    readu, 1, frame
;endfor

;film_l=600
;for j = 500, film_l - 1 do begin
;for j = 0, film_l - 1 do begin
for j = 50, film_l + 50 - 1 do begin
	if((j mod 5) eq 0) then print, j, film_l
	readu, 1, frame
	ave_arr = ave_arr + frame
endfor
close, 1
ave_arr = ave_arr/float(film_l)*1
frame = byte(ave_arr);2*byte(ave_arr)

tv, frame

WRITE_TIFF, dir + run + "_ave.tif", frame, 1, RED = R_ORIG, GREEN = G_ORIG, BLUE = B_ORIG,compression=2

print, "median value :", median(frame)

end