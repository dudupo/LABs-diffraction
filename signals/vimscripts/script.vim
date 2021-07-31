
function! Insert_in_line(line_num)
:execute ":" . a:line_num 
:execute "normal i\<ENTER>\<ESC>"
:execute ":" . a:line_num 
endfunction

function! ET()
	return "\<ENTER>\<TAB>" 
endfunction

function! Run()
	:execute "normal \<ESC>\:w\<ENTER>"
	let output = system("python3 "  . @% )
	echo output 
endfunction

function! Init_py_file()
:execute "normal idef main():" 
		\ . ET() . "pass\<ENTER>if __name__ == \"__main__\":"
		\ . ET() . "main()\<ENTER>"

	:call Insert_in_line(1)

	let import_dicts = {
		\"numpy" : "np",
		\"matplotlib.pyplot" : "plt"
	\} 
	for [key, value] in items( import_dicts )  
		:execute "normal iimport " . key . " as " . value ."\<ENTER> \<ESC>" 
	endfor
endfunction 
