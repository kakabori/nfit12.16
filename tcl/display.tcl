proc open_funclmdif_window {} {
  global flw
  if {![winfo exists .funclmdif]} {
    toplevel .funclmdif
    set flw .funclmdif
    wm title $flw "funclmdif"
    text $flw.text -state disabled -width 80 -height 10 \
                   -yscrollcommand "$flw.scroll set" -cursor arrow
    scrollbar $flw.scroll -command "$flw.text yview"
    #pack $flw.scroll $flw.text -side right -fill y
    pack $flw.scroll -side right -fill y
    pack $flw.text -side left -fill both -expand true
    show_log_file
  }
}


proc show_log_file {} {
  #if {![winfo exists .funclmdif]} {
  #  open_funclmdif_window
  #}
  global flw
  set fileId [open "log.txt" r]
  set data [read $fileId]
  close $fileId
  $flw.text config -state normal
  $flw.text insert end $data
  $flw.text see end
  $flw.text config -state disabled
}


proc insert_buffer_file {} {
  if {![winfo exists .funclmdif]} {
    open_funclmdif_window
  }
  global flw
  set fileId [open "buffer.txt" r]
  set data [read $fileId]
  close $fileId
  $flw.text config -state normal
  $flw.text insert end $data
  $flw.text see end
  $flw.text config -state disabled
  
  # Empty the content of buffer.txt
  set fileId [open "buffer.txt" w]
  close $fileId
}


proc show_funclmdif {msg} {
  global flw
  if {[winfo exists .funclmdif]} {
	  $flw.text insert end $msg
	  $flw.text see end
	} else {
	  puts "flw does not exist!"
	}
}
